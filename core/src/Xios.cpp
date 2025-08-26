/*!
 * @author  Tom Meltzer <tdm39@cam.ac.uk>
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @author  Adeleke Bankole <ab3191@cam.ac.uk>
 * @brief   XIOS interface implementation
 * @details
 *
 * Implementation of XIOS interface
 *
 * This C++ interface is designed to implement core functionality of XIOS so
 * that it can be used in nextSIM-DG.
 *
 * To enable XIOS in nextSIM-DG add the following lines to the config file.
 *   [xios]
 *   enable = true
 *
 * The start time, timestep, and output period will also be read from the
 * following config file entries. (Values shown below are the defaults, while
 * ellipses imply that no default is set.)
 *   [model]
 *   start = 1970-01-01T00:00:00Z
 *   time_step = P0-0T01:00:00
 *   [XiosInput]
 *   period = ...
 *   filename = ...
 *   field_names = ...
 *   [XiosOutput]
 *   period = ...
 *   filename = ...
 *   field_names = ...
 */
#include <boost/date_time/posix_time/time_parsers.hpp>
#if USE_XIOS

#include "StructureModule/include/ParametricGrid.hpp"
#include "include/Finalizer.hpp"
#include "include/ModelMetadata.hpp"
#include "include/ParallelNetcdfFile.hpp"
#include "include/Xios.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/format.hpp>
#include <boost/format/group.hpp>
#include <filesystem>
#include <include/xios_c_interface.hpp>
#include <iostream>
#include <mpi.h>
#include <ncDim.h>
#include <ncException.h>
#include <ncFile.h>
#include <ncGroup.h>
#include <ncVar.h>
#include <regex>
#include <string>

namespace Nextsim {

static const std::string xOutputPfx = "XiosOutput";
static const std::string xInputPfx = "XiosInput";
static const std::map<int, std::string> keyMap = { { Xios::ENABLED_KEY, "xios.enable" },
    { Xios::START_TIME_KEY, "model.start" }, { Xios::TIME_STEP_KEY, "model.time_step" },
    { Xios::OUTPUT_PERIOD_KEY, xOutputPfx + ".period" },
    { Xios::OUTPUT_FILENAME_KEY, xOutputPfx + ".filename" },
    { Xios::OUTPUT_FIELD_NAMES_KEY, xOutputPfx + ".field_names" },
    { Xios::INPUT_PERIOD_KEY, xInputPfx + ".period" },
    { Xios::INPUT_FILENAME_KEY, xInputPfx + ".filename" },
    { Xios::INPUT_FIELD_NAMES_KEY, xInputPfx + ".field_names" } };

//! Enable XIOS in the 'config'
void enableXios()
{
    std::stringstream config;
    config << "[xios]" << std::endl << "enable = true" << std::endl;
    Configurator::addStream(std::unique_ptr<std::istream>(new std::stringstream(config.str())));
}

/*!
 * Constructor: Configure an XIOS server
 *
 * @param calendartype Type of calendar to use
 */
Xios::Xios(const std::string contextid, const std::string calendartype)
{
    static bool firstTime = true;
    contextId = contextid;
    calendarType = calendartype;
    configure();
    static bool doneOnce = doOnce();

    // Create the input and output files (if found in the config)
    if (firstTime) {
        istringstream(Configured::getConfiguration(keyMap.at(INPUT_FILENAME_KEY), std::string()))
            >> inputFilename;
        if (inputFilename.length() > 0) {
            inputFileId = ((std::filesystem::path)inputFilename).replace_extension();
        }
        istringstream(Configured::getConfiguration(keyMap.at(OUTPUT_FILENAME_KEY), std::string()))
            >> outputFilename;
        if (outputFilename.length() > 0) {
            outputFileId = ((std::filesystem::path)outputFilename).replace_extension();
        }

        for (std::string fileId : { inputFileId, outputFileId }) {
            if (fileId.length() > 0) {
                createFile(fileId);

                // Set file name
                xios::CFile* file = getFile(fileId);
                cxios_set_file_name(file, fileId.c_str(), fileId.length());
                if (!cxios_is_defined_file_name(file)) {
                    throw std::runtime_error("Xios: Failed to set name for file '" + fileId + "'");
                }
            }
        }
    }
    firstTime = false;
}

bool Xios::doOnce()
{
    // Register the finalization function here
    Finalizer::registerUnique(finalize);
    return true;
}

//! Destructor
Xios::~Xios() { finalize(); }

//! Close XIOS context definition once xml config has been read and calendar settings updated
void Xios::close_context_definition()
{
    if (isEnabled) {
        cxios_context_close_definition();
    }
}

//! Finalize XIOS context
void Xios::context_finalize()
{
    if (isEnabled) {
        cxios_context_finalize();
    }
}

//! Finalize XIOS server
void Xios::finalize()
{
    if (isEnabled) {
        cxios_finalize();
    }
    isEnabled = false;
}

/*!
 * Overrides `Configure` method from `Configured`
 *
 * Configure the XIOS server if XIOS is enabled in the settings.
 *
 */
void Xios::configure()
{
    // Check if XIOS is enabled in the nextSIM-DG configuration
    istringstream(Configured::getConfiguration(keyMap.at(ENABLED_KEY), std::string()))
        >> std::boolalpha >> isEnabled;

    // Extract the start time from the model configuration
    std::string startTimeStr;
    istringstream(Configured::getConfiguration(keyMap.at(START_TIME_KEY), std::string()))
        >> startTimeStr;
    if (startTimeStr.length() == 0) {
        Logged::warning("Xios: Setting default start: 1970-01-01T00:00:00Z");
        startTimeStr = "1970-01-01T00:00:00Z";
    }
    startTime = TimePoint(startTimeStr);

    // Extract the timestep from the model configuration
    std::string timeStepStr;
    istringstream(Configured::getConfiguration(keyMap.at(TIME_STEP_KEY), std::string()))
        >> timeStepStr;
    if (timeStepStr.length() == 0) {
        Logged::warning("Xios: Setting default time_step: P0-0T01:00:00");
        timeStepStr = "P0-0T01:00:00";
    }
    timestep = Duration(timeStepStr);

    if (isEnabled) {
        configureServer();
    }
}

//! Configure calendar settings
void Xios::configureServer()
{
    // Initialize XIOS Server process and store MPI communicator
    clientId = "client";
    nullComm_F = MPI_Comm_c2f(MPI_COMM_NULL);
    cxios_init_client(clientId.c_str(), clientId.length(), &nullComm_F, &clientComm_F);

    // Initialize MPI rank and size
    clientComm = MPI_Comm_f2c(clientComm_F);
    MPI_Comm_rank(clientComm, &mpi_rank);
    MPI_Comm_size(clientComm, &mpi_size);

    // Initialize 'nextSIM-DG' context
    cxios_context_initialize(contextId.c_str(), contextId.length(), &clientComm_F);

    // Initialize calendar wrapper for 'nextSIM-DG' context
    cxios_get_current_calendar_wrapper(&clientCalendar);
    cxios_set_calendar_wrapper_type(clientCalendar, calendarType.c_str(), calendarType.length());
    cxios_set_calendar_wrapper_timestep(clientCalendar, convertDurationToXios(timestep));
    cxios_create_calendar(clientCalendar);
    cxios_update_calendar_timestep(clientCalendar);

    // Set default calendar origin
    setCalendarOrigin(TimePoint("1970-01-01T00:00:00Z")); // Unix epoch

    // Set start time from configuration file
    setCalendarStart(TimePoint(startTime));
}

/*!
 * @return size of the client MPI communicator
 */
int Xios::getClientMPISize() { return mpi_size; }

/*!
 * @return rank of the client MPI communicator
 */
int Xios::getClientMPIRank() { return mpi_rank; }

/*!
 * Verify XIOS server is initialized
 *
 * @return true when XIOS server is initialized
 */
bool Xios::isInitialized()
{
    bool init = false;
    cxios_context_is_initialized(contextId.c_str(), contextId.length(), &init);
    return init;
}

/*!
 * Return datetime as std::string using ISO 8601 format (default).
 *
 * - If `isoFormat` is true  format will be 2023-03-03T17:11:00Z
 * - If `isoFormat` is false format will be 2023-03-03 17:11:00
 *
 * @param XIOS datetime representation
 * @param isoFormat as bool
 * @return corresponding string representation
 */
std::string Xios::convertXiosDatetimeToString(const cxios_date datetime, const bool isoFormat)
{
    boost::format fmt;
    if (isoFormat) {
        fmt = boost::format("%1$4d-%2$02d-%3$02dT%4$02d:%5$02d:%6$02dZ") % datetime.year
            % datetime.month % datetime.day % datetime.hour % datetime.minute % datetime.second;
    } else {
        fmt = boost::format("%1$4d-%2$02d-%3$02d %4$02d:%5$02d:%6$02d") % datetime.year
            % datetime.month % datetime.day % datetime.hour % datetime.minute % datetime.second;
    }
    return fmt.str();
}

/*!
 * Return std::string in ISO 8601 format (default) as an XIOS datetime object.
 *
 * - If `isoFormat` is true  format will be 2023-03-03T17:11:00Z
 * - If `isoFormat` is false format will be 2023-03-03 17:11:00
 *
 * @param string representation
 * @param isoFormat as bool
 * @return corresponding XIOS datetime representation
 */
cxios_date Xios::convertStringToXiosDatetime(const std::string datetimeStr, const bool isoFormat)
{
    std::string str = datetimeStr;
    if (isoFormat) {
        str = str.replace(10, 1, " "); // replaces T with a space
        str = str.replace(19, 1, " "); // replaces Z with a space
    }
    return cxios_date_convert_from_string(str.c_str(), str.length());
}

/*!
 * Convert a C-string to a C++ `std::string`.
 *
 * @param C-string
 * @param length of C-string
 * @return C++ string version
 */
std::string Xios::convertCStrToCppStr(const char* cStr, int cStrLen)
{
    std::string cppStr(cStr, cStrLen);
    boost::algorithm::trim_right(cppStr);
    return cppStr;
}

/*!
 * Convert an XIOS duration object into a nextSIM-DG one.
 *
 * @param XIOS duration object
 * @return nextSIM-DG version
 */
Duration Xios::convertDurationFromXios(const cxios_duration duration)
{
    char cStr[cStrLen];
    cxios_duration_convert_to_string(duration, cStr, cStrLen);
    std::string durationStr = convertCStrToCppStr(cStr, cStrLen);
    boost::erase_all(durationStr, "s");
    return Duration(std::stod(durationStr));
}

/*!
 * Convert a nextSIM-DG duration object into an XIOS one.
 *
 * @param nextSIM-DG duration object
 * @return XIOS version
 */
cxios_duration Xios::convertDurationToXios(const Duration duration)
{
    return cxios_duration({ 0.0, 0.0, 0.0, 0.0, 0.0, duration.seconds() });
}

/*!
 * Set calendar origin
 *
 * @param origin
 */
void Xios::setCalendarOrigin(const TimePoint origin)
{
    cxios_date datetime = convertStringToXiosDatetime(origin.format(), true);
    cxios_set_calendar_wrapper_date_time_origin(clientCalendar, datetime);
}

/*!
 * Set calendar start date
 *
 * @param start date
 */
void Xios::setCalendarStart(const TimePoint start)
{
    cxios_date datetime = convertStringToXiosDatetime(start.format(), true);
    cxios_set_calendar_wrapper_date_start_date(clientCalendar, datetime);
}

/*!
 * Set calendar timestep
 *
 * @param timestep
 */
void Xios::setCalendarTimestep(const Duration timestep)
{
    cxios_set_calendar_wrapper_timestep(clientCalendar, convertDurationToXios(timestep));
    cxios_update_calendar_timestep(clientCalendar);
}

/*!
 * Update XIOS calendar iteration/step number to some value
 *
 * @param Step number to update to
 */
void Xios::setCalendarStep(const int stepNumber) { cxios_update_calendar(stepNumber); }

/*!
 * Increment XIOS' calendar iteration/step number by one.
 */
void Xios::incrementCalendar() { setCalendarStep(getCalendarStep() + 1); }

/*!
 * Get calendar type
 *
 * @return calendar type
 */
std::string Xios::getCalendarType()
{
    char cStr[cStrLen];
    cxios_get_calendar_wrapper_type(clientCalendar, cStr, cStrLen);
    return convertCStrToCppStr(cStr, cStrLen);
}

/*!
 * Get calendar origin
 *
 * @return calendar origin
 */
TimePoint Xios::getCalendarOrigin()
{
    if (!cxios_is_defined_calendar_wrapper_time_origin(clientCalendar)) {
        throw std::runtime_error("Xios: Calendar origin has not been set");
    }
    cxios_date calendar_origin;
    cxios_get_calendar_wrapper_date_time_origin(clientCalendar, &calendar_origin);
    return TimePoint(convertXiosDatetimeToString(calendar_origin, true));
}

/*!
 * Get calendar start date
 *
 * @return calendar start date
 */
TimePoint Xios::getCalendarStart()
{
    if (!cxios_is_defined_calendar_wrapper_start_date(clientCalendar)) {
        throw std::runtime_error("Xios: Calendar start date has not been set");
    }
    cxios_date calendar_start;
    cxios_get_calendar_wrapper_date_start_date(clientCalendar, &calendar_start);
    return TimePoint(convertXiosDatetimeToString(calendar_start, true));
}

/*!
 * Get calendar timestep
 *
 * @return calendar timestep
 */
Duration Xios::getCalendarTimestep()
{
    if (!cxios_is_defined_calendar_wrapper_timestep(clientCalendar)) {
        throw std::runtime_error("Xios: Calendar timestep has not been set");
    }
    cxios_duration calendar_timestep;
    cxios_get_calendar_wrapper_timestep(clientCalendar, &calendar_timestep);
    return convertDurationFromXios(calendar_timestep);
}

/*!
 * Get calendar step
 *
 * @return calendar step
 */
int Xios::getCalendarStep() { return clientCalendar->getCalendar()->getStep(); }

/*!
 * Get current calendar date
 *
 * @return current calendar date
 */
std::string Xios::getCurrentDate(const bool isoFormat)
{
    cxios_date xiosDate;
    cxios_get_current_date(&xiosDate);
    return convertXiosDatetimeToString(xiosDate, isoFormat);
}

/*!
 * Get the axis_definition group
 *
 * @return a pointer to the XIOS CAxisGroup object
 */
xios::CAxisGroup* Xios::getAxisGroup()
{
    const std::string groupId = "axis_definition";
    xios::CAxisGroup* group = NULL;
    cxios_axisgroup_handle_create(&group, groupId.c_str(), groupId.length());
    if (!group) {
        throw std::runtime_error("Xios: Null pointer for group 'axis_definition'");
    }
    return group;
}

/*!
 * Get the axis associated with a given ID
 *
 * @param the axis ID
 * @return a pointer to the XIOS CAxis object
 */
xios::CAxis* Xios::getAxis(const std::string axisId)
{
    bool exists;
    cxios_axis_valid_id(&exists, axisId.c_str(), axisId.length());
    if (!exists) {
        throw std::runtime_error("Xios: Undefined axis '" + axisId + "'");
    }
    xios::CAxis* axis = NULL;
    cxios_axis_handle_create(&axis, axisId.c_str(), axisId.length());
    if (!axis) {
        throw std::runtime_error("Xios: Null pointer for axis '" + axisId + "'");
    }
    return axis;
}

/*!
 * Create an axis with some ID.
 *
 * @param the axis ID
 */
void Xios::createAxis(const std::string axisId)
{
    bool exists;
    cxios_axis_valid_id(&exists, axisId.c_str(), axisId.length());
    if (exists) {
        throw std::runtime_error("Xios: Axis '" + axisId + "' already exists");
    }
    xios::CAxis* axis = NULL;
    cxios_xml_tree_add_axis(getAxisGroup(), &axis, axisId.c_str(), axisId.length());
    if (!axis) {
        throw std::runtime_error("Xios: Null pointer for axis '" + axisId + "'");
    }
    cxios_axis_valid_id(&exists, axisId.c_str(), axisId.length());
    if (!exists) {
        throw std::runtime_error("Xios: Failed to create axis '" + axisId + "'");
    }
}

/*!
 * Set the size of a given axis (the number of global points)
 *
 * @param the axis ID
 * @param the size to set
 */
void Xios::setAxisSize(const std::string axisId, const size_t size)
{
    xios::CAxis* axis = getAxis(axisId);
    if (cxios_is_defined_axis_n_glo(axis)) {
        Logged::warning("Xios: Size already set for axis '" + axisId + "'");
    }
    cxios_set_axis_n_glo(axis, (int)size);
    if (!cxios_is_defined_axis_n_glo(axis)) {
        throw std::runtime_error("Xios: Failed to set size for axis '" + axisId + "'");
    }
}

/*!
 * Set the values associated with a given axis
 *
 * @param the axis ID
 * @param the values to set
 */
void Xios::setAxisValues(const std::string axisId, std::vector<double> values)
{
    xios::CAxis* axis = getAxis(axisId);
    if (cxios_is_defined_axis_value(axis)) {
        Logged::warning("Xios: Values already set for axis '" + axisId + "'");
    }
    if (!cxios_is_defined_axis_n_glo(axis)) {
        setAxisSize(axisId, values.size());
    }
    int size = getAxisSize(axisId);
    if (size != values.size()) {
        throw std::runtime_error("Xios: Size incompatible with values for axis '" + axisId + "'");
    }
    cxios_set_axis_value(axis, values.data(), &size);
    if (!cxios_is_defined_axis_value(axis)) {
        throw std::runtime_error("Xios: Failed to set values for axis '" + axisId + "'");
    }
}

/*!
 * Get the size of a given axis (the number of global points)
 *
 * @param the axis ID
 * @return size of the corresponding axis
 */
size_t Xios::getAxisSize(const std::string axisId)
{
    xios::CAxis* axis = getAxis(axisId);
    if (!cxios_is_defined_axis_n_glo(axis)) {
        throw std::runtime_error("Xios: Undefined size for axis '" + axisId + "'");
    }
    int size;
    cxios_get_axis_n_glo(axis, &size);
    return (size_t)size;
}

/*!
 * Get the values associated with a given axis
 *
 * @param the axis ID
 * @return the corresponding values
 */
std::vector<double> Xios::getAxisValues(const std::string axisId)
{
    xios::CAxis* axis = getAxis(axisId);
    if (!cxios_is_defined_axis_value(axis)) {
        throw std::runtime_error("Xios: Undefined values for axis '" + axisId + "'");
    }
    int size = getAxisSize(axisId);
    double* values = new double[size];
    cxios_get_axis_value(axis, values, &size);
    std::vector<double> vec(values, values + size);
    delete[] values;
    return vec;
}

/*!
 * Get the domain_definition group
 *
 * @return a pointer to the XIOS CDomainGroup object
 */
xios::CDomainGroup* Xios::getDomainGroup()
{
    const std::string groupId = "domain_definition";
    xios::CDomainGroup* group = NULL;
    cxios_domaingroup_handle_create(&group, groupId.c_str(), groupId.length());
    if (!group) {
        throw std::runtime_error("Xios: Null pointer for group 'domain_definition'");
    }
    return group;
}

/*!
 * Get the domain associated with a given ID
 *
 * @return a pointer to the XIOS CDomain object
 */
xios::CDomain* Xios::getDomain(const std::string domainId)
{
    bool exists;
    cxios_domain_valid_id(&exists, domainId.c_str(), domainId.length());
    if (!exists) {
        throw std::runtime_error("Xios: Undefined domain '" + domainId + "'");
    }
    xios::CDomain* domain = NULL;
    cxios_domain_handle_create(&domain, domainId.c_str(), domainId.length());
    if (!domain) {
        throw std::runtime_error("Xios: Null pointer for domain '" + domainId + "'");
    }
    return domain;
}

/*!
 * Create axes, domains, and grids based off the provided metadata.
 *
 * @param metadata ModelMetadata object containing the partition metadata
 */
void Xios::affixModelMetadata(ModelMetadata& metadata)
{
    // Initial read of the NetCDF file to deduce the dimensions
    if (inputFilename.length() > 0) {
        try {
            netCDF::NcFilePar ncFile(inputFilename, netCDF::NcFile::read, metadata.mpiComm);

            // Dimensions and DG components
            std::multimap<std::string, netCDF::NcDim> dimMap = ncFile.getDims();
            for (auto entry : ModelArray::definedDimensions) {
                auto dimType = entry.first;
                // TODO: Account for DG
                // if (dimCompMap.count(dimType) > 0)
                //     // TODO Assertions that DG in the file equals the compile time DG in the
                //     // model. See #205
                //     continue;

                ModelArray::DimensionSpec& dimensionSpec = entry.second;
                // Find dimensions in the netCDF file by their name in the ModelArray details
                netCDF::NcDim dim = ncFile.getDim(dimensionSpec.name);
                // Also check the old name
                if (dim.isNull()) {
                    dim = ncFile.getDim(dimensionSpec.altName);
                }
                // If we didn't find a dimension with the dimensions name or altName, throw.
                if (dim.isNull()) {
                    throw std::out_of_range(
                        std::string(
                            "No netCDF dimension found corresponding to the dimension named ")
                        + dimensionSpec.name + std::string(" or ") + dimensionSpec.altName);
                }
                auto dimName = dim.getName();
                size_t localLength = 0;
                size_t start = 0;
                if (dimType == ModelArray::Dimension::X) {
                    localLength = metadata.localExtentX;
                    start = metadata.localCornerX;
                } else if (dimType == ModelArray::Dimension::Y) {
                    localLength = metadata.localExtentY;
                    start = metadata.localCornerY;
                } else if (dimType == ModelArray::Dimension::XVERTEX) {
                    localLength = metadata.localExtentX + 1;
                    start = metadata.localCornerX;
                } else if (dimType == ModelArray::Dimension::YVERTEX) {
                    localLength = metadata.localExtentY + 1;
                    start = metadata.localCornerY;
                } else {
                    localLength = dim.getSize();
                    start = 0;
                }
                ModelArray::setDimension(dimType, dim.getSize(), localLength, start);
            }

            // Create map for field types
            const std::map<std::string, ModelArray::Type> dimensionKeys = {
                { "yx", ModelArray::Type::H },
                { "ydimxdim", ModelArray::Type::H },
                { "yxdg_comp", ModelArray::Type::DG },
                { "ydimxdimdg_comp", ModelArray::Type::DG },
                { "yxdgstress_comp", ModelArray::Type::DGSTRESS },
                { "ydimxdimdgstress_comp", ModelArray::Type::DGSTRESS },
                { "ycgxcg", ModelArray::Type::CG },
                { "yvertexxvertexncoords", ModelArray::Type::VERTEX },
            };

            // Determine field types
            for (auto entry : ncFile.getVars()) {
                const std::string& fieldId = entry.first;
                netCDF::NcVar& var = entry.second;
                // Determine the type from the dimensions
                std::vector<netCDF::NcDim> varDims = var.getDims();
                std::string dimKey = "";
                for (netCDF::NcDim& dim : varDims) {
                    const std::string name = dim.getName();
                    // Skip the time_counter dim as it's handled differently
                    if (name != "time_counter") {
                        dimKey += dim.getName();
                    }
                }
                // Skip invalid dimension keys
                if (!dimensionKeys.count(dimKey)) {
                    continue;
                }
                ModelArray::Type type = dimensionKeys.at(dimKey);
                setFieldType(fieldId, type);
            }
            ncFile.close();
        } catch (const netCDF::exceptions::NcException& nce) {
            std::string ncWhat(nce.what());
            ncWhat += ": " + inputFilename;
            throw std::runtime_error(ncWhat);
        }
    }

    // Create XIOS domains
    xios::CDomain* domain = NULL;
    for (ModelArray::Type type : { ModelArray::Type::H, ModelArray::Type::VERTEX }) {
        const std::string domainId = domainIds[type];
        bool exists;
        cxios_domain_valid_id(&exists, domainId.c_str(), domainId.length());
        if (exists) {
            throw std::runtime_error("Xios: Domain '" + domainId + "' already exists");
        }
        cxios_xml_tree_add_domain(getDomainGroup(), &domain, domainId.c_str(), domainId.length());
        if (!domain) {
            throw std::runtime_error("Xios: Null pointer for domain '" + domainId + "'");
        }
        cxios_domain_valid_id(&exists, domainId.c_str(), domainId.length());
        if (!exists) {
            throw std::runtime_error("Xios: Failed to create domain '" + domainId + "'");
        }

        // Set domain type
        const std::string domainType = "rectilinear";
        if (cxios_is_defined_domain_type(domain)) {
            Logged::warning("Xios: Overwriting type for domain '" + domainId + "'");
        }
        cxios_set_domain_type(domain, domainType.c_str(), domainType.length());
        if (!cxios_is_defined_domain_type(domain)) {
            throw std::runtime_error("Xios: Failed to set type for domain '" + domainId + "'");
        }
    }

    // TODO: Automate domain setup

    // Set metadata for 'HDomain'
    const std::string hDomainId = domainIds[ModelArray::Type::H];
    domain = getDomain(hDomainId);
    cxios_set_domain_ni_glo(domain, (int)metadata.globalExtentX);
    cxios_set_domain_nj_glo(domain, (int)metadata.globalExtentY);
    cxios_set_domain_ibegin(domain, (int)metadata.localCornerX);
    cxios_set_domain_jbegin(domain, (int)metadata.localCornerY);
    cxios_set_domain_ni(domain, (int)metadata.localExtentX);
    cxios_set_domain_nj(domain, (int)metadata.localExtentY);
    const std::vector<std::string> hDomainNames = { "x", "y" };
    cxios_set_domain_dim_i_name(domain, hDomainNames[0].c_str(), hDomainNames[0].length());
    cxios_set_domain_dim_j_name(domain, hDomainNames[1].c_str(), hDomainNames[1].length());
    cxios_set_domain_lon_name(domain, hDomainNames[0].c_str(), hDomainNames[0].length());
    cxios_set_domain_lat_name(domain, hDomainNames[1].c_str(), hDomainNames[1].length());

    // Set metadata for 'VertexDomain'
    const std::string vertexDomainId = domainIds[ModelArray::Type::VERTEX];
    domain = getDomain(vertexDomainId);
    cxios_set_domain_ni_glo(domain, (int)metadata.globalExtentX + 1);
    cxios_set_domain_nj_glo(domain, (int)metadata.globalExtentY + 1);
    cxios_set_domain_ibegin(domain, (int)metadata.localCornerX);
    cxios_set_domain_jbegin(domain, (int)metadata.localCornerY);
    cxios_set_domain_ni(domain, (int)metadata.localExtentX + 1);
    cxios_set_domain_nj(domain, (int)metadata.localExtentY + 1);
    const std::vector<std::string> vertexDomainNames = { "xvertex", "yvertex" };
    cxios_set_domain_dim_i_name(
        domain, vertexDomainNames[0].c_str(), vertexDomainNames[0].length());
    cxios_set_domain_dim_j_name(
        domain, vertexDomainNames[1].c_str(), vertexDomainNames[1].length());
    cxios_set_domain_lon_name(domain, vertexDomainNames[0].c_str(), vertexDomainNames[0].length());
    cxios_set_domain_lat_name(domain, vertexDomainNames[1].c_str(), vertexDomainNames[1].length());

    // Check the domains were set up correctly
    for (ModelArray::Type type : { ModelArray::Type::H, ModelArray::Type::VERTEX }) {
        const std::string domainId = domainIds[type];
        domain = getDomain(domainId);
        if (!cxios_is_defined_domain_ni_glo(domain)) {
            throw std::runtime_error(
                "Xios: Failed to set global x-size for domain '" + domainId + "'");
        }
        if (!cxios_is_defined_domain_nj_glo(domain)) {
            throw std::runtime_error(
                "Xios: Failed to set global y-size for domain '" + domainId + "'");
        }
        if (!cxios_is_defined_domain_ibegin(domain)) {
            throw std::runtime_error(
                "Xios: Failed to set local starting x-index for domain '" + domainId + "'");
        }
        if (!cxios_is_defined_domain_jbegin(domain)) {
            throw std::runtime_error(
                "Xios: Failed to set local starting y-index for domain '" + domainId + "'");
        }
        if (!cxios_is_defined_domain_ni(domain)) {
            throw std::runtime_error(
                "Xios: Failed to set local x-size for domain '" + domainId + "'");
        }
        if (!cxios_is_defined_domain_nj(domain)) {
            throw std::runtime_error(
                "Xios: Failed to set local y-size for domain '" + domainId + "'");
        }
        if (!cxios_is_defined_domain_dim_i_name(domain)) {
            throw std::runtime_error(
                "Xios: Failed to set x-coordinate name for domain '" + domainId + "'");
        }
        if (!cxios_is_defined_domain_dim_j_name(domain)) {
            throw std::runtime_error(
                "Xios: Failed to set y-coordinate name for domain '" + domainId + "'");
        }
        if (!cxios_is_defined_domain_lat_name(domain)) {
            throw std::runtime_error(
                "Xios: Failed to set latitude name for domain '" + domainId + "'");
        }
        if (!cxios_is_defined_domain_lon_name(domain)) {
            throw std::runtime_error(
                "Xios: Failed to set longitude name for domain '" + domainId + "'");
        }
    }

    // Create axes
    std::map<ModelArray::Type, ModelArray::Dimension> typeDimMap = {
        { ModelArray::Type::DG, ModelArray::Dimension::DG },
        { ModelArray::Type::VERTEX, ModelArray::Dimension::NCOORDS },
    };
    for (auto entry : typeDimMap) {
        ModelArray::Type type = entry.first;
        ModelArray::Dimension dim = entry.second;
        const std::string axisId = axisIds[type];
        createAxis(axisId);
        setAxisSize(axisId, ModelArray::size(dim));
        xios::CAxis* axis = getAxis(axisId);
        const std::string axisName = axisNames[axisId];
        cxios_set_axis_dim_name(axis, axisName.c_str(), axisName.length());
        if (!cxios_is_defined_axis_dim_name(axis)) {
            throw std::runtime_error("Xios: Failed to set name for axis '" + axisId + "'");
        }
    }

    // Create XIOS grid 'HGrid2D' associated with HDomain
    const std::string hGridId = gridIds[ModelArray::Type::H];
    createGrid(hGridId);
    xios::CGrid* grid = getGrid(hGridId);
    cxios_xml_tree_add_domaintogrid(grid, &domain, hDomainId.c_str(), hDomainId.length());

    // Create XIOS grid 'DGGrid2D' associated with HDomain and DGAxis
    const std::string dgGridId = gridIds[ModelArray::Type::DG];
    createGrid(dgGridId);
    grid = getGrid(dgGridId);
    const std::string dgAxisId = axisIds[ModelArray::Type::DG];
    xios::CAxis* axis = getAxis(dgAxisId);
    cxios_xml_tree_add_axistogrid(grid, &axis, dgAxisId.c_str(), dgAxisId.length());
    cxios_xml_tree_add_domaintogrid(grid, &domain, hDomainId.c_str(), hDomainId.length());

    // Create XIOS grid 'VertexGrid2D' associated with VertexDomain
    const std::string vertexGridId = gridIds[ModelArray::Type::VERTEX];
    const std::string vertexAxisId = axisIds[ModelArray::Type::VERTEX];
    createGrid(vertexGridId);
    grid = getGrid(vertexGridId);
    cxios_xml_tree_add_axistogrid(grid, &axis, vertexAxisId.c_str(), vertexAxisId.length());
    cxios_xml_tree_add_domaintogrid(grid, &domain, vertexDomainId.c_str(), vertexDomainId.length());
}

/*!
 * Get the grid_definition group
 *
 * @return a pointer to the XIOS CGridGroup object
 */
xios::CGridGroup* Xios::getGridGroup()
{
    const std::string groupId = "grid_definition";
    xios::CGridGroup* group = NULL;
    cxios_gridgroup_handle_create(&group, groupId.c_str(), groupId.length());
    if (!group) {
        throw std::runtime_error("Xios: Null pointer for group 'grid_definition'");
    }
    return group;
}

/*!
 * Get the grid associated with a given ID
 *
 * @param the grid ID
 * @return a pointer to the XIOS CGrid object
 */
xios::CGrid* Xios::getGrid(const std::string gridId)
{
    bool exists;
    cxios_grid_valid_id(&exists, gridId.c_str(), gridId.length());
    if (!exists) {
        throw std::runtime_error("Xios: Undefined grid '" + gridId + "'");
    }
    xios::CGrid* grid = NULL;
    cxios_grid_handle_create(&grid, gridId.c_str(), gridId.length());
    if (!grid) {
        throw std::runtime_error("Xios: Null pointer for grid '" + gridId + "'");
    }
    return grid;
}

/*!
 * Create a grid with some ID
 *
 * @param the grid ID
 */
void Xios::createGrid(const std::string gridId)
{
    bool exists;
    cxios_grid_valid_id(&exists, gridId.c_str(), gridId.length());
    if (exists) {
        throw std::runtime_error("Xios: Grid '" + gridId + "' already exists");
    }
    xios::CGrid* grid = NULL;
    cxios_xml_tree_add_grid(getGridGroup(), &grid, gridId.c_str(), gridId.length());
    if (!grid) {
        throw std::runtime_error("Xios: Null pointer for grid '" + gridId + "'");
    }
    cxios_grid_valid_id(&exists, gridId.c_str(), gridId.length());
    if (!exists) {
        throw std::runtime_error("Xios: Failed to create grid '" + gridId + "'");
    }
    cxios_set_grid_name(grid, gridId.c_str(), gridId.length());
    if (!cxios_is_defined_grid_name(grid)) {
        throw std::runtime_error("Xios: Failed to set name for grid '" + gridId + "'");
    }
}

/*!
 * Associate an axis with a grid
 *
 * @param the grid ID
 * @param the axis ID
 */
void Xios::gridAddAxis(const std::string gridId, const std::string axisId)
{
    xios::CAxis* axis = getAxis(axisId);
    cxios_xml_tree_add_axistogrid(getGrid(gridId), &axis, axisId.c_str(), axisId.length());
}

/*!
 * Get all axis IDs associated with a given grid
 *
 * @param the grid ID
 * @return all axis IDs associated with the grid
 */
std::vector<std::string> Xios::getGridAxisIds(const std::string gridId)
{
    return getGrid(gridId)->getAxisList();
}

/*!
 * Get the field_definition group
 *
 * @return a pointer to the XIOS CFieldGroup object
 */
xios::CFieldGroup* Xios::getFieldGroup()
{
    const std::string groupId = "field_definition";
    xios::CFieldGroup* group = NULL;
    cxios_fieldgroup_handle_create(&group, groupId.c_str(), groupId.length());
    if (!group) {
        throw std::runtime_error("Xios: Null pointer for group 'field_definition'");
    }
    return group;
}

/*!
 * Get the field associated with a given ID
 *
 * @param the field ID
 * @return a pointer to the XIOS CField object
 */
xios::CField* Xios::getField(const std::string fieldId)
{
    bool exists;
    cxios_field_valid_id(&exists, fieldId.c_str(), fieldId.length());
    if (!exists) {
        throw std::runtime_error("Xios: Undefined field '" + fieldId + "'");
    }
    xios::CField* field = NULL;
    cxios_field_handle_create(&field, fieldId.c_str(), fieldId.length());
    if (!field) {
        throw std::runtime_error("Xios: Null pointer for field '" + fieldId + "'");
    }
    return field;
}

// Extract the field_names entry from the XiosInput or XiosOutput section of the config
std::set<std::string> Xios::configGetFieldNames(const bool reading)
{
    std::string fieldsStr;
    if (reading) {
        istringstream(Configured::getConfiguration(keyMap.at(INPUT_FIELD_NAMES_KEY), std::string()))
            >> fieldsStr;
    } else {
        istringstream(
            Configured::getConfiguration(keyMap.at(OUTPUT_FIELD_NAMES_KEY), std::string()))
            >> fieldsStr;
    }
    std::set<std::string> fieldNames;
    if (fieldsStr.length() > 0) {
        const char delim = ',';
        std::istringstream iss(fieldsStr);
        std::string item;
        while (std::getline(iss, item, delim)) {
            fieldNames.insert(item);
        }
    }
    return fieldNames;
}

// Check whether a fieldId exists in a string of field names separated by commas, as determined by
// the map key
bool Xios::configCheckField(const std::string fieldId, const bool reading)
{
    std::set<std::string> fieldNames = configGetFieldNames(reading);
    return fieldNames.find(fieldId) != fieldNames.end();
}

/*!
 * Create a field with some ID
 *
 * @param the field ID
 */
void Xios::createField(const std::string fieldId)
{
    // Check if the field already exists
    bool exists;
    cxios_field_valid_id(&exists, fieldId.c_str(), fieldId.length());
    if (exists) {
        throw std::runtime_error("Xios: Field '" + fieldId + "' already exists");
    }

    // Check that the field is in the XiosOutput or XiosInput config
    bool readAccess = configCheckField(fieldId, true);
    bool writeAccess = configCheckField(fieldId, false);
    if (!(readAccess || writeAccess)) {
        throw std::runtime_error("Xios: Field '" + fieldId
            + "' cannot be found in the XiosInput or XiosOutput config sections");
    }

    // Determine whether the field has read access
    if (readAccess && writeAccess) {
        throw std::runtime_error("Xios: Field '" + fieldId
            + "' found in both the XiosInput and XiosOutput config sections");
        // TODO: Refactor to allow a field to be both read and written
    }

    // Attempt to create the field
    xios::CField* field = NULL;
    cxios_xml_tree_add_field(getFieldGroup(), &field, fieldId.c_str(), fieldId.length());
    if (!field) {
        throw std::runtime_error("Xios: Null pointer for field '" + fieldId + "'");
    }
    cxios_field_valid_id(&exists, fieldId.c_str(), fieldId.length());
    if (!exists) {
        throw std::runtime_error("Xios: Failed to create field '" + fieldId + "'");
    }
}

/*!
 * Set the operation for a field with a given ID
 *
 * @param the field ID
 * @param operation to set
 */
void Xios::setFieldOperation(const std::string fieldId, const std::string operation)
{
    xios::CField* field = getField(fieldId);
    if (cxios_is_defined_field_operation(field)) {
        Logged::warning("Xios: Overwriting operation for field '" + fieldId + "'");
    }
    cxios_set_field_operation(field, operation.c_str(), operation.length());
    if (!cxios_is_defined_field_operation(field)) {
        throw std::runtime_error("Xios: Failed to set operation for field '" + fieldId + "'");
    }
}

/*!
 * Set the grid reference for a field with a given ID
 *
 * @param the field ID
 * @param grid reference to set
 */
void Xios::setFieldGridRef(const std::string fieldId, const std::string gridRef)
{
    xios::CField* field = getField(fieldId);
    if (cxios_is_defined_field_grid_ref(field)) {
        Logged::warning("Xios: Overwriting grid reference for field '" + fieldId + "'");
    }
    cxios_set_field_grid_ref(field, gridRef.c_str(), gridRef.length());
    if (!cxios_is_defined_field_grid_ref(field)) {
        throw std::runtime_error("Xios: Failed to set grid reference for field '" + fieldId + "'");
    }
}

/*!
 * Set the read access for a field with a given ID
 *
 * @param the field ID
 * @param read access to set
 */
void Xios::setFieldReadAccess(const std::string fieldId, const bool readAccess)
{
    xios::CField* field = getField(fieldId);
    if (cxios_is_defined_field_read_access(field)) {
        Logged::warning("Xios: Overwriting read access for field '" + fieldId + "'");
    }
    cxios_set_field_read_access(field, readAccess);
    if (!cxios_is_defined_field_read_access(field)) {
        throw std::runtime_error("Xios: Failed to set read access for field '" + fieldId + "'");
    }
}

/*!
 * Set the frequency offset for a field with a given ID
 *
 * @param the field ID
 * @param frequency offset to set
 */
void Xios::setFieldFreqOffset(const std::string fieldId, const Duration freqOffset)
{
    xios::CField* field = getField(fieldId);
    if (cxios_is_defined_field_freq_offset(field)) {
        Logged::warning("Xios: Overwriting frequency offset for field '" + fieldId + "'");
    }
    cxios_set_field_freq_offset(field, convertDurationToXios(freqOffset));
    if (!cxios_is_defined_field_freq_offset(field)) {
        throw std::runtime_error(
            "Xios: Failed to set frequency offset for field '" + fieldId + "'");
    }
}

/*!
 * Get the operation associated with a field with a given ID
 *
 * @param the field ID
 * @return operation used for the corresponding field
 */
std::string Xios::getFieldOperation(const std::string fieldId)
{
    xios::CField* field = getField(fieldId);
    if (!cxios_is_defined_field_operation(field)) {
        throw std::runtime_error("Xios: Undefined operation for field '" + fieldId + "'");
    }
    char cStr[cStrLen];
    cxios_get_field_operation(field, cStr, cStrLen);
    return convertCStrToCppStr(cStr, cStrLen);
}

/*!
 * Get the grid reference associated with a field with a given ID
 *
 * @param the field ID
 * @return grid reference used for the corresponding field
 */
std::string Xios::getFieldGridRef(const std::string fieldId)
{
    xios::CField* field = getField(fieldId);
    if (!cxios_is_defined_field_grid_ref(field)) {
        throw std::runtime_error("Xios: Undefined grid reference for field '" + fieldId + "'");
    }
    char cStr[cStrLen];
    cxios_get_field_grid_ref(field, cStr, cStrLen);
    return convertCStrToCppStr(cStr, cStrLen);
}

/*!
 * Get the read access associated with a field with a given ID
 *
 * @param the field ID
 * @return read access used for the corresponding field
 */
bool Xios::getFieldReadAccess(const std::string fieldId)
{
    xios::CField* field = getField(fieldId);
    if (!cxios_is_defined_field_read_access(field)) {
        throw std::runtime_error("Xios: Undefined read access for field '" + fieldId + "'");
    }
    bool readAccess;
    cxios_get_field_read_access(field, &readAccess);
    return readAccess;
}

/*!
 * Get the frequency offset associated with a field with a given ID
 *
 * @param the field ID
 * @return frequency offset used for the corresponding field
 */
Duration Xios::getFieldFreqOffset(const std::string fieldId)
{
    xios::CField* field = getField(fieldId);
    if (!cxios_is_defined_field_freq_offset(field)) {
        throw std::runtime_error("Xios: Undefined frequency offset for field '" + fieldId + "'");
    }
    cxios_duration duration;
    cxios_get_field_freq_offset(field, &duration);
    char cStr[cStrLen];
    cxios_duration_convert_to_string(duration, cStr, cStrLen);
    return convertDurationFromXios(duration);
}

/*!
 * Get the field type associated with a field with a given ID
 *
 * @param the field ID
 * @return ModelArray::Type used for the corresponding field
 */
ModelArray::Type Xios::getFieldType(const std::string fieldId) { return fieldTypes[fieldId]; }

/*!
 * Set the field type associated with a field with a given ID
 *
 * @param the field ID
 * @param ModelArray::Type used for the corresponding field
 */
void Xios::setFieldType(const std::string fieldId, ModelArray::Type type)
{
    fieldTypes[fieldId] = type;
    setFieldGridRef(fieldId, gridIds[type]);
}

/*!
 * Get the file_definition group
 *
 * @return a pointer to the XIOS CFileGroup object
 */
xios::CFileGroup* Xios::getFileGroup()
{
    const std::string groupId = "file_definition";
    xios::CFileGroup* group = NULL;
    cxios_filegroup_handle_create(&group, groupId.c_str(), groupId.length());
    if (!group) {
        throw std::runtime_error("Xios: Null pointer for group 'file_definition'");
    }
    return group;
}

/*!
 * Get the file associated with a given ID
 *
 * @param the file ID
 * @return a pointer to the XIOS CFile object
 */
xios::CFile* Xios::getFile(const std::string fileId)
{
    bool exists;
    cxios_file_valid_id(&exists, fileId.c_str(), fileId.length());
    if (!exists) {
        throw std::runtime_error("Xios: Undefined file '" + fileId + "'");
    }
    xios::CFile* file = NULL;
    cxios_file_handle_create(&file, fileId.c_str(), fileId.length());
    if (!file) {
        throw std::runtime_error("Xios: Null pointer for file '" + fileId + "'");
    }
    return file;
}

/*!
 * Create a file with some ID
 *
 * @param the file ID
 */
void Xios::createFile(const std::string fileId)
{
    xios::CFile* file = NULL;
    bool exists;
    cxios_file_valid_id(&exists, fileId.c_str(), fileId.length());
    if (exists) {
        throw std::runtime_error("Xios: File '" + fileId + "' already exists");
    }
    cxios_xml_tree_add_file(getFileGroup(), &file, fileId.c_str(), fileId.length());
    if (!file) {
        throw std::runtime_error("Xios: Null pointer for file '" + fileId + "'");
    }
    cxios_file_valid_id(&exists, fileId.c_str(), fileId.length());
    if (!exists) {
        throw std::runtime_error("Xios: Failed to create file '" + fileId + "'");
    }

    // Determine whether the file is configured for reading or writing
    bool readAccess = ((inputFileId.length() > 0) && (inputFileId == fileId));
    bool writeAccess = ((outputFileId.length() > 0) && (outputFileId == fileId));

    // Check that the filename is not in both the XiosOutput and XiosInput config sections
    if (readAccess && writeAccess) {
        throw std::runtime_error("Xios: File '" + fileId
            + "' found in both the XiosInput and XiosOutput config sections");
        // TODO: Refactor to allow a field to be both read and written
    }

    // Terminate early for special unit test cases, for which IDs start with 'unittest'
    if (fileId.rfind("unittest", 0) == 0) {
        Logged::warning("Xios: Special 'unittest' ID found; skipping automated setup. Are you sure "
                        "you want to do this?");
        return;
    }

    // Check that the filename is in the XiosOutput or XiosInput config section
    if (!(readAccess || writeAccess)) {
        throw std::runtime_error("Xios: File '" + fileId
            + "' cannot be found in the XiosInput or XiosOutput config sections");
    }

    // Set the file mode and some defaults
    if (readAccess) {
        setFileMode(fileId, "read");
    } else {
        setFileMode(fileId, "write");
    }
    setFileType(fileId, "one_file");
    setFileParAccess(fileId, "collective");

    // Set the input or output period based on the model configuration
    std::string periodStr;
    if (readAccess) {
        istringstream(Configured::getConfiguration(keyMap.at(INPUT_PERIOD_KEY), std::string()))
            >> periodStr;
    } else {
        istringstream(Configured::getConfiguration(keyMap.at(OUTPUT_PERIOD_KEY), std::string()))
            >> periodStr;
    }
    if (periodStr.length() > 0) {
        setFileOutputFreq(fileId, Duration(periodStr));
    }

    // Create all fields found in the config based off the field names found in the
    // XiosInput.field_names or XiosOutput.field_names entries in the config.
    for (std::string fieldId : configGetFieldNames(readAccess)) {
        createField(fieldId);
        fileAddField(fileId, fieldId);
        setFieldReadAccess(fieldId, readAccess);

        // Set field name
        xios::CField* field = getField(fieldId);
        cxios_set_field_name(field, fieldId.c_str(), fieldId.length());
        if (!cxios_is_defined_field_name(field)) {
            throw std::runtime_error("Xios: Failed to set name for field '" + fieldId + "'");
        }
    }
}

/*!
 * Set the type of a file with a given ID
 *
 * @param the file ID
 * @param file type to set
 */
void Xios::setFileType(const std::string fileId, const std::string fileType)
{
    xios::CFile* file = getFile(fileId);
    if (cxios_is_defined_file_type(file)) {
        Logged::warning("Xios: Overwriting type for file '" + fileId + "'");
    }
    cxios_set_file_type(file, fileType.c_str(), fileType.length());
    if (!cxios_is_defined_file_type(file)) {
        throw std::runtime_error("Xios: Failed to set type for file '" + fileId + "'");
    }
}

/*!
 * Set the output frequency of a file with a given ID
 *
 * @param the file ID
 * @param output frequency to set
 */
void Xios::setFileOutputFreq(const std::string fileId, const Duration freq)
{
    xios::CFile* file = getFile(fileId);
    if (cxios_is_defined_file_output_freq(file)) {
        Logged::warning("Xios: Overwriting output frequency for file '" + fileId + "'");
    }
    cxios_set_file_output_freq(file, convertDurationToXios(freq));
    if (!cxios_is_defined_file_output_freq(file)) {
        throw std::runtime_error("Xios: Failed to set output frequency for file '" + fileId + "'");
    }
}

/*!
 * Set the split frequency of a file with a given ID
 *
 * @param the file ID
 * @param split frequency to set
 */
void Xios::setFileSplitFreq(const std::string fileId, const Duration freq)
{
    xios::CFile* file = getFile(fileId);
    if (cxios_is_defined_file_split_freq(file)) {
        Logged::warning("Xios: Split frequency already set for file '" + fileId + "'");
    }
    cxios_set_file_split_freq(file, convertDurationToXios(freq));
    if (!cxios_is_defined_file_split_freq(file)) {
        throw std::runtime_error("Xios: Failed to set split frequency for file '" + fileId + "'");
    }
}

/*!
 * Set the mode of a file with a given ID
 *
 * @param the file ID
 * @param file mode to set
 */
void Xios::setFileMode(const std::string fileId, const std::string mode)
{
    xios::CFile* file = getFile(fileId);
    if (cxios_is_defined_file_mode(file)) {
        Logged::warning("Xios: Overwriting mode for file '" + fileId + "'");
    }
    cxios_set_file_mode(file, mode.c_str(), mode.length());
    if (!cxios_is_defined_file_mode(file)) {
        throw std::runtime_error("Xios: Failed to set mode for file '" + fileId + "'");
    }
}

/*!
 * Set the parallel access mode of a file with a given ID
 *
 * @param the file ID
 * @param parallel access mode to set
 */
void Xios::setFileParAccess(const std::string fileId, const std::string parAccess)
{
    xios::CFile* file = getFile(fileId);
    if (cxios_is_defined_file_par_access(file)) {
        Logged::warning("Xios: Overwriting parallel access for file '" + fileId + "'");
    }
    cxios_set_file_par_access(file, parAccess.c_str(), parAccess.length());
    if (!cxios_is_defined_file_par_access(file)) {
        throw std::runtime_error("Xios: Failed to set parallel access for file '" + fileId + "'");
    }
}

/*!
 * Get the type of a file with a given ID
 *
 * @param the file ID
 * @return type of the corresponding file
 */
std::string Xios::getFileType(const std::string fileId)
{
    xios::CFile* file = getFile(fileId);
    if (!cxios_is_defined_file_type(file)) {
        throw std::runtime_error("Xios: Undefined type for file '" + fileId + "'");
    }
    char cStr[cStrLen];
    cxios_get_file_type(file, cStr, cStrLen);
    return convertCStrToCppStr(cStr, cStrLen);
}

/*!
 * Get the output frequency of a file with a given ID
 *
 * @param the file ID
 * @return the corresponding output frequency
 */
Duration Xios::getFileOutputFreq(const std::string fileId)
{
    xios::CFile* file = getFile(fileId);
    if (!cxios_is_defined_file_output_freq(file)) {
        throw std::runtime_error("Xios: Undefined output frequency for file '" + fileId + "'");
    }
    cxios_duration duration;
    cxios_get_file_output_freq(file, &duration);
    return convertDurationFromXios(duration);
}

/*!
 * Get the split frequency of a file with a given ID
 *
 * @param the file ID
 * @return split frequency of the corresponding file
 */
Duration Xios::getFileSplitFreq(const std::string fileId)
{
    xios::CFile* file = getFile(fileId);
    if (!cxios_is_defined_file_split_freq(file)) {
        throw std::runtime_error("Xios: Undefined split frequency for file '" + fileId + "'");
    }
    cxios_duration duration;
    cxios_get_file_split_freq(file, &duration);
    return convertDurationFromXios(duration);
}

/*!
 * Get the mode of a file with a given ID
 *
 * @param the file ID
 * @return mode of the corresponding file
 */
std::string Xios::getFileMode(const std::string fileId)
{
    xios::CFile* file = getFile(fileId);
    if (!cxios_is_defined_file_mode(file)) {
        throw std::runtime_error("Xios: Undefined mode for file '" + fileId + "'");
    }
    char cStr[cStrLen];
    cxios_get_file_mode(file, cStr, cStrLen);
    std::string mode(cStr, cStrLen);
    boost::algorithm::trim_right(mode);
    return mode;
}

/*!
 * Get the parallel access mode of a file with a given ID
 *
 * @param the file ID
 * @return parallel access mode of the corresponding file
 */
std::string Xios::getFileParAccess(const std::string fileId)
{
    xios::CFile* file = getFile(fileId);
    if (!cxios_is_defined_file_par_access(file)) {
        throw std::runtime_error("Xios: Undefined parallel access for file '" + fileId + "'");
    }
    char cStr[cStrLen];
    cxios_get_file_par_access(file, cStr, cStrLen);
    std::string parAccess(cStr, cStrLen);
    boost::algorithm::trim_right(parAccess);
    return parAccess;
}

/*!
 * Get all field IDs associated with a given file
 *
 * @param the file ID
 * @return all field IDs associated with the file
 */
std::vector<std::string> Xios::fileGetFieldIds(const std::string fileId)
{
    std::vector<xios::CField*> fields = getFile(fileId)->getAllFields();
    std::vector<std::string> fieldIds(fields.size());
    for (int i = 0; i < fields.size(); i++) {
        fieldIds[i] = fields[i]->getId();
    }
    return fieldIds;
}

/*!
 * Associate a field with a file
 *
 * @param the file ID
 * @param the field ID
 */
void Xios::fileAddField(const std::string fileId, const std::string fieldId)
{
    xios::CField* field = getField(fieldId);
    cxios_xml_tree_add_fieldtofile(getFile(fileId), &field, fieldId.c_str(), fieldId.length());
}

/*!
 * Send a field to the XIOS server to be written to file
 *
 * @param field name
 * @param reference to the ModelArray containing the data to be written
 */
void Xios::write(const std::string fieldId, ModelArray& modelarray)
{
    const bool readAccess = false;
    std::set<std::string> fieldNames = configGetFieldNames(readAccess);
    if (fieldNames.find(fieldId) == fieldNames.end()) {
        throw std::runtime_error(
            "Xios::write: field " + fieldId + " has not been configured for writing with XIOS.");
    }
    if (modelarray.nDimensions() != 2) {
        throw std::invalid_argument("Only ModelArrays of dimension 2 are supported");
    }
    auto dims = modelarray.dimensions();
    auto type = modelarray.getType();
    if (type == ModelArray::Type::H) {
        cxios_write_data_k82(
            fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0], dims[1], -1);
    } else if (type == ModelArray::Type::VERTEX) {
        cxios_write_data_k83(fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0],
            dims[1], ModelArray::size(ModelArray::Dimension::NCOORDS), -1);
    } else if (type == ModelArray::Type::DG) {
        cxios_write_data_k83(fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0],
            dims[1], ModelArray::size(ModelArray::Dimension::DG), -1);
    } else {
        throw std::invalid_argument("Only HFields, VertexFields, and DGFields are supported");
    }
}

/*!
 * Receive field from XIOS server that has been read from file.
 *
 * @param field name
 * @param reference to the ModelArray containing the data to be written
 */
void Xios::read(const std::string fieldId, ModelArray& modelarray)
{
    const bool readAccess = true;
    std::set<std::string> fieldNames = configGetFieldNames(readAccess);
    if (fieldNames.find(fieldId) == fieldNames.end()) {
        throw std::runtime_error(
            "Xios::read: field " + fieldId + " has not been configured for reading with XIOS.");
    }
    if (modelarray.nDimensions() != 2) {
        throw std::invalid_argument("Only ModelArrays of dimension 2 are supported");
    }
    auto dims = modelarray.dimensions();
    auto type = modelarray.getType();
    if (type == ModelArray::Type::H) {
        cxios_read_data_k82(
            fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0], dims[1]);
    } else if (type == ModelArray::Type::VERTEX) {
        cxios_read_data_k83(fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0],
            dims[1], ModelArray::size(ModelArray::Dimension::NCOORDS));
    } else if (type == ModelArray::Type::DG) {
        cxios_read_data_k83(fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0],
            dims[1], ModelArray::size(ModelArray::Dimension::DG));
    } else {
        throw std::invalid_argument("Only HFields, VertexFields, and DGFields are supported");
    }
}
}

#endif
