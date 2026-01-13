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
 * that it can be used in nextSIM-DG. Run
 * ```sh
 * ./nextsimdg --help-config
 * ```
 * to see the configuration options available for XIOS.
 */
#include <boost/date_time/posix_time/time_parsers.hpp>
#if USE_XIOS

#include "StructureModule/include/ParametricGrid.hpp"
#include "include/Finalizer.hpp"
#include "include/ModelMPI.hpp"
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

#ifndef CGDEGREE
#define CGDEGREE 2 // Define to prevent errors from static analysis tools
#error "CG degree (CGDEGREE) not defined" // But throw an error anyway
#endif

namespace Nextsim {

static const std::string xOutputPfx = "XiosOutput";
static const std::string xInputPfx = "XiosInput";
static const std::string xDiagnosticPfx = "XiosDiagnostic";
static const std::string xForcingPfx = "XiosForcing";
static const std::map<int, std::string> keyMap = { { Xios::ENABLED_KEY, "xios.enable" },
    { Xios::OUTPUT_FIELD_NAMES_KEY, xOutputPfx + ".field_names" },
    { Xios::INPUT_FIELD_NAMES_KEY, xInputPfx + ".field_names" },
    { Xios::DIAGNOSTIC_PERIOD_KEY, xDiagnosticPfx + ".period" },
    { Xios::DIAGNOSTIC_FILE_KEY, xDiagnosticPfx + ".filename" },
    { Xios::DIAGNOSTIC_FIELD_NAMES_KEY, xDiagnosticPfx + ".field_names" },
    { Xios::FORCING_PERIOD_KEY, xForcingPfx + ".period" },
    { Xios::FORCING_FILE_KEY, xForcingPfx + ".filename" },
    { Xios::FORCING_FIELD_NAMES_KEY, xForcingPfx + ".field_names" } };

Xios::HelpMap& Xios::getHelpText(HelpMap& map, bool getAll)
{
    map["Xios"] = {
        { keyMap.at(ENABLED_KEY), ConfigType::BOOLEAN, {}, "", "",
            "Boolean option to toggle whether XIOS is enabled in the build. This should not need "
            "to be modifed by the user. Build nextSIM-DG with XIOS support with the CMake argument "
            "-DENABLE_XIOS=ON, passing the path to your XIOS installation with "
            "-Dxios_DIR=/path/to/xios." },
    };
    map["XiosOutput"] = {
        { keyMap.at(OUTPUT_FIELD_NAMES_KEY), ConfigType::STRING, {}, "restart%Y-%m-%dT%H:%M:%SZ.nc",
            "", "Comma-separated list of field names to be written to the output file." },
    };
    map["XiosInput"] = {
        { keyMap.at(INPUT_FIELD_NAMES_KEY), ConfigType::STRING, {}, "", "",
            "Comma-separated list of field names to be read from the input file." },
    };
    map["XiosDiagnostic"] = {
        { keyMap.at(DIAGNOSTIC_PERIOD_KEY), ConfigType::STRING, {}, "0", "",
            "The period between diagnostics file outputs expected in a file to be "
            "read, formatted as an ISO8601 duration (P prefix) or number of "
            "seconds. A value of zero assumes no intermediate diagnostics files." },
        { keyMap.at(DIAGNOSTIC_FILE_KEY), ConfigType::STRING, {}, "diagnostic%Y-%m-%dT%H:%M:%SZ.nc",
            "", "The file name to be used for diagnostics." },
        { keyMap.at(DIAGNOSTIC_FIELD_NAMES_KEY), ConfigType::STRING, {}, "", "",
            "Comma-separated list of field names to be read from the diagnostics "
            "file." },
    };
    map["XiosForcing"] = {
        { keyMap.at(FORCING_PERIOD_KEY), ConfigType::STRING, {}, "0", "",
            "The period between forcing file outputs expected in a file to be "
            "read, formatted as an ISO8601 duration (P prefix) or number of "
            "seconds. A value of zero assumes no intermediate forcing files." },
        { keyMap.at(FORCING_FILE_KEY), ConfigType::STRING, {}, "", "",
            "The file name to be used for forcings." },
        { keyMap.at(FORCING_FIELD_NAMES_KEY), ConfigType::STRING, {}, "", "",
            "Comma-separated list of field names to be read from the forcings "
            "file." },
    };

    return map;
}

Xios::HelpMap& Xios::getHelpRecursive(HelpMap& map, bool getAll)
{
    getHelpText(map, getAll);
    return map;
}

//! Constructor for the XIOS handler
Xios::Xios()
{
    // Check if XIOS is enabled in the nextSIM-DG configuration
    istringstream(Configured::getConfiguration(keyMap.at(ENABLED_KEY), std::string()))
        >> std::boolalpha >> isEnabled;

    if (isEnabled) {
        configure();
    }
    static bool doneOnce = doOnce();
}

//! Configure XIOS client
void Xios::setupClient()
{
    // Initialize XIOS Server process and store MPI communicator
    nullComm_F = MPI_Comm_c2f(MPI_COMM_NULL);
    cxios_init_client(clientId.c_str(), clientId.length(), &nullComm_F, &clientComm_F);

    // Initialize MPI rank and size
    clientComm = MPI_Comm_f2c(clientComm_F);
    MPI_Comm_rank(clientComm, &mpi_rank);
    MPI_Comm_size(clientComm, &mpi_size);
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
        postprocessOutputFiles();
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
    if (isEnabled) {
        setupClient();
        setupContext();
        setupCalendar();
        setupFiles();
        setupFields();
    }
}

//! Initialize the XIOS context with ID contextId
void Xios::setupContext()
{
    // Initialize the XIOS context 'nextSIM-DG'
    cxios_context_initialize(contextId.c_str(), contextId.length(), &clientComm_F);

    // Verify the XIOS context was created properly
    bool exists;
    cxios_context_valid_id(&exists, contextId.c_str(), contextId.length());
    if (!exists) {
        throw std::runtime_error("Xios: context '" + contextId + "' was not created");
    }

    // Verify the XIOS context has been initialized properly
    bool init;
    cxios_context_is_initialized(contextId.c_str(), contextId.length(), &init);
    if (!init) {
        throw std::runtime_error("Xios: context '" + contextId + "' not initialized");
    }

    // Verify the correct context ID is being used
    xios::CContext* context = NULL;
    cxios_context_get_current(&context);
    char cStr[cStrLen];
    cxios_context_get_id(context, cStr, cStrLen);
    if (convertCStrToCppStr(cStr, cStrLen) != contextId) {
        throw std::runtime_error(
            "Xios: current context ID does not match expected ID '" + contextId + "'");
    }
}

//! Initialize calendar wrapper for the context
// NOTE: The calendar itself is set up in iodef.xml
void Xios::setupCalendar()
{
    cxios_get_current_calendar_wrapper(&clientCalendar);

    // Set timestep from configuration file
    ModelMetadata& metadata = ModelMetadata::getInstance();
    cxios_set_calendar_wrapper_timestep(
        clientCalendar, convertDurationToXios(metadata.stepLength()));
    cxios_update_calendar_timestep(clientCalendar);

    // Verify the timestep was set correctly
    if (!cxios_is_defined_calendar_wrapper_timestep(clientCalendar)) {
        throw std::runtime_error("Xios: Calendar timestep has not been set");
    }

    // Set start time from configuration file
    setCalendarStart(metadata.startTime());
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
std::string Xios::convertXiosDatetimeToString(const cxios_date& datetime, const bool isoFormat)
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
cxios_date Xios::convertStringToXiosDatetime(std::string datetimeStr, const bool isoFormat)
{
    if (isoFormat) {
        datetimeStr = datetimeStr.replace(10, 1, " "); // replaces T with a space
        datetimeStr = datetimeStr.replace(19, 1, " "); // replaces Z with a space
    }
    return cxios_date_convert_from_string(datetimeStr.c_str(), datetimeStr.length());
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
Duration Xios::convertDurationFromXios(const cxios_duration& duration)
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
cxios_duration Xios::convertDurationToXios(const Duration& duration)
{
    return cxios_duration({ 0.0, 0.0, 0.0, 0.0, 0.0, duration.seconds() });
}

/*!
 * Set calendar start date
 *
 * @param start date
 */
void Xios::setCalendarStart(const TimePoint& start)
{
    cxios_date datetime = convertStringToXiosDatetime(start.format(), true);
    cxios_set_calendar_wrapper_date_start_date(clientCalendar, datetime);
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
TimePoint Xios::getCurrentDate()
{
    cxios_date xiosDate;
    cxios_get_current_date(&xiosDate);
    return TimePoint(convertXiosDatetimeToString(xiosDate, true));
}

/*!
 * Get the axis associated with a given ID
 *
 * @param the axis ID
 * @return a pointer to the XIOS CAxis object
 */
xios::CAxis* Xios::getAxis(const std::string& axisId)
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
 * Get the size of a given axis (the number of global points)
 *
 * @param the axis ID
 * @return size of the corresponding axis
 */
size_t Xios::getAxisSize(const std::string& axisId)
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
xios::CDomain* Xios::getDomain(const std::string& domainId)
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
 * @brief   Create XIOS domains associated with each ModelArray type
 *
 * @details This function sets up the XIOS domains for each field type based on the
 *          configuration in the domainIds map and in the ModelMetadata class.
 */
void Xios::setupDomains()
{
    auto& metadata = ModelMetadata::getInstance();

    ModelArray::setNComponents(ModelArray::Type::VERTEX, ModelArray::nCoords);
    ModelArray::setNComponents(ModelArray::Type::DG, getAxisSize("DGAxis"));
    ModelArray::setNComponents(ModelArray::Type::DGSTRESS, getAxisSize("DGSAxis"));
    for (const auto& [type, domainId] : domainIds) {
        bool exists;
        cxios_domain_valid_id(&exists, domainId.c_str(), domainId.length());
        if (exists) {
            continue;
        }

        // Create the domain
        xios::CDomain* domain = NULL;
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

        // Set domain extents based on model metadata
        size_t counter = 0;
        for (ModelArray::Dimension& dim : ModelArray::typeDimensions[type]) {
            if (counter == 0) {
                const std::string dimName = "x";
                int ni_glo;
                int ni;
                int ibegin;
                if (dim == ModelArray::Dimension::X) {
                    ni_glo = metadata.getGlobalExtentX();
                    ni = metadata.getLocalExtentX();
                    ibegin = metadata.getLocalCornerX();
                } else if (dim == ModelArray::Dimension::XVERTEX) {
                    ni_glo = metadata.getGlobalExtentX() + 1;
                    ni = metadata.getLocalExtentX() + 1;
                    ibegin = metadata.getLocalCornerX();
                } else if (dim == ModelArray::Dimension::XCG) {
                    ni_glo = CGDEGREE * metadata.getGlobalExtentX() + 1;
                    ni = CGDEGREE * metadata.getLocalExtentX() + 1;
                    ibegin = CGDEGREE * metadata.getLocalCornerX();
                } else {
                    throw std::runtime_error(
                        "Xios: Could not set domain extents based on dimension '" + dimName + "'");
                }
                cxios_set_domain_ni_glo(domain, ni_glo);
                if (!cxios_is_defined_domain_ni_glo(domain)) {
                    throw std::runtime_error(
                        "Xios: Failed to set global x-size for domain '" + domainId + "'");
                }
                cxios_set_domain_ni(domain, ni);
                if (!cxios_is_defined_domain_ni(domain)) {
                    throw std::runtime_error(
                        "Xios: Failed to set local x-size for domain '" + domainId + "'");
                }
                cxios_set_domain_ibegin(domain, ibegin);
                if (!cxios_is_defined_domain_ibegin(domain)) {
                    throw std::runtime_error(
                        "Xios: Failed to set local starting x-index for domain '" + domainId + "'");
                }
                std::vector<double> lonvalue;
                for (int i = 0; i < ni; i++) {
                    lonvalue.push_back(ibegin + i);
                }
                cxios_set_domain_lonvalue_1d(domain, lonvalue.data(), &ni);
                if (!cxios_is_defined_domain_lonvalue_1d(domain)) {
                    throw std::runtime_error(
                        "Xios: Failed to set local x-indices for domain '" + domainId + "'");
                }
                cxios_set_domain_dim_i_name(domain, dimName.c_str(), dimName.length());
                if (!cxios_is_defined_domain_dim_i_name(domain)) {
                    throw std::runtime_error(
                        "Xios: Failed to set x-coordinate name for domain '" + domainId + "'");
                }
                cxios_set_domain_lon_name(domain, dimName.c_str(), dimName.length());
                if (!cxios_is_defined_domain_lon_name(domain)) {
                    throw std::runtime_error(
                        "Xios: Failed to set longitude name for domain '" + domainId + "'");
                }
            } else if (counter == 1) {
                const std::string dimName = "y";
                int nj_glo;
                int nj;
                int jbegin;
                if (dim == ModelArray::Dimension::Y) {
                    nj_glo = metadata.getGlobalExtentY();
                    nj = metadata.getLocalExtentY();
                    jbegin = metadata.getLocalCornerY();
                } else if (dim == ModelArray::Dimension::YVERTEX) {
                    nj_glo = metadata.getGlobalExtentY() + 1;
                    nj = metadata.getLocalExtentY() + 1;
                    jbegin = metadata.getLocalCornerY();
                } else if (dim == ModelArray::Dimension::YCG) {
                    nj_glo = CGDEGREE * metadata.getGlobalExtentY() + 1;
                    nj = CGDEGREE * metadata.getLocalExtentY() + 1;
                    jbegin = CGDEGREE * metadata.getLocalCornerY();
                } else {
                    throw std::runtime_error(
                        "Xios: Could not set domain extents based on dimension '" + dimName + "'");
                }
                cxios_set_domain_nj_glo(domain, nj_glo);
                if (!cxios_is_defined_domain_nj_glo(domain)) {
                    throw std::runtime_error(
                        "Xios: Failed to set global y-size for domain '" + domainId + "'");
                }
                cxios_set_domain_nj(domain, nj);
                if (!cxios_is_defined_domain_nj(domain)) {
                    throw std::runtime_error(
                        "Xios: Failed to set local y-size for domain '" + domainId + "'");
                }
                cxios_set_domain_jbegin(domain, jbegin);
                if (!cxios_is_defined_domain_jbegin(domain)) {
                    throw std::runtime_error(
                        "Xios: Failed to set local starting y-index for domain '" + domainId + "'");
                }
                std::vector<double> latvalue;
                for (int j = 0; j < nj; j++) {
                    latvalue.push_back(jbegin + j);
                }
                cxios_set_domain_latvalue_1d(domain, latvalue.data(), &nj);
                if (!cxios_is_defined_domain_latvalue_1d(domain)) {
                    throw std::runtime_error(
                        "Xios: Failed to set local y-indices for domain '" + domainId + "'");
                }
                cxios_set_domain_dim_j_name(domain, dimName.c_str(), dimName.length());
                if (!cxios_is_defined_domain_dim_j_name(domain)) {
                    throw std::runtime_error(
                        "Xios: Failed to set y-coordinate name for domain '" + domainId + "'");
                }
                cxios_set_domain_lat_name(domain, dimName.c_str(), dimName.length());
                if (!cxios_is_defined_domain_lat_name(domain)) {
                    throw std::runtime_error(
                        "Xios: Failed to set latitude name for domain '" + domainId + "'");
                }
            } else {
                throw std::runtime_error(
                    "Xios: More than 2 dimensions were associated with a domain.");
            }
            counter++;
        }
    }
}

/*!
 * @brief   Create XIOS grids for each ModelArray type
 *
 * @details This function sets up the XIOS grids for each field type based on the configuration
 *          in the gridIds, axisIds, and domainIds maps.
 */
void Xios::setupGrids()
{
    // Create XIOS grid associated with domain and possibly axis
    for (const auto& [type, gridId] : gridIds) {
        xios::CGrid* grid = getGrid(gridId);
        const std::string& domainId = domainIds[type];
        xios::CDomain* domain = getDomain(domainId);
        cxios_xml_tree_add_domaintogrid(grid, &domain, domainId.c_str(), domainId.length());
    }
}

/*!
 * Get the grid associated with a given ID
 *
 * @param the grid ID
 * @return a pointer to the XIOS CGrid object
 */
xios::CGrid* Xios::getGrid(const std::string& gridId)
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
xios::CField* Xios::getField(const std::string& fieldId)
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

// Split a string into a set by some delimiter.
std::set<std::string> str2set(const std::string& asStr, const char& delim = ',')
{
    std::set<std::string> asSet;
    if (asStr.length() > 0) {
        const char delim = ',';
        std::istringstream iss(asStr);
        std::string item;
        while (std::getline(iss, item, delim)) {
            asSet.insert(item);
        }
    }
    return asSet;
}

// Extract the field_names entry from the XiosInput section of the config.
std::set<std::string> Xios::configGetInputRestartFieldNames()
{
    return str2set(Configured::getConfiguration(keyMap.at(INPUT_FIELD_NAMES_KEY), std::string()));
}

/*!
 * Extract the field_names entry from the XiosForcing section of the config.
 */
std::set<std::string> Xios::configGetForcingFieldNames()
{
    return str2set(Configured::getConfiguration(keyMap.at(FORCING_FIELD_NAMES_KEY), std::string()));
}

// Extract the field_names entry from the XiosOutput section of the config.
std::set<std::string> Xios::configGetOutputRestartFieldNames()
{
    return str2set(Configured::getConfiguration(keyMap.at(OUTPUT_FIELD_NAMES_KEY), std::string()));
}

// Extract the field_names entry from the XiosDiagnostic section of the config.
std::set<std::string> Xios::configGetDiagnosticFieldNames()
{
    return str2set(
        Configured::getConfiguration(keyMap.at(DIAGNOSTIC_FIELD_NAMES_KEY), std::string()));
}

// Check whether a fieldId exists in a string of field names separated by commas, as determined by
// the map key
bool Xios::configCheckField(const std::string& fieldId, const bool& readAccess)
{
    std::set<std::string> fieldNames;
    if (readAccess) {
        fieldNames = configGetInputRestartFieldNames();
        std::set<std::string> forcingFieldNames = configGetForcingFieldNames();
        fieldNames.insert(forcingFieldNames.begin(), forcingFieldNames.end());
    } else {
        fieldNames = configGetOutputRestartFieldNames();
        std::set<std::string> diagnosticFieldNames = configGetDiagnosticFieldNames();
        fieldNames.insert(diagnosticFieldNames.begin(), diagnosticFieldNames.end());
    }
    return fieldNames.find(fieldId) != fieldNames.end();
}

/*!
 * Create a field with some ID
 *
 * @param the field ID
 */
void Xios::createField(const std::string& fieldId)
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

    // Set the operation type
    std::string operation;
    if (configGetInputRestartFieldNames().count(fieldId) > 0) {
        // Restarts are read "once"
        operation = "once";
    } else if (configGetDiagnosticFieldNames().count(fieldId) > 0) {
        // Diagonstics are averaged over the diagnostic output period
        operation = "average";
    } else {
        // Otherwise, read/write all timesteps without post-processing
        operation = "instant";
    }
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
void Xios::setFieldGridRef(const std::string& fieldId, const std::string& gridRef)
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
void Xios::setFieldReadAccess(const std::string& fieldId, const bool& readAccess)
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
void Xios::setFieldFreqOffset(const std::string& fieldId, const Duration& freqOffset)
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
 * Get the grid reference associated with a field with a given ID
 *
 * @param the field ID
 * @return grid reference used for the corresponding field
 */
std::string Xios::getFieldGridRef(const std::string& fieldId)
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
bool Xios::getFieldReadAccess(const std::string& fieldId)
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
Duration Xios::getFieldFreqOffset(const std::string& fieldId)
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
ModelArray::Type Xios::getFieldType(const std::string& fieldId) { return fieldTypes[fieldId]; }

/*!
 * Set the field type associated with a field with a given ID
 *
 * @param the field ID
 * @param ModelArray::Type used for the corresponding field
 */
void Xios::setFieldType(const std::string& fieldId, const ModelArray::Type& fieldType)
{
    fieldTypes[fieldId] = fieldType;
    setFieldGridRef(fieldId, gridIds[fieldType]);
}

/*!
 * @brief   Do an initial read of input files to deduce field dimensions.
 *
 * @details This function will read the dimension information from any NetCDF input files (restarts
 *          and/or forcings) and set dimensions appropriately. It will then set the field type of
 *          each input field.
 */
void Xios::setupFields()
{
    ModelMetadata& metadata = ModelMetadata::getInstance();

    for (const std::string& filename : { metadata.initialFileName, forcingFilename }) {
        if (filename.empty()) {
            break;
        }
        metadata.setDimensionsFromFile(filename);

        // Create map for field types
        const std::map<std::string, ModelArray::Type> dimensionKeys = {
            { "ydimxdim", ModelArray::Type::H },
            { "y_dimx_dim", ModelArray::Type::H },
            { "ydimxdimdg_comp", ModelArray::Type::DG },
            { "y_dimx_dimdg_comp", ModelArray::Type::DG },
            { "ydimxdimdgstress_comp", ModelArray::Type::DGSTRESS },
            { "y_dimx_dimdgstress_comp", ModelArray::Type::DGSTRESS },
            { "y_cgx_cg", ModelArray::Type::CG },
            { "yvertexxvertexncoords", ModelArray::Type::VERTEX },
            { "y_vertexx_vertexncoords", ModelArray::Type::VERTEX },
        };

        // Determine field types
        std::set<std::string> configFieldIds;
        if (filename == metadata.initialFileName) {
            configFieldIds = configGetInputRestartFieldNames();
        } else {
            configFieldIds = configGetForcingFieldNames();
        }
        try {
            auto& modelMPI = ModelMPI::getInstance();
            netCDF::NcFilePar ncFile(filename, netCDF::NcFile::read, modelMPI.getComm());

            for (auto& [fieldId, var] : ncFile.getVars()) {
                // Only consider fields that appear in the config
                if (configFieldIds.count(fieldId) == 0) {
                    continue;
                }
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
                setFieldType(fieldId, dimensionKeys.at(dimKey));
            }
            ncFile.close();
        } catch (const netCDF::exceptions::NcException& nce) {
            std::string ncWhat(nce.what());
            ncWhat += ": " + filename;
            throw std::runtime_error(ncWhat);
        }
    }
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
xios::CFile* Xios::getFile(const std::string& fileId)
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
void Xios::createFile(const std::string& fileId)
{
    if (!(fileId == outputFileId || fileId == inputFileId || fileId == diagnosticFileId
            || fileId == forcingFileId)) {
        throw std::runtime_error("Xios::createFile: Invalid fileId '" + fileId + "'");
    }

    // Deduce the field type
    int ioType = -1;
    for (const auto& [ioTypeOther, fileIdOther] : fileMap) {
        if (fileId == fileIdOther) {
            ioType = ioTypeOther;
            break;
        }
    }
    if (ioType == -1) {
        throw std::runtime_error(
            "Xios::createFile: Could not deduce file type for file '" + fileId + "'");
    }

    // Create the file
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

    // Set file name
    cxios_set_file_name(file, fileId.c_str(), fileId.length());
    if (!cxios_is_defined_file_name(file)) {
        throw std::runtime_error("Xios: Failed to set name for file '" + fileId + "'");
    }

    // Determine whether the file is configured for reading or writing
    bool readAccess = (ioType == INPUT_RESTART || ioType == FORCING);

    // Set the file mode
    std::string fileMode;
    if (readAccess) {
        fileMode = "read";
    } else {
        fileMode = "write";
    }
    if (cxios_is_defined_file_mode(file)) {
        Logged::warning("Xios: Overwriting mode for file '" + fileId + "'");
    }
    cxios_set_file_mode(file, fileMode.c_str(), fileMode.length());
    if (!cxios_is_defined_file_mode(file)) {
        throw std::runtime_error("Xios: Failed to set mode for file '" + fileId + "'");
    }

    // Set the file type to one_file
    const std::string fileType = "one_file";
    if (cxios_is_defined_file_type(file)) {
        Logged::warning("Xios: Overwriting type for file '" + fileId + "'");
    }
    cxios_set_file_type(file, fileType.c_str(), fileType.length());
    if (!cxios_is_defined_file_type(file)) {
        throw std::runtime_error("Xios: Failed to set type for file '" + fileId + "'");
    }

    // Set the file parallel access to collective
    const std::string parAccess = "collective";
    if (cxios_is_defined_file_par_access(file)) {
        Logged::warning("Xios: Overwriting parallel access for file '" + fileId + "'");
    }
    cxios_set_file_par_access(file, parAccess.c_str(), parAccess.length());
    if (!cxios_is_defined_file_par_access(file)) {
        throw std::runtime_error("Xios: Failed to set parallel access for file '" + fileId + "'");
    }

    // Get the fieldIds
    std::set<std::string> fieldIds;
    if (ioType == INPUT_RESTART) {
        fieldIds = configGetInputRestartFieldNames();
    } else if (ioType == OUTPUT_RESTART) {
        fieldIds = configGetOutputRestartFieldNames();
    } else if (ioType == FORCING) {
        fieldIds = configGetForcingFieldNames();
    } else if (ioType == DIAGNOSTIC) {
        fieldIds = configGetDiagnosticFieldNames();
    }

    // Determine the file output frequency
    cxios_duration outputFreq;
    ModelMetadata& metadata = ModelMetadata::getInstance();
    if (ioType == INPUT_RESTART || ioType == OUTPUT_RESTART) {
        outputFreq = convertDurationToXios(metadata.restartPeriod);
    } else {
        std::string periodStr;
        if (ioType == FORCING) {
            periodStr = Configured::getConfiguration(keyMap.at(FORCING_PERIOD_KEY), std::string());
        } else {
            periodStr
                = Configured::getConfiguration(keyMap.at(DIAGNOSTIC_PERIOD_KEY), std::string());
            cxios_set_file_split_freq(file, convertDurationToXios(Duration(periodStr)));
        }
        if (periodStr.empty() || periodStr == "0") {
            outputFreq = convertDurationToXios(metadata.runLength());
        } else {
            outputFreq = convertDurationToXios(Duration(periodStr));
        }
    }

    // Set the file output frequency
    if (cxios_is_defined_file_output_freq(file)) {
        Logged::warning("Xios: Overwriting output frequency for file '" + fileId + "'");
    }
    cxios_set_file_output_freq(file, outputFreq);
    if (!cxios_is_defined_file_output_freq(file)) {
        throw std::runtime_error("Xios: Failed to set output frequency for file '" + fileId + "'");
    }

    // Set the file split frequency to coincide with the output frequency for output files
    if (!readAccess) {
        if (cxios_is_defined_file_split_freq(file)) {
            Logged::warning("Xios: Split frequency already set for file '" + fileId + "'");
        }
        cxios_set_file_split_freq(file, outputFreq);
        if (!cxios_is_defined_file_split_freq(file)) {
            throw std::runtime_error(
                "Xios: Failed to set split frequency for file '" + fileId + "'");
        }

        // Set format string for file splitting, converting characters as expected by XIOS
        std::string split_freq_format;
        if (ioType == OUTPUT_RESTART) {
            split_freq_format = outputFormatStr;
        } else {
            split_freq_format = diagnosticFormatStr;
        }
        for (const auto& [from, to] : formatStrMap) {
            auto pos = split_freq_format.find(from);
            if (pos != std::string::npos) {
                split_freq_format.replace(pos, from.length(), to);
            }
        }
        cxios_set_file_split_freq_format(
            file, split_freq_format.c_str(), split_freq_format.length());
        if (!cxios_is_defined_file_split_freq_format(file)) {
            throw std::runtime_error(
                "Xios: Failed to set split frequency format for file '" + fileId + "'");
        }
    }

    // Loop over field_names entries in the config
    for (const std::string& fieldId : fieldIds) {
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
 * Get all field IDs associated with a given file
 *
 * @param the file ID
 * @return all field IDs associated with the file
 */
std::vector<std::string> Xios::fileGetFieldIds(const std::string& fileId)
{
    std::vector<xios::CField*> fields = getFile(fileId)->getAllFields();
    std::vector<std::string> fieldIds(fields.size());
    for (size_t i = 0; i < fields.size(); i++) {
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
void Xios::fileAddField(const std::string& fileId, const std::string& fieldId)
{
    xios::CField* field = getField(fieldId);
    cxios_xml_tree_add_fieldtofile(getFile(fileId), &field, fieldId.c_str(), fieldId.length());
}

/*!
 * Set up files based on the configuration.
 */
void Xios::setupFiles()
{
    auto& metadata = ModelMetadata::getInstance();

    // Get restart file IDs from the configuration
    inputFileId = ((std::filesystem::path)metadata.initialFileName).filename().replace_extension();
    outputFileId = ((std::filesystem::path)metadata.finalFileName).filename().replace_extension();
    if (outputFileId.find("%") != std::string::npos) {
        outputFormatStr = outputFileId.substr(outputFileId.find("%"), outputFileId.find(".nc"));
        outputFileId.erase(outputFileId.find("%"), outputFileId.length());
    }
    if (!inputFileId.empty() && inputFileId == outputFileId) {
        throw std::runtime_error("Xios::setupFiles: Input and restart file names must differ.");
    }

    // Get forcing file name and ID from the configuration
    forcingFilename = Configured::getConfiguration(keyMap.at(FORCING_FILE_KEY), std::string());
    forcingFileId = ((std::filesystem::path)forcingFilename).filename().replace_extension();
    if (!forcingFileId.empty()) {
        if (inputFileId == forcingFileId) {
            throw std::runtime_error("Xios::setupFiles: Input and forcing file names must differ.");
        }
        if (outputFileId == forcingFileId) {
            throw std::runtime_error(
                "Xios::setupFiles: Restart and forcing file names must differ.");
        }
    }

    // Get diagnostic file name and ID from the configuration
    diagnosticFilename
        = Configured::getConfiguration(keyMap.at(DIAGNOSTIC_FILE_KEY), std::string());
    diagnosticFileId = ((std::filesystem::path)diagnosticFilename).filename().replace_extension();
    if (diagnosticFileId.find("%") != std::string::npos) {
        diagnosticFormatStr
            = diagnosticFileId.substr(diagnosticFileId.find("%"), diagnosticFileId.find(".nc"));
        diagnosticFileId.erase(diagnosticFileId.find("%"), diagnosticFileId.length());
    }
    if (!diagnosticFileId.empty()) {
        if (inputFileId == diagnosticFileId) {
            throw std::runtime_error(
                "Xios::setupFiles: Input and diagnostic file names must differ.");
        }
        if (outputFileId == diagnosticFileId) {
            throw std::runtime_error(
                "Xios::setupFiles: Restart and diagnostic file names must differ.");
        }
        if (forcingFileId == diagnosticFileId) {
            throw std::runtime_error(
                "Xios::setupFiles: Forcing and diagnostic file names must differ.");
        }
    }

    // Create files for any non-empty file IDs
    for (const std::string& fileId :
        { outputFileId, inputFileId, diagnosticFileId, forcingFileId }) {
        if (!fileId.empty()) {
            createFile(fileId);
        }
    }
}

/*!
 * @brief   Postprocess output files after the simulation has completed.
 *
 * @details If only a single domain was written to file, rename the x and y dimensions and variables
 *          to x_dim and y_dim, respectively, for compatibility with other model components.
 */
void Xios::postprocessOutputFiles()
{
    // Count how many domains were written out
    int sum = 0;
    for (const auto& [domainId, written] : domainWritten) {
        if (written) {
            sum++;
        }
    }

    // If a single domain was written then modify x and y dimensions and variables in the output
    // files
    if (sum == 1 && mpi_rank == 0) {

        // Consider both restart files and diagnostic files
        for (std::string fileId : { outputFileId, diagnosticFileId }) {
            bool exists;
            cxios_file_valid_id(&exists, fileId.c_str(), fileId.length());
            if (!exists) {
                continue;
            }

            // Get file output frequency
            xios::CFile* file = getFile(fileId);
            if (!cxios_is_defined_file_output_freq(file)) {
                throw std::runtime_error(
                    "Xios: Undefined output frequency for file '" + fileId + "'");
            }
            cxios_duration duration;
            cxios_get_file_output_freq(file, &duration);
            Duration step = convertDurationFromXios(duration);

            // Deduce the format string
            std::string formatStr;
            if (fileId == outputFileId) {
                formatStr = outputFormatStr;
            } else {
                formatStr = diagnosticFormatStr;
            }

            // Loop over the output window splits
            ModelMetadata& metadata = ModelMetadata::getInstance();
            TimePoint time = metadata.startTime();
            TimePoint endTime = metadata.stopTime();
            while (time < endTime) {

                // Compute the end time of the window, subtracting 1 second to avoid overlap
                TimePoint nextTime = time + step - Duration(1);

                // Generate the filename used by XIOS
                std::string filename = fileId + "_" + time.format(formatStr) + "-"
                    + nextTime.format(formatStr) + ".nc";
                std::cout << "DEBUG: filename=" << filename << std::endl;

                // Increment the time then check if the file exists
                time += step;
                if (!std::filesystem::exists(filename)) {
                    continue;
                }
                std::cout << "DEBUG: file exists" << std::endl;

                try {
                    // Open the netCDF file for both reading and writing
                    netCDF::NcFile ncFile(filename, netCDF::NcFile::write);

                    // Rename the x and y dimensions with x_dim and y_dim, respectively
                    ncFile.getDim("x").rename("x_dim");
                    ncFile.getDim("y").rename("y_dim");

                    // Rename the x and y variables with x_dim and y_dim, respectively
                    ncFile.getVar("x").rename("x_dim");
                    ncFile.getVar("y").rename("y_dim");

                    // Ensure changes are flushed to disk before closing
                    ncFile.sync();
                } catch (const netCDF::exceptions::NcException& e) {
                    std::cerr << "Error processing NetCDF file: " << e.what() << std::endl;
                }
            }
        }
    }
}

/*!
 * Send a field to the XIOS server to be written to file
 *
 * @param field name
 * @param reference to the ModelArray containing the data to be written
 */
void Xios::write(const std::string& fieldId, const ModelArray& modelarray)
{
    if (getFieldReadAccess(fieldId)) {
        throw std::runtime_error("Xios::write: field " + fieldId
            + " is not configured for writing, but is being written to file.");
    };
    if (modelarray.nDimensions() != 2) {
        throw std::invalid_argument("Only ModelArrays of dimension 2 are supported");
    }
    auto& dims = modelarray.dimensions();
    const ModelArray::Type& type = modelarray.getType();
    domainWritten[domainIds[type]] = true;
    if ((type == ModelArray::Type::H) || (type == ModelArray::Type::CG)) {
        cxios_write_data_k82(
            fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0], dims[1], -1);
    } else if (type == ModelArray::Type::VERTEX) {
        cxios_write_data_k83(fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0],
            dims[1], ModelArray::size(ModelArray::Dimension::NCOORDS), -1);
    } else if (type == ModelArray::Type::DG) {
        cxios_write_data_k83(fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0],
            dims[1], ModelArray::size(ModelArray::Dimension::DG), -1);
    } else if (type == ModelArray::Type::DGSTRESS) {
        cxios_write_data_k83(fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0],
            dims[1], ModelArray::size(ModelArray::Dimension::DGSTRESS), -1);
    } else {
        throw std::invalid_argument(
            "Only HFields, VertexFields, DGFields, DGSFields, and CGFields are supported");
    }
}

/*!
 * Receive field from XIOS server that has been read from file.
 *
 * @param field name
 * @param reference to the ModelArray containing the data to be written
 */
void Xios::read(const std::string& fieldId, ModelArray& modelarray)
{
    if (!getFieldReadAccess(fieldId)) {
        throw std::runtime_error("Xios::read: field " + fieldId
            + " is not configured for reading, but is being read from file.");
    };
    if (modelarray.nDimensions() != 2) {
        throw std::invalid_argument("Only ModelArrays of dimension 2 are supported");
    }
    auto& dims = modelarray.dimensions();
    const ModelArray::Type& type = modelarray.getType();
    if ((type == ModelArray::Type::H) || (type == ModelArray::Type::CG)) {
        cxios_read_data_k82(
            fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0], dims[1]);
    } else if (type == ModelArray::Type::VERTEX) {
        cxios_read_data_k83(fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0],
            dims[1], ModelArray::size(ModelArray::Dimension::NCOORDS));
    } else if (type == ModelArray::Type::DG) {
        cxios_read_data_k83(fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0],
            dims[1], ModelArray::size(ModelArray::Dimension::DG));
    } else if (type == ModelArray::Type::DGSTRESS) {
        cxios_read_data_k83(fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0],
            dims[1], ModelArray::size(ModelArray::Dimension::DGSTRESS));
    } else {
        throw std::invalid_argument(
            "Only HFields, VertexFields, DGFields, DGSFields, and CGFields are supported");
    }
}
}

#endif
