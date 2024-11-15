/*!
 * @file    ParaGridIO_Xios.cpp
 *
 * @date    15 Nov 2024
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 */
#ifdef USE_XIOS
#include "include/ParaGridIO_Xios.hpp"
#include "include/NZLevels.hpp"
#include "include/gridNames.hpp"

#include <include/xios_c_interface.hpp>

#include <filesystem>

namespace Nextsim {

ParaGridIO::ParaGridIO(ParametricGrid& grid, const std::string contextId)
    : IParaGridIO(grid)
    , xiosHandler(contextId)
{
    // TODO-JGW: Implement this method
}

ParaGridIO::~ParaGridIO() = default;

// TODO-JGW: Move to ParaGridIO API below
/*!
 * Send a field to the XIOS server to be written to file
 *
 * @param field name
 * @param reference to the ModelArray containing the data to be written
 */
void ParaGridIO::write(const std::string fieldId, ModelArray& modelarray)
{
    auto ndim = modelarray.nDimensions();
    auto dims = modelarray.dimensions();
    if (ndim == 2) {
        cxios_write_data_k82(
            fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0], dims[1], -1);
    } else if (ndim == 3) {
        cxios_write_data_k83(
            fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0], dims[1], dims[2], -1);
    } else if (ndim == 4) {
        cxios_write_data_k84(fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0],
            dims[1], dims[2], dims[3], -1);
    } else {
        throw std::invalid_argument("Only ModelArrays of dimension 2, 3, or 4 are supported");
    }
}

// TODO-JGW: Move to ParaGridIO API below
/*!
 * Receive field from XIOS server that has been read from file.
 *
 * @param field name
 * @param reference to the ModelArray containing the data to be written
 */
void ParaGridIO::read(const std::string fieldId, ModelArray& modelarray)
{
    auto ndim = modelarray.nDimensions();
    auto dims = modelarray.dimensions();
    if (ndim == 2) {
        cxios_read_data_k82(
            fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0], dims[1]);
    } else if (ndim == 3) {
        cxios_read_data_k83(
            fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0], dims[1], dims[2]);
    } else if (ndim == 4) {
        cxios_read_data_k84(fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0],
            dims[1], dims[2], dims[3]);
    } else {
        throw std::invalid_argument("Only ModelArrays of dimension 2, 3, or 4 are supported");
    }
}

/*!
 * Setup the Axis, Domain, Grid, Field, and File attributes required by XIOS based on the
 * ModelState, ModelMetadata, and FilePath provided.
 *
 * @params state The model state and configuration object.
 * @params metadata The model metadata (principally the initial file creation model time).
 * @params filePath The path for the restart file.
 */
void ParaGridIO::setupXios(
    const ModelState& state, const ModelMetadata& meta, const std::string& filePath)
{
    if (_xiosSetup) {
        return;
    }

    // Setup XIOS Domain attribute with special domainId 'xy_domain' that creates grid_2D along
    // with it
    xiosHandler.createDomain("xy_domain");
    xiosHandler.setDomainType("xy_domain", "rectilinear");
    for (auto entry : ModelArray::definedDimensions) {
        auto dimType = entry.first;
        ModelArray::DimensionSpec& dimensionSpec = entry.second;
        const std::string name = dimensionSpec.name;
        std::cout << "DEBUG setupXios: found dimension with name: " << name << std::endl;
        if (name == "xdim") {
            std::cout << "DEBUG setupXios: setting x-attrs of xy_domain: " << name << std::endl;
            xiosHandler.setDomainGlobalXSize("xy_domain", dimensionSpec.globalLength);
            xiosHandler.setDomainLocalXSize("xy_domain", dimensionSpec.localLength);
            xiosHandler.setDomainLocalXStart("xy_domain", dimensionSpec.start);
        } else if (name == "ydim") {
            std::cout << "DEBUG setupXios: setting y-attrs of xy_domain: " << name << std::endl;
            xiosHandler.setDomainGlobalYSize("xy_domain", dimensionSpec.globalLength);
            xiosHandler.setDomainLocalYSize("xy_domain", dimensionSpec.localLength);
            xiosHandler.setDomainLocalYStart("xy_domain", dimensionSpec.start);
        } else if (name == "zdim") {
            // Setup XIOS Axis attribute with special axisId 'z_axis' that creates grid_3D along
            // with it
            std::cout << "DEBUG setupXios: setting attrs of z_axis: " << name << std::endl;
            xiosHandler.createAxis("z_axis");
            xiosHandler.setAxisSize("z_axis", dimensionSpec.globalLength);
            if (dimensionSpec.globalLength != dimensionSpec.localLength) {
                throw std::runtime_error("ParaGridIO_Xios: Inconsistent dimensionSpec for z-axis");
            }
        }
        // TODO-JGW: What about *vertex, *_cg, dg_comp, dgstress_comp, ncoords?
    }

    // Setup XIOS File attribute
    std::string fileId = ((std::filesystem::path)filePath).replace_extension();
    xiosHandler.createFile(fileId);
    xiosHandler.setFileType(fileId, "one_file");
    xiosHandler.setFileMode(fileId, "write");
    Duration timestep = xiosHandler.getCalendarTimestep();
    xiosHandler.setFileOutputFreq(fileId, timestep); // TODO-JGW: Set actual output freq
    // xiosHandler.setFileSplitFreq(fileId, Duration(...));

    // Loop over fields in the ModelState and create XIOS Field attributes for each
    for (auto entry : state.data) {
        std::string fieldId = entry.first;
        ModelArray field = entry.second;

        // Set local x- and y-values of the XIOS Domain
        if (fieldId == "x") {
            std::vector<double> xValues {};
            std::cout << "DEBUG setupXios: x field data:" << std::endl;
            // TODO-JGW: Look up how ordering works
            for (int i = meta.localCornerX; i < meta.localCornerX + meta.localExtentX; i += 2) {
                std::cout << field.data()(i, 0) << ", " << std::endl;
                xValues.push_back(field.data()(i, 0));
            }
            xiosHandler.setDomainLocalXValues("xy_domain", xValues);
            continue;
        } else if (fieldId == "y") {
            std::vector<double> yValues {};
            std::cout << "DEBUG setupXios: y field data:" << std::endl;
            // TODO-JGW: Look up how ordering works
            for (int i = meta.localCornerY; i < meta.localCornerY + meta.localExtentY; i += 2) {
                std::cout << field.data()(i, 0) << ", " << std::endl;
                yValues.push_back(field.data()(i, 0));
            }
            xiosHandler.setDomainLocalYValues("xy_domain", yValues);
            continue;
        } else if (fieldId == "z") {
            // xiosHandler.setAxisValues("z_axis", ...); // TODO-JGW
            continue;
        }
        std::cout << "DEBUG setupXios: Creating field with fieldId=" << fieldId << std::endl;

        // Setup XIOS Field attribute and associate it with the File
        xiosHandler.createField(fieldId);
        xiosHandler.setFieldOperation(fieldId, "instant");
        xiosHandler.setFieldGridRef(fieldId, "grid_2D");
        // xiosHandler.setFieldGridRef(fieldId, "grid_3D"); // TODO-JGW: Account for 3D fields
        xiosHandler.setFieldReadAccess(fieldId, false);
        xiosHandler.fileAddField(fileId, fieldId);
    }

    // Mark XIOS setup complete
    xiosHandler.close_context_definition();
    _xiosSetup = true;
}

/*!
 * Retrieves the ModelState from a restart file of the parametric_grid type.
 * @param filePath The file path containing the file to be read.
 */
ModelState ParaGridIO::getModelState(const std::string& filePath, ModelMetadata& metadata)
{
    // TODO-JGW: Use setupXios
    ModelState state;

    if (!std::filesystem::exists(filePath)) {
        throw std::invalid_argument("ParaGridIO_Xios: File " + filePath + " does not exist");
    }
    std::string fileId = ((std::filesystem::path)filePath).replace_extension();
    xiosHandler.createFile(fileId);
    xiosHandler.setFileType(fileId, "one_file");
    xiosHandler.setFileMode(fileId, "read");

    // Determine Axes and Domains
    std::cout << "DEBUG getModelState: entering for loop" << std::endl;
    for (auto entry : ModelArray::definedDimensions) {
        auto dimType = entry.first;
        // if (dimCompMap.count(dimType) > 0)
        //     // TODO Assertions that DG in the file equals the compile time DG in the model. See
        //     // #205
        //     continue;

        // TODO: WIP
        ModelArray::DimensionSpec& dimensionSpec = entry.second;
        const std::string name = dimensionSpec.name;
        const std::string altName = dimensionSpec.altName;
        const size_t globalLength = dimensionSpec.globalLength;
        const size_t localLength = dimensionSpec.localLength;
        const size_t start = dimensionSpec.start;
        std::cout << "DEBUG getModelState: name: " << name << std::endl;
        std::cout << "DEBUG getModelState: altName: " << altName << std::endl;

        // // Find dimensions in the netCDF file by their name in the ModelArray details
        // netCDF::NcDim dim = dataGroup.getDim(dimensionSpec.name);
        // // Also check the old name
        // if (dim.isNull()) {
        //     dim = dataGroup.getDim(dimensionSpec.altName);
        // }
        // // If we didn't find a dimension with the dimensions name or altName, throw.
        // if (dim.isNull()) {
        //     throw std::out_of_range(
        //         std::string("No netCDF dimension found corresponding to the dimension named ")
        //         + dimensionSpec.name + std::string(" or ") + dimensionSpec.altName);
        // }
        if (dimType == ModelArray::Dimension::Z) {
            // A special case, as the number of levels in the file might not be
            // the number that the selected ice thermodynamics requires.
            const size_t nz = NZLevels::get();
            ModelArray::setDimension(dimType, nz, nz, 0);
            xiosHandler.createAxis("z_axis");
            xiosHandler.setAxisSize("z_axis", nz);
            // TODO: xiosHandler.setAxisValues("z_axis", ...);
        } else {
            // auto dimName = dim.getName();
            size_t globalLength = 0;
            size_t localLength = 0;
            size_t start = 0;
            if (dimType == ModelArray::Dimension::X) {
                globalLength = metadata.globalExtentX;
                localLength = metadata.localExtentX;
                start = metadata.localCornerX;
                xiosHandler.createDomain("xy_domain"); // TODO-JGW: Check for existence
                xiosHandler.setDomainGlobalXSize("xy_domain", globalLength);
                xiosHandler.setDomainLocalXSize("xy_domain", localLength);
                xiosHandler.setDomainLocalXStart("xy_domain", start);
            } else if (dimType == ModelArray::Dimension::Y) {
                globalLength = metadata.globalExtentY;
                localLength = metadata.localExtentY;
                start = metadata.localCornerY;
                // xiosHandler.createDomain("xy_domain"); // TODO-JGW: Check for existence
                xiosHandler.setDomainGlobalYSize("xy_domain", globalLength);
                xiosHandler.setDomainLocalYSize("xy_domain", localLength);
                xiosHandler.setDomainLocalYStart("xy_domain", start);
            } else if (dimType == ModelArray::Dimension::XVERTEX) {
                globalLength = metadata.globalExtentX + 1;
                localLength = metadata.localExtentX + 1;
                start = metadata.localCornerX;
                // TODO-JGW: set up XIOS attributes
            } else if (dimType == ModelArray::Dimension::YVERTEX) {
                globalLength = metadata.globalExtentY + 1;
                localLength = metadata.localExtentY + 1;
                start = metadata.localCornerY;
                // TODO-JGW: set up XIOS attributes
            } else {
                // localLength = dim.getSize();
                start = 0;
                throw std::runtime_error("dimType not yet accounted for"); // TODO-JGW
            }
            ModelArray::setDimension(dimType, globalLength, localLength, start);
        }
    }
    std::cout << "DEBUG getModelState: exiting for loop" << std::endl << std::endl;
    // TODO: Set fields?
    // TODO: Call read method
    throw std::runtime_error("XIOS implementation of getModelState incomplete"); // TODO-JGW
    return state;
}

ModelState ParaGridIO::readForcingTimeStatic(
    const std::set<std::string>& forcings, const TimePoint& time, const std::string& filePath)
{
    ModelState state;
    throw std::runtime_error("XIOS implementation of readForcingTimeStatic incomplete"); // TODO-JGW
    return state;
}

/*!
 * @brief Writes the ModelState to a given file location from the provided
 * model data and metadata.
 *
 * @params state The model state and configuration object.
 * @params metadata The model metadata (principally the initial file
 * creation model time).
 * @params filePath The path for the restart file.
 */
void ParaGridIO::dumpModelState(
    const ModelState& state, const ModelMetadata& meta, const std::string& filePath)
{
    // Setup the XIOS context if it hasn't been already
    setupXios(state, meta, filePath);

    std::set<std::string> restartFields = { hiceName, ciceName, hsnowName, ticeName, sstName,
        sssName, maskName, coordsName, xName, yName, longitudeName, latitudeName, gridAzimuthName,
        uName, vName, damageName }; // TODO and others
    // If the above fields are found in the supplied ModelState, output them
    std::cout << "DEBUG dumpModelState: time=" << meta.time() << std::endl;
    std::cout << "DEBUG dumpModelState: entering for loop" << std::endl;
    for (auto entry : state.data) {
        if (restartFields.count(entry.first)) {
            std::string fieldId = entry.first;
            ModelArray field = entry.second;
            std::cout << "DEBUG dumpModelState: fieldId=" << fieldId << " has name "
                      << xiosHandler.getFieldName(fieldId) << std::endl;

            // Send data to XIOS to be written to disk
            // FIXME: Doesn't throw an error, but doesn't write a file either. May be because Axes,
            // Domains and Grids weren't created but may also be because contexts haven't been
            // handled correctly?
            std::cout << "DEBUG dumpModelState: writing fieldId=" << fieldId << std::endl;
            write(fieldId, field);
            std::cout << "DEBUG dumpModelState: wrote fieldId=" << fieldId << std::endl;
        }
    }
    std::cout << "DEBUG dumpModelState: exiting for loop" << std::endl;
    throw std::runtime_error("XIOS implementation of dumpModelState incomplete"); // TODO-JGW
}

/*!
 * @brief Reads forcings from a ParameticGrid flavoured file.
 *
 * @param forcings The names of the forcings required.
 * @param time The time for which to get the forcings.
 * @param filePath Path to the file to read.
 */
ModelState ParaGridIO::readForcingTime(
    const std::set<std::string>& forcings, const TimePoint& time, const std::string& filePath)
{
    ModelState ms;
    throw std::runtime_error("XIOS implementation of readForcingTime incomplete"); // TODO-JGW
    return ms;
}

/*!
 * @brief Writes diagnostic data to a file.
 *
 * @param state The state to write to the file.
 * @param time The time of the passed data.
 * @param filePath Path of the file to write to.
 */
void ParaGridIO::writeDiagnosticTime(
    const ModelState& state, const ModelMetadata& meta, const std::string& filePath)
{
    throw std::runtime_error("XIOS implementation of writeDiagnosticTime incomplete"); // TODO-JGW
}
};
#endif /* USE_XIOS */
