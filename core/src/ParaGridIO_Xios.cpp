/*!
 * @file    ParaGridIO_Xios.cpp
 *
 * @date    12 Nov 2024
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
 * Retrieves the ModelState from a restart file of the parametric_grid type.
 * @param filePath The file path containing the file to be read.
 */
ModelState ParaGridIO::getModelState(const std::string& filePath, ModelMetadata& metadata)
{
    ModelState state;

    if (!std::filesystem::exists(filePath)) {
        throw std::invalid_argument("ParaGridIO_Xios: File " + filePath + " does not exist");
    }
    std::string fileId = ((std::filesystem::path)filePath).replace_extension();
    xiosHandler.createFile(fileId);
    xiosHandler.setFileType(fileId, "one_file");
    xiosHandler.setFileMode(fileId, "read");

    // Determine Axes and Domains
    std::cout << "enter for loop" << std::endl;
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
        std::cout << "name: " << name << std::endl;
        std::cout << "altName: " << altName << std::endl;

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
    std::cout << "exit for loop" << std::endl << std::endl;
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
    // Setup XIOS File attribute
    std::string fileId = ((std::filesystem::path)filePath).replace_extension();
    xiosHandler.createFile(fileId);
    xiosHandler.setFileType(fileId, "one_file");
    xiosHandler.setFileMode(fileId, "write");
    // xiosHandler.setFileOutputFreq(fieldId, timestep)
    // xiosHandler.setFileSplitFreq(fieldId, Duration(...))

    // NOTE: ModelState.data provides a map from strings to ModelArrays

    // TODO: Deduce Axes, Domains and Grids

    std::set<std::string> restartFields = { hiceName, ciceName, hsnowName, ticeName, sstName,
        sssName, maskName, coordsName, xName, yName, longitudeName, latitudeName, gridAzimuthName,
        uName, vName, damageName }; // TODO and others
    // If the above fields are found in the supplied ModelState, output them
    for (auto entry : state.data) {
        if (restartFields.count(entry.first)) {
            std::string fieldId = entry.first;
            std::cout << "fieldId: " << fieldId << std::endl;
            ModelArray field = entry.second;

            // Setup XIOS Field attribute and associate it with the File
            xiosHandler.createField(fieldId);
            xiosHandler.setFieldOperation(fieldId, "instant");
            // xiosHandler.setFieldGridRef(fieldId, ...);
            xiosHandler.setFieldReadAccess(fieldId, false);
            xiosHandler.fileAddField(fileId, fieldId);

            // Send data to XIOS to be written to disk
            // FIXME: Doesn't throw an error, but doesn't write a file either. May be because Axes,
            // Domains and Grids weren't created but may also be because contexts haven't been
            // handled correctly?
            write(fieldId, field);
        }
    }
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
