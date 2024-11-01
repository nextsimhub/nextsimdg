/*!
 * @file    ParaGridIO_Xios.cpp
 *
 * @date    01 Nov 2024
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 */
#ifdef USE_XIOS
#include "include/ParaGridIO_Xios.hpp"
#include <include/xios_c_interface.hpp>

#include <filesystem>

namespace Nextsim {

ParaGridIO::ParaGridIO(ParametricGrid& grid)
    : IParaGridIO(grid)
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
#ifndef USE_MPI
ModelState ParaGridIO::getModelState(const std::string& filePath)
{
    ModelState state;
    throw std::invalid_argument("XIOS cannot be used without MPI");
    return state;
#else
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
    throw std::runtime_error("XIOS implementation of getModelState incomplete"); // TODO-JGW
    return state;
#endif
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
    std::string fileId = ((std::filesystem::path)filePath).replace_extension();
    xiosHandler.createFile(fileId);
    xiosHandler.setFileType(fileId, "one_file");
    xiosHandler.setFileMode(fileId, "write");
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
