/*!
 * @file   ParaGridIO_Xios.cpp
 *
 * @date   27 Aug, 2024
 * @author Joe Wallwork <jw2423@cam.ac.uk>
 */
#ifdef USE_XIOS
#include "include/ParaGridIO_Xios.hpp"
#include <include/xios_c_interface.hpp>

namespace Nextsim {
ParaGridIO::~ParaGridIO() = default;

// TODO: Move to ParaGridIO API below
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

// TODO: Move to ParaGridIO API below
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
#ifdef USE_MPI
ModelState ParaGridIO::getModelState(const std::string& filePath, ModelMetadata& metadata)
{
    // TODO: Implement
    ModelState ms;
    return ms;
}
#else
ModelState ParaGridIO::getModelState(const std::string& filePath)
{
    // TODO: Implement
    ModelState ms;
    return ms;
}
#endif

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
    // TODO: Implement
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
    // TODO: Implement
    ModelState ms;
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
    // TODO: Implement
}
};
#endif /* USE_XIOS */
