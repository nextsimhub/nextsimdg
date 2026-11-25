/*!
 * @author  Einar Olason <einar.olason@nersc.no>
 */

#include "include/PlainInput.hpp"

#include <boost/log/sinks/text_file_backend.hpp>
#include <ncDim.h>
#include <ncVar.h>

namespace Nextsim {

void PlainInput::initialise(const TimePoint& time, const std::string& pathSpecIn,
    const std::string& ncLonNameIn, const std::string& ncLatNameIn, const std::string& ncTimeNameIn,
    const std::set<std::string>& forcingsIn,
    const std::set<std::pair<std::string, std::string>>& vectorsIn, const ModelArray& modelLonsIn,
    const ModelArray& modelLatsIn)
{
    forcings = forcingsIn;
    filePath = pathSpecIn;
    initTime = time;

    // Deduce the frequency from the first two elements of the time vector
    try {
        // Open the file
        const netCDF::NcFile ncFile(filePath, netCDF::NcFile::read);
        // Read the time variable
        const netCDF::NcVar timeVar = ncFile.getVar(ncTimeNameIn);
        // Get the time axis as a vector. We use double here, because Duration expects double.
        std::vector<double> timeVec(ncFile.getDim(ncTimeNameIn).getSize());
        timeVar.getVar(timeVec.data());

        if (timeVec.size() <= 1)
            throw std::runtime_error("PlainInput::initialise: Time axis in " + filePath
                + " too short (" + std::to_string(timeVec.size()) + ").\n");

        deltaTime = timeVec[1] - timeVec[0];

    } catch (const netCDF::exceptions::NcException& nce) {
        std::string ncWhat(nce.what());
        ncWhat += ": " + filePath;
        throw std::runtime_error("PlainInput::initialise: " + ncWhat + "\n");
    }
}

void PlainInput::update(const TimePoint& time)
{
    if (std::fmod((time - initTime).seconds(), deltaTime) == 0.0_ft) {
        ModelState state = ParaGridIO::readForcingTimeStatic(forcings, time, filePath);
#ifdef USE_KOKKOS
        std::map<std::string, DeviceViewMA> deviceState;
        for (auto& [name, field] : state.data) {
            deviceState.emplace(name,
                makeKokkosDeviceView<MakeViewOptions::ALWAYS_COPY>(
                    "PlainInput " + name, field.getDataRef()));
        }
        forcingState = std::move(deviceState);
#else
        forcingState = std::move(state.data);
#endif
    }
}

void PlainInput::getField(const std::string& fieldName, ModelArrayAuto& destination) const
{
    try {
#ifdef USE_KOKKOS
        Kokkos::deep_copy(destination.deviceView(), forcingState.at(fieldName));
#else
        destination = forcingState.at(fieldName);
#endif
    } catch (const std::out_of_range& e) {
        throw std::runtime_error("PlainInput::getField: Field '" + fieldName
            + "' not found in file '" + filePath + "'.\n");
    }
}

[[nodiscard]] bool PlainInput::isValid(const std::string& fieldName) const
{
    return forcingState.count(fieldName);
}

}
