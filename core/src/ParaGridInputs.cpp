/*!
 * @author  Einar Olason on 21/07/2026.
 */

#include <ncDim.h>
#include <ncFile.h>
#include <ncVar.h>

#include "include/NetCDFUtils.hpp"
#include "include/ParaGridInputs.hpp"

namespace Nextsim {

void ParaGridInputs::configure(const std::string& filePathIn, const std::string& ncLonNameIn,
    const std::string& ncLatNameIn, const std::string& ncTimeNameIn,
    const std::set<std::string>& forcingsIn,
    const std::set<std::pair<std::string, std::string>>& vectorsIn, const ModelArray& modelLonsIn,
    const ModelArray& modelLatsIn)
{
    filePath = filePathIn;
    forcings = forcingsIn;
    vectors = vectorsIn;

    ncTimeName = ncTimeNameIn;
    ncLatName = ncLatNameIn;
    ncLonName = ncLonNameIn;

    modelLons = modelLonsIn;
    modelLats = modelLatsIn;
}

void ParaGridInputs::update(const TimePoint& time)
{
    currentTime = time;

    if (time < timeRange.before || time > timeRange.after) {
        std::map<std::string, std::vector<FloatType>> rawDataBefore, rawDataAfter;
        readRawForcing(rawDataBefore, rawDataAfter);

        rotateInputVectors(rawDataBefore);
        rotateInputVectors(rawDataAfter);

        forcingStateBefore = interpolateSpatially(rawDataBefore);
        forcingStateAfter = interpolateSpatially(rawDataAfter);
    }
}

ModelArray ParaGridInputs::getField(const std::string& fieldName)
{
    ModelArray ma;
    ma.reinitialize();

    const double frac = (currentTime - timeRange.before).seconds()
        / (timeRange.after - timeRange.before).seconds();

    for (size_t i = 0; i < ma.size(); ++i) {
        ma[i] = frac * forcingStateAfter.data[fieldName][i]
            + (1. - frac) * forcingStateBefore.data[fieldName][i];
    }

    return ma;
}

ModelState ParaGridInputs::interpolateSpatially(
    const std::map<std::string, std::vector<FloatType>>& rawData)
{
    ModelState state;
    for (const auto& [name, data] : rawData) {
        state.data[name].reinitialize();
        for (size_t i = 0; i < data.size(); ++i) {
            state.data[name][i] = data[i]; // TODO: Actually interpolate
        }
    }

    return state;
}

void ParaGridInputs::rotateInputVectors(std::map<std::string, std::vector<FloatType>>& rawData)
{
    auto rawDataCopy = rawData;
    for (auto& pair : vectors) {
        vectorRotationLogic(rawDataCopy[pair.first], rawDataCopy[pair.second], rawData[pair.first],
            rawDataCopy[pair.second]);
    }
}

void ParaGridInputs::vectorRotationLogic(const std::vector<FloatType>& vectorIn1st,
    const std::vector<FloatType>& vectorIn2nd, std::vector<FloatType>& vectorOut1st,
    std::vector<FloatType>& vectorOut2nd)
{
    // Nothing done yet
    vectorOut1st = vectorIn1st;
    vectorOut2nd = vectorIn2nd;
}

void ParaGridInputs::readRawForcing(std::map<std::string, std::vector<FloatType>>& rawDataBefore,
    std::map<std::string, std::vector<FloatType>>& rawDataAfter)
{
    const std::string fileNameBefore = currentTime.format(filePath);
    std::string fileNameAfter = fileNameBefore;
    size_t targetTIndexAfter, targetTIndexBefore;
    try {
        netCDF::NcFile ncFile(fileNameBefore, netCDF::NcFile::read);

        // Read the time axis
        netCDF::NcDim timeDim = ncFile.getDim(ncTimeName);
        // Read the time variable
        netCDF::NcVar timeVar = ncFile.getVar(ncTimeName);
        // Get the time axis as a vector
        std::vector<double> timeVec(timeDim.getSize());
        timeVar.getVar(timeVec.data());

        // Time units nonsense
        std::string unitStr, timeUnit, sinceKeyWord, timeOrigin;
        timeVar.getAtt("units").getValues(unitStr);
        std::stringstream ss(unitStr);
        ss >> timeUnit >> sinceKeyWord >> timeOrigin;

        // Get a TimePoint for the origin
        const auto timePointOrigin = TimePoint(timeOrigin);

        // Multiply the time axis to get seconds
        double multiplier;
        if (timeUnit == "seconds")
            multiplier = 1.;
        else if (timeUnit == "minutes")
            multiplier = 60.;
        else if (timeUnit == "hours")
            multiplier = 3600.;
        else if (timeUnit == "days")
            multiplier = 24. * 3600.;
        else
            throw std::runtime_error("ParaGridInputs::readRawForcing(): unsupported time unit "
                + timeUnit + " in '" + unitStr + "'.\n");

        std::transform(timeVec.begin(), timeVec.end(), timeVec.begin(),
            [multiplier](const double t) { return t * multiplier; });

        // Because currentTime can't be captured in the lambda function
        const TimePoint& time = currentTime;
        // Get the index of the first TimePoint greater than the target.
        targetTIndexAfter = std::find_if(timeVec.begin(), timeVec.end(),
                                [time, timePointOrigin](const double t) {
                                    return TimePoint(timePointOrigin, Duration(t)) > time;
                                })
            - timeVec.begin();

        targetTIndexBefore = targetTIndexAfter - 1;
        timeRange.before = TimePoint(timePointOrigin, Duration(timeVec[targetTIndexBefore]));
        timeRange.after = TimePoint(timePointOrigin, Duration(timeVec[targetTIndexAfter]));

        /* We need to check if targetTIndexAfter is actually pointing at the right time, or if we
         * need to go on to the next file */
        if (TimePoint(timePointOrigin, Duration(timeVec[targetTIndexAfter])) < currentTime) {
            targetTIndexBefore = targetTIndexAfter;
            timeRange.before = TimePoint(timePointOrigin, Duration(timeVec[targetTIndexBefore]));
            timeRange.after = TimePoint(timePointOrigin, Duration(timeVec[targetTIndexAfter]))
                + Duration(timeVec[1] - timeVec[0]);
            targetTIndexAfter = 0;
            fileNameAfter = timeRange.after.format(filePath);
        }

        if (targetTIndexAfter < 0 || targetTIndexBefore < 0 || targetTIndexAfter >= timeVec.size()
            || targetTIndexBefore >= timeVec.size())
            throw std::out_of_range(
                "ParaGridInputs::readRawForcing::Target time index is out of range "
                "- how could this happen?\n");

        ncFile.close();
    } catch (const netCDF::exceptions::NcException& nce) {
        std::string ncWhat(nce.what());
        ncWhat += ": " + fileNameBefore;
        throw std::runtime_error(ncWhat);
    }

    rawDataBefore = readRawData(fileNameBefore, targetTIndexBefore);
    rawDataAfter = readRawData(fileNameAfter, targetTIndexAfter);
}

std::map<std::string, std::vector<FloatType>> ParaGridInputs::readRawData(
    const std::string& fileName, const size_t timeIndex = 0) const
{
    std::map<std::string, std::vector<FloatType>> data;
    try {
        netCDF::NcFile ncFile(fileName, netCDF::NcFile::read);

        auto availableForcings = ncFile.getVars();
        for (const std::string& varName : forcings) {
            // Don't try to read non-existent data
            if (!availableForcings.count(varName)) {
                continue;
            }
            netCDF::NcVar var = ncFile.getVar(varName);

            std::vector<netCDF::NcDim> dims = var.getDims();
            std::vector<size_t> start(dims.size(), 0);
            std::vector<size_t> count(dims.size());
            for (int i = 0; i < count.size(); ++i)
                count[i] = dims[i].getSize();

            // Pick the time slice if we have a time axis
            for (int i = 0; i < dims.size(); ++i) {
                if (dims[0].getName() == ncTimeName) {
                    start[0] = timeIndex;
                    count[0] = 1;
                }
            }

            data[varName] = std::vector<FloatType>(
                std::accumulate(count.begin(), count.end(), 1, std::multiplies<size_t>()));
            readNetCDFVar(var, start, count, data[varName].data());
        }
        ncFile.close();
    } catch (const netCDF::exceptions::NcException& nce) {
        std::string ncWhat(nce.what());
        ncWhat += ": " + filePath;
        throw std::runtime_error(ncWhat);
    }
    return data;
}

}