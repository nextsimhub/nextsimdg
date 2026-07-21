/*!
 * @author  Einar Olason <einar.olason@nersc.no>
 */

#ifndef NEXTSIM_DG_PARAGRIDINPUTS_HPP
#define NEXTSIM_DG_PARAGRIDINPUTS_HPP

#include "ModelArray.hpp"
#include "ModelState.hpp"
#include "Time.hpp"

#include <set>
#include <string>

namespace Nextsim {

class ParaGridInputs {

public:
    ParaGridInputs() = default;

    void configure(const std::string& filePathIn, const std::string& ncLonNameIn,
        const std::string& ncLatNameIn, const std::string& ncTimeNameIn,
        const std::set<std::string>& forcingsIn,
        const std::set<std::pair<std::string, std::string>>& vectorsIn,
        const ModelArray& modelLonsIn, const ModelArray& modelLatsIn);

    void update(const TimePoint& time);
    [[nodiscard]] ModelArray getField(const std::string& fieldName);
    [[nodiscard]] bool isValid(const std::string& fieldName) const
    {
        return forcingStateBefore.data.count(fieldName) && forcingStateAfter.data.count(fieldName);
    }

private:
    [[nodiscard]] ModelState interpolateSpatially(
        const std::map<std::string, std::vector<FloatType>>& rawData);
    void rotateInputVectors(std::map<std::string, std::vector<FloatType>>& rawData);

    std::string filePath, ncTimeName, ncLonName, ncLatName;
    std::set<std::string> forcings;
    std::set<std::pair<std::string, std::string>> vectors;
    ModelArray modelLons, modelLats;

    struct {
        TimePoint before = std::numeric_limits<TimePoint>::min(),
                  after = std::numeric_limits<TimePoint>::min();
    } timeRange;

    TimePoint currentTime;

    ModelState forcingStateBefore, forcingStateAfter;

    void vectorRotationLogic(const std::vector<FloatType>& vectorIn1st,
        const std::vector<FloatType>& vectorIn2nd, std::vector<FloatType>& vectorOut1st,
        std::vector<FloatType>& vectorOut2nd);

    void readRawForcing(std::map<std::string, std::vector<FloatType>>& rawDataBefore,
        std::map<std::string, std::vector<FloatType>>& rawDataAfter);

    [[nodiscard]] std::map<std::string, std::vector<FloatType>> readRawData(
        const std::string& fileName, size_t timeIndex) const;
};

}

#endif // NEXTSIM_DG_PARAGRIDINPUTS_HPP
