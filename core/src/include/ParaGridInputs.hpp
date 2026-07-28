/*!
 * @author  Einar Olason <einar.olason@nersc.no>
 */

#ifndef NEXTSIM_DG_PARAGRIDINPUTS_HPP
#define NEXTSIM_DG_PARAGRIDINPUTS_HPP

#include "ModelArray.hpp"
#include "ModelState.hpp"
#include "Time.hpp"

#include <regex>
#include <set>
#include <string>

namespace Nextsim {

class ParaGridInputs {

public:
    ParaGridInputs() = default;

    void setData(const TimePoint& time, const std::string& filePathIn,
        const std::string& ncLonNameIn, const std::string& ncLatNameIn,
        const std::string& ncTimeNameIn, const std::set<std::string>& forcingsIn,
        const std::set<std::pair<std::string, std::string>>& vectorsIn,
        const ModelArray& modelLonsIn, const ModelArray& modelLatsIn);

    void update(const TimePoint& time);
    [[nodiscard]] ModelArray getField(const std::string& fieldName);
    [[nodiscard]] bool isValid(const std::string& fieldName) const
    {
        return forcingStateBefore.data.count(fieldName) && forcingStateAfter.data.count(fieldName);
    }

private:
    typedef struct {
        std::vector<size_t> dims;
        std::map<std::string, std::vector<FloatType>> data;
    } RawDataMap;

    void setWeights();
    void setWeights1D();
    void setWeights2D();
    bool recursiveBisectSearch(size_t k, FloatType targetLon, FloatType targetLat, size_t i,
        size_t ii, size_t j, size_t jj);
    void orthographicProjection(FloatType lon, FloatType lat, FloatType lon0, FloatType lat0,
        FloatType& x, FloatType& y) const;
    [[nodiscard]] bool pointInBoundingBox(
        const std::vector<FloatType>& xCorners, const std::vector<FloatType>& yCorners);
    [[nodiscard]] bool findLocalCoordinates(size_t k, FloatType x00, FloatType y00, FloatType x10,
        FloatType y10, FloatType x01, FloatType y01, FloatType x11, FloatType y11);
    [[nodiscard]] ModelState interpolateSpatially(const RawDataMap& rawData);
    void rotateInputVectors(RawDataMap& rawData);

    std::string filePath, ncTimeName, ncLonName, ncLatName;
    ModelArray modelLons, modelLats;
    RawDataMap forcingLonLats;
    std::set<std::string> forcings;
    std::set<std::pair<std::string, std::string>> vectors;
    std::vector<size_t> ij00, ij01, ij10, ij11;
    ModelArray xi, eta;

    struct {
        TimePoint before = std::numeric_limits<TimePoint>::min(),
                  after = std::numeric_limits<TimePoint>::min();
    } timeRange;

    TimePoint currentTime;

    ModelState forcingStateBefore, forcingStateAfter;

    void vectorRotationLogic(const std::vector<FloatType>& vectorIn1st,
        const std::vector<FloatType>& vectorIn2nd, std::vector<FloatType>& vectorOut1st,
        std::vector<FloatType>& vectorOut2nd);

    void readRawForcing(RawDataMap& rawDataBefore, RawDataMap& rawDataAfter);

    [[nodiscard]] RawDataMap readRawData(
        const TimePoint& time, const std::set<std::string>& fields, size_t timeIndex = 0) const;

    [[nodiscard]] double wrapLon(const double lon) const
    {
        double wrappedLon = std::fmod(lon + 180., 360.);
        if (wrappedLon < 0.0)
            wrappedLon += 360.;

        wrappedLon -= 180.;
        return wrappedLon;
    }

    [[nodiscard]] std::string formatFileName(const std::string& filePathIn, const TimePoint& time,
        const std::string& variable = "NONE") const
    {
        const std::regex pattern("%v");
        const std::string temp = std::regex_replace(filePathIn, pattern, variable);
        return time.format(temp);
    }
};

}

#endif // NEXTSIM_DG_PARAGRIDINPUTS_HPP
