/*!
 * @file DevStep.hpp
 *
 * @date Jan 12, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef DEVSTEP_HPP
#define DEVSTEP_HPP

#include "include/IModelStep.hpp"
#include "include/IStructure.hpp"

#include <string>

namespace Nextsim {

//! A class providing a simple implementation of Iterator.
class DevStep : public IModelStep {
public:
    DevStep();
    virtual ~DevStep() = default;

    // Member functions inherited from IModelStep
    void writeRestartFile(const std::string& filePath) override {};

    void setData(PrognosticData& pDat) override { pData = &pDat; }
    void setMetadata(ModelMetadata& metadata) override { mData = &metadata; }

    /*!
     * Sets the period with which restart files are created.
     *
     * @param restartPeriod The Duration between writing out restart files.
     * @param filename The file name pattern to be used for the restart files.
     *        The string will be used as a format for the time string, based on
     *        the syntax of strftime.
     */
    void setRestartDetails(const Duration& restartPeriod, const std::string& fileName);
    // Member functions inherited from Iterant
    void init() override;
    void start(const TimePoint& startTime) override;
    void iterate(const TimestepTime& dt) override;
    void stop(const TimePoint& stopTime) override {};

private:
    PrognosticData* pData;
    ModelMetadata* mData;
    Duration m_restartPeriod;
    TimePoint lastOutput;
    std::string m_restartFileName;
};

} /* namespace Nextsim */

#endif /* DEVSTEP_HPP */
