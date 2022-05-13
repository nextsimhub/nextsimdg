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
    DevStep() = default;
    virtual ~DevStep() = default;

    // Member functions inherited from IModelStep
    void writeRestartFile(const std::string& filePath) override {};

    void setInitialData(PrognosticData& pDat) override { pData = &pDat; }

    // Member functions inherited from Iterant
    void init() override {};
    void start(const TimePoint& startTime) override {};
    void iterate(const TimestepTime& dt) override;
    void stop(const TimePoint& stopTime) override {};

private:
    PrognosticData* pData;
};

} /* namespace Nextsim */

#endif /* DEVSTEP_HPP */
