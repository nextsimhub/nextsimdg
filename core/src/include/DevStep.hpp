/*!
 * @file DevStep.hpp
 *
 * @date Jan 12, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_DEVSTEP_HPP
#define CORE_SRC_INCLUDE_DEVSTEP_HPP

#include "include/IModelStep.hpp"
#include "include/IStructure.hpp"

#include <string>

namespace Nextsim {

class DevStep : public IModelStep {
public:
    DevStep() = default;
    virtual ~DevStep() = default;

    // Member functions inherited from IModelStep
    void writeRestartFile(const std::string& filePath) override {};

    void setInitialData(IStructure& dataStructure) override { pStructure = &dataStructure; };

    // Member functions inherited from Iterant
    void init() override {};
    void start(const Iterator::TimePoint& startTime) override {};
    void iterate(const Iterator::Duration& dt) override;
    void stop(const Iterator::TimePoint& stopTime) override {};

private:
    IStructure* pStructure;
};

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_DEVSTEP_HPP */
