/*!
 * @file IModelStep.hpp
 *
 * @date Jan 18, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_IMODELSTEP_HPP
#define CORE_SRC_INCLUDE_IMODELSTEP_HPP

#include "include/IStructure.hpp"
#include "include/Iterator.hpp"

namespace Nextsim {

class IModelStep : public Iterator::Iterant {
public:
    IModelStep() = default;
    virtual ~IModelStep() = default;

    void setInitFile(const std::string& filePath) { initialRestartFilePath = filePath; };
    virtual void writeRestartFile(const std::string& filePath) = 0;

    virtual void setInitialData(IStructure& dataStructure) = 0;

    // Member functions inherited from Iterant
    virtual void init() = 0;
    virtual void start(const Iterator::TimePoint& startTime) = 0;
    virtual void iterate(const Iterator::Duration& dt) = 0;
    virtual void stop(const Iterator::TimePoint& stopTime) = 0;

protected:
    std::string initialRestartFilePath;
};

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_IMODELSTEP_HPP */
