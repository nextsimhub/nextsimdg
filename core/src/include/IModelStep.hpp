/*!
 * @file IModelStep.hpp
 *
 * @date Jan 18, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IMODELSTEP_HPP
#define IMODELSTEP_HPP

#include "include/Iterator.hpp"
#include "include/PrognosticData.hpp"

namespace Nextsim {

class IModelStep : public Iterator::Iterant {
public:
    IModelStep() = default;
    virtual ~IModelStep() = default;

    void setInitFile(const std::string& filePath) { initialRestartFilePath = filePath; };
    virtual void writeRestartFile(const std::string& filePath) = 0;

    virtual void setInitialData(PrognosticData& data) = 0;

    // Member functions inherited from Iterant
    virtual void init() = 0;
    virtual void start(const TimePoint& startTime) = 0;
    virtual void iterate(const TimestepTime& dt) = 0;
    virtual void stop(const TimePoint& stopTime) = 0;

protected:
    std::string initialRestartFilePath;
};

} /* namespace Nextsim */

#endif /* IMODELSTEP_HPP */
