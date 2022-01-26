/*!
 * @file IPrognosticUpdater.hpp
 *
 * @date Jan 21, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_IPROGNOSTICUPDATER_HPP
#define CORE_SRC_INCLUDE_IPROGNOSTICUPDATER_HPP

namespace {
class IPrognosticUpdater {
public:
    virtual ~IPrognosticUpdater() = default;

    virtual double updatedIceThickness() const = 0;
    virtual double updatedSnowThickness() const = 0;
    virtual double updatedIceConcentration() const = 0;
    virtual const std::vector<double>& updatedIceTemperatures() const = 0;
};

}

#endif /* CORE_SRC_INCLUDE_IPROGNOSTICUPDATER_HPP */
