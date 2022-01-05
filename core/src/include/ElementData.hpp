/*!
 * @file ElementData.hpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_ELEMENTDATA_HPP
#define SRC_INCLUDE_ELEMENTDATA_HPP

#include "ExternalData.hpp"
#include "ModuleLoader.hpp"
#include "PrognosticData.hpp"
#include "include/IPhysics1d.hpp"
#include "include/PhysicsData.hpp"

namespace Nextsim {

//! A class to be used when a non-specific class derived from BaseElementData
//! is needed.
class UnusedData : public BaseElementData {
};

/*!
 * @brief The class which holds all the data for a single element of the model.
 *
 * @details Inherits from PrognosticData, PhysicsData and ExternalData. The
 * physics implementation is provided as a per-instance module.
 */
class ElementData : public PrognosticData,
                    public PhysicsData,
                    public ExternalData,
                    public UnusedData,
                    public Configured<ElementData> {
public:
    ElementData()
    {
        m_physicsImplData = std::move(ModuleLoader::getLoader().getInstance<IPhysics1d>());
    }
    ~ElementData() = default;

    using PrognosticData::operator=;
    using PhysicsData::operator=;
    using ExternalData::operator=;

    //! Configures the PrognosticData and physics implementation aspects of the
    //!  object.
    void configure() override { PrognosticData::configure(); }

    void updateDerivedData(
        const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys)
    {
        m_physicsImplData->updateDerivedData(prog, exter, phys);
    }

    void calculate(
        const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys)
    {
        m_physicsImplData->calculate(prog, exter, phys);
    }

private:
    std::unique_ptr<IPhysics1d> m_physicsImplData;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_ELEMENTDATA_HPP */
