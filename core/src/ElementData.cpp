/*!
 * @file ElementData.cpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ElementData.hpp"

namespace Nextsim {
ElementData::ElementData()
{
    m_physicsImplData = std::move(ModuleLoader::getLoader().getInstance<IPhysics1d>());
}

//! Copy constructor
ElementData::ElementData(const ElementData& src)
{
    *this = static_cast<PrognosticData>(src);
    *this = static_cast<ExternalData>(src);
    *this = static_cast<PhysicsData>(src);
    this->m_physicsImplData = std::move(ModuleLoader::getLoader().getInstance<IPhysics1d>());
    *(this->m_physicsImplData) = *(src.m_physicsImplData);
}

//! Move constructor
ElementData::ElementData(ElementData&& src)
{
    *this = static_cast<PrognosticData>(src);
    *this = static_cast<ExternalData>(src);
    *this = static_cast<PhysicsData>(src);
    this->m_physicsImplData = std::move(src.m_physicsImplData);
}

//! Copy assignment operator
ElementData& ElementData::operator=(ElementData other)
{
    if (this == &other)
        return *this;

    *this = static_cast<PrognosticData>(other);
    *this = static_cast<ExternalData>(other);
    *this = static_cast<PhysicsData>(other);
    this->m_physicsImplData = std::move(ModuleLoader::getLoader().getInstance<IPhysics1d>());
    *(this->m_physicsImplData) = *(other.m_physicsImplData);

    return *this;
}

//! Move assignment operator
ElementData& ElementData::operator=(ElementData&& other)
{
    if (this == &other)
        return *this;

    *this = static_cast<PrognosticData>(other);
    *this = static_cast<ExternalData>(other);
    *this = static_cast<PhysicsData>(other);
    this->m_physicsImplData = std::move(other.m_physicsImplData);

    return *this;
}

//! Configures the PrognosticData and physics implementation aspects of the
//!  object.
void ElementData::configure()
{
    PrognosticData::configure();
    Nextsim::tryConfigure(&ModuleLoader::getLoader().getImplementation<IPhysics1d>());
}

void ElementData::updateDerivedData(
    const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys)
{
    m_physicsImplData->updateDerivedData(prog, exter, phys);
}

void ElementData::calculate(
    const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys)
{
    m_physicsImplData->calculate(prog, exter, phys);
}

} /* namespace Nextsim */
