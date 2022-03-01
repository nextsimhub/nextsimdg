/*!
 * @file ElementData.cpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ElementData.hpp"
#include "include/Module.hpp"

namespace Nextsim {
ElementData::ElementData()
    : ElementData(1)
{
}

ElementData::ElementData(int nIceLayers)
    : PrognosticElementData(nIceLayers)
    , PhysicsData(nIceLayers)
{
    m_physicsImplData = std::move(Module::getInstance<IPhysics1d>());
}
//! Copy constructor
ElementData::ElementData(const ElementData& src)
{
    *this = static_cast<PrognosticElementData>(src);
    *this = static_cast<ExternalData>(src);
    *this = static_cast<PhysicsData>(src);
    this->m_physicsImplData = std::move(Module::getInstance<IPhysics1d>());
    *(this->m_physicsImplData) = *(src.m_physicsImplData);
}

//! Move constructor
ElementData::ElementData(ElementData&& src)
{
    *this = static_cast<PrognosticElementData>(src);
    *this = static_cast<ExternalData>(src);
    *this = static_cast<PhysicsData>(src);
    this->m_physicsImplData = std::move(src.m_physicsImplData);
}

//! Copy assignment operator
ElementData& ElementData::operator=(ElementData other)
{
    if (this == &other)
        return *this;

    *this = static_cast<PrognosticElementData>(other);
    *this = static_cast<ExternalData>(other);
    *this = static_cast<PhysicsData>(other);
    this->m_physicsImplData = std::move(Module::getInstance<IPhysics1d>());
    *(this->m_physicsImplData) = *(other.m_physicsImplData);

    return *this;
}

//! Move assignment operator
ElementData& ElementData::operator=(ElementData&& other)
{
    if (this == &other)
        return *this;

    *this = static_cast<PrognosticElementData>(other);
    *this = static_cast<ExternalData>(other);
    *this = static_cast<PhysicsData>(other);
    this->m_physicsImplData = std::move(other.m_physicsImplData);

    return *this;
}

//! Configures the PrognosticElementData and physics implementation aspects
//! of the object.
void ElementData::configure()
{
    PrognosticElementData::configure();
    Nextsim::tryConfigure(&Module::getImplementation<IPhysics1d>());
}

void ElementData::updateDerivedData(
    const PrognosticElementData& prog, const ExternalData& exter, PhysicsData& phys)
{
    m_physicsImplData->updateDerivedData(prog, exter, phys);
}

void ElementData::calculate(
    const PrognosticElementData& prog, const ExternalData& exter, PhysicsData& phys)
{
    m_physicsImplData->calculate(prog, exter, phys);
}

} /* namespace Nextsim */
