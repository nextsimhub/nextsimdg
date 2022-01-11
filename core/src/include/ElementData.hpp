/*!
 * @file ElementData.hpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_ELEMENTDATA_HPP
#define SRC_INCLUDE_ELEMENTDATA_HPP

#include "Configured.hpp"
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
    //! Copy constructor
    ElementData(const ElementData& src)
    {
        *this = static_cast<PrognosticData>(src);
        *this = static_cast<ExternalData>(src);
        *this = static_cast<PhysicsData>(src);
        this->m_physicsImplData = std::move(ModuleLoader::getLoader().getInstance<IPhysics1d>());
        *(this->m_physicsImplData) = *(src.m_physicsImplData);
    }
    //! Move constructor
    ElementData(ElementData&& src)
    {
        *this = static_cast<PrognosticData>(src);
        *this = static_cast<ExternalData>(src);
        *this = static_cast<PhysicsData>(src);
        this->m_physicsImplData = std::move(src.m_physicsImplData);
    }
    ~ElementData() = default;

    using PrognosticData::operator=;
    using PhysicsData::operator=;
    using ExternalData::operator=;

    //! Copy assignment operator
    ElementData& operator=(ElementData other)
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
    ElementData& operator=(ElementData&& other)
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
    void configure() override
    {
        PrognosticData::configure();
        Nextsim::tryConfigure(&ModuleLoader::getLoader().getImplementation<IPhysics1d>());
    }

    void updateDerivedData(const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys)
    {
        m_physicsImplData->updateDerivedData(prog, exter, phys);
    }

    void calculate(const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys)
    {
        m_physicsImplData->calculate(prog, exter, phys);
    }

private:
    std::unique_ptr<IPhysics1d> m_physicsImplData;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_ELEMENTDATA_HPP */
