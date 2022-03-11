/*!
 * @file ElementData.hpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_ELEMENTDATA_HPP
#define SRC_INCLUDE_ELEMENTDATA_HPP

#include "Configured.hpp"
#include "ExternalData.hpp"
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
    ElementData();
    ElementData(int nIceLayers);

    //! Copy constructor
    ElementData(const ElementData& src);

    //! Move constructor
    ElementData(ElementData&& src);

    ~ElementData() = default;

    using PrognosticData::operator=;
    using PhysicsData::operator=;
    using ExternalData::operator=;

    //! Copy assignment operator
    ElementData& operator=(ElementData other);

    //! Move assignment operator
    ElementData& operator=(ElementData&& other);

    //! Configures the PrognosticData and physics implementation aspects of the
    //!  object.
    void configure() override;

    void updateDerivedData(
        const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys);

    void calculate(const PrognosticData& prog, const ExternalData& exter, PhysicsData& phys);

private:
    std::unique_ptr<IPhysics1d> m_physicsImplData;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_ELEMENTDATA_HPP */
