/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef SRC_INCLUDE_IFREEZINGPOINT_HPP_
#define SRC_INCLUDE_IFREEZINGPOINT_HPP_

#include "include/KernelAlternatives.hpp"
#include "include/ModelArrayAccessor.hpp"
#include "include/ModelComponent.hpp"

namespace Nextsim {

//! The interface class for calculation of the freezing point of seawater.
class IFreezingPoint : public ModelComponent {
public:
    virtual ~IFreezingPoint() = default;

    void setData(const ModelState::DataMap& state) override { }

    /*!
     * @brief A virtual function that calculates the freezing point of
     * seawater.
     *
     * @details Freezing point in ˚C of water with the given salinity at
     * standard pressure.
     *
     * @param sss Sea surface salinity [PSU]
     */
    virtual FloatType operator()(FloatType sss) const = 0;

    /*!
     * @brief A virtual function that calculates the freezing point of
     * seawater for the whole field.
     *
     * @details Freezing point in ˚C of water with the given salinity at
     * standard pressure.
     *
     * @param tf Target for the freezing point [˚C].
     * @param sss Sea surface salinity [PSU]
     */
    virtual void update(ModelArrayAuto& tf, const ConstModelArrayAuto& sss) const = 0;
};

//! Helper to generate implementations of IFreezingPoint based on just a scalar function calculate.
// Base needs to inherit from IFreezingPoint and implement getName().
template <typename Base> class FreezingPointImpl : public IFreezingPoint {
public:
    std::string getName() const { return m_impl.getName(); }

    FloatType operator()(FloatType sss) const override { return m_impl.calculate(sss); }
    void update(ModelArrayAuto& tf, const ConstModelArrayAuto& sss) const override
    {
        // we can't capture this directly because copying polymorphic objects between host and device is undefined behaviour
        const Base& impl = m_impl;
        ModelComponent::overElementsAuto(
            OVER_ELEMENTS_LAMBDA(const ElementIndex i) { tf[i] = impl.calculate(sss[i]); });
    }

private:
    Base m_impl;
};
}
#endif /* SRC_INCLUDE_IFREEZINGPOINT_HPP_ */
