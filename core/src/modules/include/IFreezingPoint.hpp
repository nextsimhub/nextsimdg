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
template <typename Base> class FreezingPointImpl : public Base {
public:
    static_assert(
        std::is_base_of_v<IFreezingPoint, Base>, "Base needs to inherit from IFreezingPoint");

    FloatType operator()(FloatType sss) const override { return Base::calculate(sss); }
    void update(ModelArrayAuto& tf, const ConstModelArrayAuto& sss) const override
    {
        ModelComponent::overElementsAuto(
            OVER_ELEMENTS_CLASS_LAMBDA(const ElementIndex i) { tf[i] = Base::calculate(sss[i]); });
    }
};
}
#endif /* SRC_INCLUDE_IFREEZINGPOINT_HPP_ */
