/*!
 * @file HiblerConcentration.hpp
 *
 * @date Nov 11, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_HIBLERCONCENTRATION_HPP
#define SRC_INCLUDE_HIBLERCONCENTRATION_HPP

#include "include/Configured.hpp"
#include "IConcentrationModel.hpp"

namespace Nextsim {

//! The implementation class of Hibler's model of ice concentration.
class HiblerConcentration : public IConcentrationModel, public Configured<HiblerConcentration> {
public:
    HiblerConcentration() = default;
    virtual ~HiblerConcentration() = default;

    //! Configure the parameters of the model.
    void configure() override;
    enum {
        H0_KEY,
        PHIM_KEY,
    };

    /*!
     * @brief Calculates the amount of freezing during the timestep from the
     * Hibler model.
     *
     * @param prog PrognosticData for this element (constant).
     * @param phys PhysicsData for this element.
     * @param nsphys Nextsim physics implementation data for this element.
     */
    double freeze(const PrognosticData&, PhysicsData&, NextsimPhysics&) const override;
    /*!
     * @brief Calculates the amount of melting during the timestep from the
     * Hibler model.
     *
     * @param prog PrognosticData for this element (constant).
     * @param phys PhysicsData for this element.
     * @param nsphys Nextsim physics implementation data for this element.
     */
    double melt(const PrognosticData&, PhysicsData&, NextsimPhysics&) const override;

    /*!
     * @brief Sets the value of the h0 parameter.
     *
     * @param h0_in The value of the h0 parameter to be set.
     */
    inline static void setH0(double h0_in) { h0 = h0_in; };

private:
    static double h0;
    static double phiM;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_HIBLERCONCENTRATION_HPP */
