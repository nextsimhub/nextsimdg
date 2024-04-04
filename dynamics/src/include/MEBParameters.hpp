/*!
 * @file MEBParameters.hpp
 * @date 18 Jan 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#ifndef __MEBPARAMETERS_HPP
#define __MEBPARAMETERS_HPP

#include "DynamicsParameters.hpp"

namespace Nextsim {

class MEBParameters : public DynamicsParameters {
public:
    // MEB
    double compaction_param; //!< Compation parameter
    double nu0; //!< \param Poisson's ratio
    double young; //!< \param Young's modulus
    double P0; //! < \param Constant to parametrize Pmax
    double damage_timescale = 1.; //<! Damage timescale
    double undamaged_time_relaxation_sigma = 1e5; //!< \param lambda
    int exponent_relaxation_sigma;
    double c0; //! \param
    double sigma_c0; //! \param

    double time_relaxation_damage; //!< 25 days in seconds
    double compression_factor; //! \param Max pressure for damaged converging ice
    double exponent_compression_factor; //! \param Power of ice thickness in the pressure coefficient

    // TODO missing 45\deg it goes to Compresssion
    double tan_phi; //!< \param tan_phi (double) Internal friction coefficient (mu)
    double sin_phi; //!< \param sin_phi (double) Internal friction coefficient (mu)

    double compr_strength; //! \param compr_strength (double) Maximum compressive strength [N/m2]
    double C_lab; //! \param C_lab (double) Test [Pa]

    MEBParameters()
    {
        // MEB
        compaction_param = -20.; //!< Compation parameter
        nu0 = 1. / 3.; //!< \param Poisson's ratio
        // young = 1e9; //!< \param Young's modulus Compression
        young = 5.96e8; //!< \param Young's modulus
        P0 = 10.e3; //! < \param Constant to parametrize Pmax

        undamaged_time_relaxation_sigma = 1e7; //!< \param lambda
        exponent_relaxation_sigma = 5;
        c0 = 10.e3; //! \param
        sigma_c0 = 50.e3; //! \param

        damage_timescale = 1.; //<! Damage timescale
        time_relaxation_damage = 2160000.; //!< 25 days in seconds
        compression_factor = 10e3; //! \param Max pressure for damaged converging ice
        exponent_compression_factor = 1.5; //! \param Power of ice thickness in the pressure coefficient
        tan_phi = 0.7; //!< \param tan_phi (double) Internal friction coefficient (mu)
        sin_phi = 0.7; //!< \param sin_phi (double) Internal friction coefficient (mu)

        compr_strength = 1e10; //! \param compr_strength (double) Maximum compressive strength [N/m2]
        C_lab = 2.0e6; //! \param C_lab (double) Test [Pa]
    }
};

} /* namespace Nextsim */

#endif /* __MEB_HPP */
