/*!
 * @file mMEB.hpp
 * @date 1 Jun 2022
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#ifndef __MMEB_HPP
#define __MMEB_HPP

#include "MEB.hpp"
#include "codeGenerationDGinGauss.hpp"
#include "dgVector.hpp"

namespace Nextsim {

/*!
 * This namespace collects the routines required for the MEB solver
 */
namespace mMEB {

    inline constexpr double SQR(double x) { return x * x; }

    //! WIP new discretization scheme
    template <int DGstress, int DGtracer>
    void mMEBStressRelaxation(const Mesh& mesh, CellVector<DGstress>& S11, CellVector<DGstress>& S12,
        CellVector<DGstress>& S22, CellVector<DGtracer>& D, CellVector<DGtracer>& A, const double dt_momentum)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {

            double const expC = std::exp(RefScale::compaction_param * (1. - A(i, 0)));
            double const elasticity = RefScale::young * (1. - D(i, 0)) * expC;

            double const sigma_s = std::hypot((S11(i, 0) - S22(i, 0)) / 2., S12(i, 0));
            // update sigma_n
            double const sigma_n = 0.5 * (S11(i, 0) + S22(i, 0));
            double const C_fix = RefScale::C_lab * std::sqrt(0.1 / mesh.hx);

            // d critical Eqn. (29)
            double dcrit;
            if (sigma_n < -RefScale::compr_strength)
                dcrit = -RefScale::compr_strength / sigma_n;
            else
                dcrit = C_fix / (sigma_s + RefScale::tan_phi * sigma_n);

            double const td = mesh.hx * std::sqrt(2. * (1. + RefScale::nu0) * RefScale::rho_ice)
                / std::sqrt(elasticity);

            // Calculate the adjusted level of damage
            if ((0. < dcrit) && (dcrit < 1.)) {
                // Eqn. (34)
                D(i, 0) += (1.0 - D(i, 0)) * (1.0 - dcrit) * dt_momentum / td;

                S11.row(i) -= S11.row(i) * (1. - dcrit) * dt_momentum / td;
                S12.row(i) -= S12.row(i) * (1. - dcrit) * dt_momentum / td;
                S22.row(i) -= S22.row(i) * (1. - dcrit) * dt_momentum / td;
            }
            // Relax damage
            D(i, 0) = std::max(0., D(i, 0) - dt_momentum / RefScale::time_relaxation_damage);
        }
    }

    //! WIP new discretization scheme
    template <int DGstress, int DGtracer>
    void mMEBDamageUpdate(const Mesh& mesh, CellVector<DGstress>& S11, CellVector<DGstress>& S12,
        CellVector<DGstress>& S22, CellVector<DGtracer>& D, CellVector<DGtracer>& A, const double dt_momentum)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {

            double const expC = std::exp(RefScale::compaction_param * (1. - A(i, 0)));
            double const elasticity = RefScale::young * (1. - D(i, 0)) * expC;

            double const sigma_s = std::hypot((S11(i, 0) - S22(i, 0)) / 2., S12(i, 0));
            // update sigma_n
            double const sigma_n = 0.5 * (S11(i, 0) + S22(i, 0));
            double const C_fix = RefScale::C_lab * std::sqrt(0.1 / mesh.hx);

            // d critical Eqn. (29)
            double dcrit;
            if (sigma_n < -RefScale::compr_strength)
                dcrit = -RefScale::compr_strength / sigma_n;
            else
                dcrit = C_fix / (sigma_s + RefScale::tan_phi * sigma_n);

            double const td = mesh.hx * std::sqrt(2. * (1. + RefScale::nu0) * RefScale::rho_ice)
                / std::sqrt(elasticity);
            // Calculate the adjusted level of damage
            if ((0. < dcrit) && (dcrit < 1.)) {
                // Eqn. (34)
                D(i, 0) += (1.0 - D(i, 0)) * (1.0 - dcrit) * dt_momentum / td;

                //S11.row(i) -= S11.row(i) * (1. - dcrit) * dt_momentum / td;
                //S12.row(i) -= S12.row(i) * (1. - dcrit) * dt_momentum / td;
                //S22.row(i) -= S22.row(i) * (1. - dcrit) * dt_momentum / td;
            }
            // Relax damage
            //D(i, 0) = std::max(0., D(i, 0) - dt_momentum / RefScale::time_relaxation_damage);
        }
    }

    //! WIP new discretization scheme
    template <int DGstress, int DGtracer>
    void mMEBStressUpdate(const Mesh& mesh, CellVector<DGstress>& S11, CellVector<DGstress>& S12,
        CellVector<DGstress>& S22, const CellVector<DGstress>& E11, const CellVector<DGstress>& E12,
        const CellVector<DGstress>& E22, const CellVector<DGtracer>& H,
        const CellVector<DGtracer>& A, CellVector<DGtracer>& D, const double dt_adv,
        const double alpha, CellVector<DGstress>& S11_mmeb, CellVector<DGstress>& S12_mmeb,
        CellVector<DGstress>& S22_mmeb)
    {

        //! Stress Update
#pragma omp parallel for
        for (size_t i = 0; i < mesh.n; ++i) {

            // There's no ice so we set sigma to 0 and carry on
            const double min_c = 0.1;
            if (A(i, 0) <= min_c) {
                D(i, 0) = 0.0;
                S11.row(i) *= 0.0;
                S12.row(i) *= 0.0;
                S22.row(i) *= 0.0;
                continue;
            }

            //! - Updates the internal stress
            double sigma_n = 0.5 * (S11(i, 0) + S22(i, 0));
            double const expC = std::exp(RefScale::compaction_param * (1. - A(i, 0)));
            double const time_viscous = RefScale::undamaged_time_relaxation_sigma
                * std::pow((1. - D(i, 0)) * expC, RefScale::exponent_relaxation_sigma - 1.);

            //! Plastic failure tildeP
            double tildeP;
            if (sigma_n < 0.) {
                double const Pmax = RefScale::compression_factor
                    * pow(H(i, 0), RefScale::exponent_compression_factor)
                    * exp(RefScale::compaction_param * (1.0 - A(i, 0)));
                tildeP = std::min(1., -Pmax / sigma_n);
            } else {
                tildeP = 0.;
            }

            // \lambda / (\lambda + dt*(1.+tildeP)) Eqn. 32
            //double const multiplicator
            //    = std::min(1. - 1e-12, time_viscous / (time_viscous + dt_adv * (1. - tildeP)));

            //double const multiplicator
            //    = std::min(1. - 1e-12, 1. / (1. + alpha + dt_adv / time_viscous * (1. - tildeP)));

            double const elasticity = RefScale::young * (1. - D(i, 0)) * expC;

            double const Dunit_factor = 1. / (1. - (RefScale::nu0 * RefScale::nu0));

            double X = (1. - tildeP) / time_viscous;
            //double Z = (alpha + 1) / dt_adv + X; //alpha + 1. / dt_adv + X; //
            //double Y = alpha + 1. + dt_adv * X;

            /*
            S11.row(i) *= alpha / Y;
            S12.row(i) *= alpha / Y;
            S22.row(i) *= alpha / Y;

            S11.row(i) += 1. / Y * S11_mmeb.row(i) + 1. / Y * dt_adv * elasticity * (1 / (1 + RefScale::nu0) * E11.row(i) + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i)));
            S12.row(i) += 1. / Y * S12_mmeb.row(i) + 1. / Y * dt_adv * elasticity * 1 / (1 + RefScale::nu0) * E12.row(i);
            S22.row(i) += 1. / Y * S22_mmeb.row(i) + 1. / Y * dt_adv * elasticity * (1 / (1 + RefScale::nu0) * E22.row(i) + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i)));
            */

            // Alternative
            //S11.row(i) = 1. / Z * (1. / dt_adv * (alpha * S11.row(i) + S11_mmeb.row(i)) + elasticity * (1 / (1 + RefScale::nu0) * E11.row(i) + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i))));
            //S12.row(i) = 1. / Z * (1. / dt_adv * (alpha * S12.row(i) + S12_mmeb.row(i)) + elasticity * 1 / (1 + RefScale::nu0) * E12.row(i));
            //S22.row(i) = 1. / Z * (1. / dt_adv * (alpha * S22.row(i) + S22_mmeb.row(i)) + elasticity * (1 / (1 + RefScale::nu0) * E22.row(i) + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i))));

            // Step by step
            /*
            S11.row(i) *= (1.0 - 1.0 / alpha);
            S12.row(i) *= (1.0 - 1.0 / alpha);
            S22.row(i) *= (1.0 - 1.0 / alpha);

            S11.row(i) += 1. / alpha * elasticity * (1 / (1 + RefScale::nu0) * E11.row(i) + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i)));
            S12.row(i) += 1. / alpha * elasticity * 1 / (1 + RefScale::nu0) * E12.row(i);
            S22.row(i) += 1. / alpha * elasticity * (1 / (1 + RefScale::nu0) * E22.row(i) + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i)));
            */

            //std::cout << 1 / dt_adv << std::endl;std::exit(0);

            // Step by step
            S11.row(i) *= (1. - (X - 1. / dt_adv) / (alpha + 1. / dt_adv));
            S12.row(i) *= (1. - (X - 1. / dt_adv) / (alpha + 1. / dt_adv));
            S22.row(i) *= (1. - (X - 1. / dt_adv) / (alpha + 1. / dt_adv));

            //add eplitic parts
            S11.row(i) += 1. / (alpha + 1. / dt_adv) * (1. / dt_adv * S11_mmeb.row(i));
            S12.row(i) += 1. / (alpha + 1. / dt_adv) * (1. / dt_adv * S12_mmeb.row(i));
            S22.row(i) += 1. / (alpha + 1. / dt_adv) * (1. / dt_adv * S22_mmeb.row(i));
            //add elasticity
            S11.row(i) += 1. / (alpha + 1. / dt_adv) * elasticity * (1. / (1. + RefScale::nu0) * E11.row(i) + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i)));
            S12.row(i) += 1. / (alpha + 1. / dt_adv) * elasticity * 1. / (1. + RefScale::nu0) * E12.row(i);
            S22.row(i) += 1. / (alpha + 1. / dt_adv) * elasticity * (1. / (1. + RefScale::nu0) * E22.row(i) + Dunit_factor * RefScale::nu0 * (E11.row(i) + E22.row(i)));
        }
    }

} /* namespace mMEB */

} /* namespace Nextsim */

#endif /* __MMEB_HPP */
