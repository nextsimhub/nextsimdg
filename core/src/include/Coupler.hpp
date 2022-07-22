/*!
 * @file Coupler.hpp
 *
 * @date Jul 22, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef COUPLER_HPP
#define COUPLER_HPP

#include "include/Configured.hpp"
#include "include/ModelMetadata.hpp"

namespace Nextsim {

class CouplerData;

class Coupler : public Configured<Coupler> {
public:
    //! Performs the configuration of the coupler based on configuration files
    //! and arguments.
    void configure() override;
    /*!
     * @brief Performs the initialization based on the metadata derived from
     * the restart file.
     * @param meta the ModelMetadata object describing the current model.
     */
    void initialize(const ModelMetadata& meta);
    //! Perform any shutdown tasks required by the coupling system.
    void terminate(); // May need to include metadata argument?
protected:
    CouplerData* pData;
};

} /* namespace Nextsim */

#endif /* COUPLER_HPP */
