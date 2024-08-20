/*!
 * @file OASISCoupledOcean.hpp
 *
 * @date Sep 26, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef OASISCOUPLEDOCEAN_HPP
#define OASISCOUPLEDOCEAN_HPP

#include "include/IFreezingPoint.hpp"
#include "include/IOceanBoundary.hpp"
#include "include/Module.hpp"
#include "include/OASISCoupled.hpp"
#include "include/SlabOcean.hpp"

namespace Nextsim {

//* Ocean boundary data values that are hardcoded.
class OASISCoupledOcean : public IOceanBoundary, OASISCoupled {
public:
    OASISCoupledOcean()
        : IOceanBoundary()
    {
    }
    ~OASISCoupledOcean() { OASISCoupled::~OASISCoupled(); }

    std::string getName() const override { return "OASISCoupledOcean"; }
    void updateBefore(const TimestepTime& tst) override;
    void updateAfter(const TimestepTime& tst) override;
    void setMetadata(const ModelMetadata& metadata) override;

private:
    int bundleSize = 1;

    SlabOcean slabOcean;

    void updateTf(size_t i, const TimestepTime& tst)
    {
        tf[i] = Module::getImplementation<IFreezingPoint>()(sss[i]);
    }

    enum couplingIdIn { SST, SSS, UOCEAN, VOCEAN, SSH, MLD }; // MLD is optional
    enum couplingIdOut { TAUX, TAUY, TAUMOD, EMP, QNOSUN, QSW, SFLX, CICE };

    const std::vector<std::string> cplStringsIn
        = { "I_SST", "I_SSS", "I_Uocn", "I_Vocn", "I_SSH", "IFrcQsr" };
    const std::vector<std::string> cplStringsOut
        = { "I_taux", "I_tauy", "I_taumod", "I_fwflux", "I_rsnos", "I_rsso", "I_sfi", "I_conc" };

    std::vector<int> cplIdIn, cplIdOut;
};

} /* namespace Nextsim */

#endif /* OASISCOUPLEDOCEAN_HPP */
