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

    // A map to relate the strings in the namcouple file to the numbers def_var spits out
    std::map<std::string, int> couplingId;
    std::string SST = "I_SST";
    std::string SSS = "I_SSS";
    std::string UOCEAN = "I_Uocn";
    std::string VOCEAN = "I_Vocn";
    std::string SSH = "I_SSH";
    std::string MLD = "I_MLD"; // This one is optional
    std::string TAUX = "I_taux";
    std::string TAUY = "I_tauy";
    std::string TAUMOD = "I_taumod";
    std::string EMP = "I_fwflux";
    std::string QNOSUN = "I_rsnos";
    std::string QSW = "I_rsso";
    std::string SFLX = "I_sfi";
    std::string CICE = "I_conc";

    const std::vector<std::string> cplStringsIn
        = { SST, SSS, UOCEAN, VOCEAN, SSH  }; // MLD can be added to this one
    const std::vector<std::string> cplStringsOut
        = { TAUX, TAUY, TAUMOD, EMP, QNOSUN, QSW, SFLX, CICE };

};

} /* namespace Nextsim */

#endif /* OASISCOUPLEDOCEAN_HPP */
