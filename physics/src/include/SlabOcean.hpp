/*!
 * @file SlabOcean.hpp
 *
 * @date 03 Jun 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SLABOCEAN_HPP
#define SLABOCEAN_HPP

#include "include/Configured.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"

namespace Nextsim {

/*!
 * A class to store and update the slab ocean that mediates between imposed
 * ocean boundary conditions and the sea ice. For more details, see §2.3 of
 * P. Rampal & al., "neXtSIM: a new Lagrangian sea ice model", The Cryosphere,
 * 10, 1055—1073 (2016)
 */
class SlabOcean : public ModelComponent, public Configured<SlabOcean> {
public:
    SlabOcean(ModelArrayReferenceStore& coupingArrays)
        : qdw(ModelArray::Type::H)
        , fdw(ModelArray::Type::H)
        , sstSlab(ModelArray::Type::H)
        , sssSlab(ModelArray::Type::H)
        , sstExt(getStore())
        , sssExt(getStore())
        , sst(getStore())
        , sss(getStore())
        , cpml(getStore())
        , qswNet(coupingArrays)
        , qNoSun(coupingArrays)
        , fwFlux(coupingArrays)
        , sFlux(coupingArrays)
    {
    }

    enum {
        TIMET_KEY,
        TIMES_KEY,
    };

    void configure() override;
    ConfigMap getConfiguration() const override;

    ModelState getStatePrognostic() const override;
    ModelState getStateDiagnostic() const override;

    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    void setData(const ModelState::DataMap& ms) override;
    std::string getName() const override { return "SlabOcean"; }

    void update(const TimestepTime&);

    static const double defaultRelaxationTime; // A default value for the relaxation time in s.

private:
    // Owned shared fields
    HField qdw;
    HField fdw;
    HField sstSlab;
    HField sssSlab;

    // Input fields
    ModelArrayRef<Protected::EXT_SST> sstExt;
    ModelArrayRef<Protected::EXT_SSS> sssExt;
    ModelArrayRef<Protected::SST> sst;
    ModelArrayRef<Protected::SSS> sss;
    ModelArrayRef<Protected::ML_BULK_CP> cpml;
    ModelArrayRef<CouplingFields::Q_SS_SW, RO> qswNet;
    ModelArrayRef<CouplingFields::Q_SS_NO_SW, RO> qNoSun;
    ModelArrayRef<CouplingFields::FWFLUX, RO> fwFlux;
    ModelArrayRef<CouplingFields::SFLUX, RO> sFlux;
    // TODO ModelArrayRef to assimilation flux

    double relaxationTimeT = defaultRelaxationTime;
    double relaxationTimeS = defaultRelaxationTime;

    double dt;

    void updateElement(size_t i, const TimestepTime& tst);
};

} /* namespace Nextsim */

#endif /* SLABOCEAN_HPP */
