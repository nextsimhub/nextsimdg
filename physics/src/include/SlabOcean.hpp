/*!
 * @file SlabOcean.hpp
 *
 * @date 7 Sep 2023
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
    SlabOcean()
        : qdw(ModelArray::Type::H)
        , fdw(ModelArray::Type::H)
        , sstSlab(ModelArray::Type::H)
        , sssSlab(ModelArray::Type::H)
        , sstExt(getStore())
        , sssExt(getStore())
        , sst(getStore())
        , sss(getStore())
        , mld(getStore())
        , cpml(getStore())
        , emp(getStore())
        , cice(getStore())
        , qio(getStore())
        , qow(getStore())
        , newIce(getStore())
        , deltaHice(getStore())
        , deltaSmelt(getStore())
    {
    }

    enum {
        TIMET_KEY,
        TIMES_KEY,
    };

    void configure() override;
    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    void setData(const ModelState::DataMap& ms) override;
    ModelState getState() const override;
    ModelState getState(const OutputLevel&) const override;
    std::string getName() const override { return "SlabOcean"; }

    std::unordered_set<std::string> hFields() const override;
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
    ModelArrayRef<Protected::MLD> mld;
    ModelArrayRef<Protected::ML_BULK_CP> cpml;
    ModelArrayRef<Protected::EVAP_MINUS_PRECIP> emp;
    ModelArrayRef<Protected::C_ICE> cice;
    ModelArrayRef<Shared::Q_IO, RW> qio;
    ModelArrayRef<Shared::Q_OW, RW> qow;
    // TODO ModelArrayRef to assimilation flux
    ModelArrayRef<Shared::NEW_ICE, RW> newIce;
    ModelArrayRef<Shared::DELTA_HICE, RW> deltaHice;
    ModelArrayRef<Shared::HSNOW_MELT, RW> deltaSmelt;

    static const std::string sstSlabName;
    static const std::string sssSlabName;

    double relaxationTimeT = defaultRelaxationTime;
    double relaxationTimeS = defaultRelaxationTime;
};

} /* namespace Nextsim */

#endif /* SLABOCEAN_HPP */
