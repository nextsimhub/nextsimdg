/*!
 * @file SlabOcean.hpp
 *
 * @date 27 Jan 2025
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
        : fdw(ModelArray::Type::H)
        , qdw(ModelArray::Type::H)
        , cice(getStore())
        , cpml(getStore())
        , deltaHice(getStore())
        , deltaSmelt(getStore())
        , emp(getStore())
        , mld(getStore())
        , newIce(getStore())
        , qio(getStore())
        , qow(getStore())
        , sss(getStore())
        , sssExt(getStore())
        , sst(getStore())
        , sstExt(getStore())
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

    void update(const TimestepTime&);

    static const double defaultRelaxationTime; // A default value for the relaxation time in s.

private:
    // Owned shared fields
    HField fdw;
    HField qdw;

    // Input fields
    ModelArrayRef<Protected::EVAP_MINUS_PRECIP> emp;
    ModelArrayRef<Protected::EXT_SSS> sssExt;
    ModelArrayRef<Protected::EXT_SST> sstExt;
    ModelArrayRef<Protected::MLD> mld;
    ModelArrayRef<Protected::ML_BULK_CP> cpml;
    ModelArrayRef<Shared::C_ICE> cice;
    ModelArrayRef<Shared::DELTA_HICE, RW> deltaHice;
    ModelArrayRef<Shared::HSNOW_MELT, RW> deltaSmelt;
    ModelArrayRef<Shared::NEW_ICE, RW> newIce;
    ModelArrayRef<Shared::Q_IO, RW> qio;
    ModelArrayRef<Shared::Q_OW, RW> qow;
    ModelArrayRef<Shared::SSS, RW> sss;
    ModelArrayRef<Shared::SST, RW> sst;
    // TODO ModelArrayRef to assimilation flux

    static const std::string sstSlabName;
    static const std::string sssSlabName;

    double relaxationTimeT = defaultRelaxationTime;
    double relaxationTimeS = defaultRelaxationTime;

    void updateElements(size_t i, const TimestepTime& tst);
};

} /* namespace Nextsim */

#endif /* SLABOCEAN_HPP */
