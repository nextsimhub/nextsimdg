/*!
 * @file SlabOcean.hpp
 *
 * @date 27 Jan 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SLABOCEAN_HPP
#define SLABOCEAN_HPP

#include "include/Configured.hpp"
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
        , sstExt(getProtectedArray())
        , sssExt(getProtectedArray())
        , sst(getProtectedArray())
        , sss(getProtectedArray())
        , mld(getProtectedArray())
        , cpml(getProtectedArray())
        , emp(getProtectedArray())
        , cice(getProtectedArray())
        , qio(getSharedArray())
        , qow(getSharedArray())
        , newIce(getSharedArray())
        , deltaHice(getSharedArray())
        , deltaSmelt(getSharedArray())
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
    ModelArrayRef<ProtectedArray::EXT_SST, MARConstBackingStore> sstExt;
    ModelArrayRef<ProtectedArray::EXT_SSS, MARConstBackingStore> sssExt;
    ModelArrayRef<ProtectedArray::SST, MARConstBackingStore> sst;
    ModelArrayRef<ProtectedArray::SSS, MARConstBackingStore> sss;
    ModelArrayRef<ProtectedArray::MLD, MARConstBackingStore> mld;
    ModelArrayRef<ProtectedArray::ML_BULK_CP, MARConstBackingStore> cpml;
    ModelArrayRef<ProtectedArray::EVAP_MINUS_PRECIP, MARConstBackingStore> emp;
    ModelArrayRef<ProtectedArray::C_ICE, MARConstBackingStore> cice;
    ModelArrayRef<SharedArray::Q_IO, MARBackingStore> qio;
    ModelArrayRef<SharedArray::Q_OW, MARBackingStore> qow;
    // TODO ModelArrayRef to assimilation flux
    ModelArrayRef<SharedArray::NEW_ICE, MARBackingStore> newIce;
    ModelArrayRef<SharedArray::DELTA_HICE, MARBackingStore> deltaHice;
    ModelArrayRef<SharedArray::HSNOW_MELT, MARBackingStore> deltaSmelt;

    static const std::string sstSlabName;
    static const std::string sssSlabName;

    double relaxationTimeT = defaultRelaxationTime;
    double relaxationTimeS = defaultRelaxationTime;
};

} /* namespace Nextsim */

#endif /* SLABOCEAN_HPP */
