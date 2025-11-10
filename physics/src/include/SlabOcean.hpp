/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SLABOCEAN_HPP
#define SLABOCEAN_HPP

#include "include/Configured.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArrayAccessor.hpp"
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
    SlabOcean(ModelArrayStore& couplingArrays)
        : qdwAccessor(getStore(), RO, ModelArray::Type::H)
        , fdwAccessor(getStore(), RO, ModelArray::Type::H)
        , sstSlabAccessor(getStore(), RO, ModelArray::Type::H)
        , sssSlabAccessor(getStore(), RO, ModelArray::Type::H)
        , sstExtAccessor(getStore())
        , sssExtAccessor(getStore())
        , sstAccessor(getStore())
        , sssAccessor(getStore())
        , cpmlAccessor(getStore())
        , qswNetAccessor(couplingArrays)
        , qNoSunAccessor(couplingArrays)
        , fwFluxAccessor(couplingArrays)
        , sFluxAccessor(couplingArrays)
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
    ModelArrayAccessor<Protected::SLAB_QDW, RW> qdwAccessor;
    ModelArrayAccessor<Protected::SLAB_FDW, RW> fdwAccessor;
    ModelArrayAccessor<Protected::SLAB_SST, RW> sstSlabAccessor;
    ModelArrayAccessor<Protected::SLAB_SSS, RW> sssSlabAccessor;

    // Input fields
    ModelArrayAccessor<Protected::EXT_SST> sstExtAccessor;
    ModelArrayAccessor<Protected::EXT_SSS> sssExtAccessor;
    ModelArrayAccessor<Protected::SST> sstAccessor;
    ModelArrayAccessor<Protected::SSS> sssAccessor;
    ModelArrayAccessor<Protected::ML_BULK_CP> cpmlAccessor;
    ModelArrayAccessor<CouplingFields::Q_SS_SW, RO> qswNetAccessor;
    ModelArrayAccessor<CouplingFields::Q_SS_NO_SW, RO> qNoSunAccessor;
    ModelArrayAccessor<CouplingFields::FWFLUX, RO> fwFluxAccessor;
    ModelArrayAccessor<CouplingFields::SFLUX, RO> sFluxAccessor;
    // TODO ModelArrayRef to assimilation flux

    double relaxationTimeT = defaultRelaxationTime;
    double relaxationTimeS = defaultRelaxationTime;

    double dt;

    void updateElement(size_t i, const TimestepTime& tst);
};

} /* namespace Nextsim */

#endif /* SLABOCEAN_HPP */
