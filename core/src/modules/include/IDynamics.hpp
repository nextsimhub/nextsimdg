/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IDYNAMICS_HPP
#define IDYNAMICS_HPP

#include "include/ModelArrayAccessor.hpp"
#include "include/ModelComponent.hpp"
#include "include/Time.hpp"
#ifdef USE_XIOS
#include "include/Xios.hpp"
#endif
#include "include/gridNames.hpp"

#include <limits>

namespace Nextsim {
class IDynamics : public ModelComponent {
public:
    IDynamics()
        : IDynamics(false)
    {
    }
    IDynamics(bool usesDamageIn)
        : uiceAccessor(getStore(), RO, ModelArray::Type::H)
        , viceAccessor(getStore(), RO, ModelArray::Type::H)
        , tauxAccessor(getStore(), RO, ModelArray::Type::H)
        , tauyAccessor(getStore(), RO, ModelArray::Type::H)
        , shearAccessor(getStore(), RO, ModelArray::Type::H)
        , divergenceAccessor(getStore(), RO, ModelArray::Type::H)
        , sigmaIAccessor(getStore(), RO, ModelArray::Type::H)
        , sigmaIIAccessor(getStore(), RO, ModelArray::Type::H)
        , damageAccessor(getStore())
        , uwindAccessor(getStore())
        , vwindAccessor(getStore())
        , uoceanAccessor(getStore())
        , voceanAccessor(getStore())
        , sshAccessor(getStore())
        , m_usesDamage(usesDamageIn)
        , hiceDGAccessor(getStore())
        , ciceDGAccessor(getStore())
        , hsnowDGAccessor(getStore())
    {

#ifdef USE_XIOS
        // Set XIOS field types
        Xios& xiosHandler = Xios::getInstance();

        // Advected fields
        xiosHandler.setPrognosticFieldType(hiceName, ModelArray::AdvectionType);
        xiosHandler.setPrognosticFieldType(ciceName, ModelArray::AdvectionType);
        xiosHandler.setPrognosticFieldType(damageName, ModelArray::Type::H);
        xiosHandler.setPrognosticFieldType(hsnowName, ModelArray::AdvectionType);
#endif
    }
    virtual ~IDynamics() = default;

    ModelState getStatePrognostic() const override
    {
        ModelState state = { {
                                 { uName, uiceAccessor.getHostRO() },
                                 { vName, viceAccessor.getHostRO() },
                             },
            getConfiguration() };

        if (m_usesDamage) {
            ModelState::DataMap damageState = { { damageName, damageAccessor.getHostRO() } };
            state.merge(damageState);
        }

        return state;
    }

    ModelState getStateDiagnostic() const override
    {
        ModelState state = { {
                                 { uIOStressName, tauxAccessor.getHostRO() },
                                 { vIOStressName, tauyAccessor.getHostRO() },
                                 { uName, uiceAccessor.getHostRO() },
                                 { vName, viceAccessor.getHostRO() },
                                 { shearName, shearAccessor.getHostRO() },
                                 { divergenceName, divergenceAccessor.getHostRO() },
                                 { sigmaIName, sigmaIAccessor.getHostRO() },
                                 { sigmaIIName, sigmaIIAccessor.getHostRO() },
                             },
            {} };
        return state.merge(getStatePrognostic());
    }

    std::string getName() const override { return "IDynamics"; }
    virtual void setData(const ModelState::DataMap& ms) override
    {
        HField& uice = uiceAccessor.getHostRW();
        uice.reinitialize();
        HField& vice = viceAccessor.getHostRW();
        vice.reinitialize();
        HField& shear = shearAccessor.getHostRW();
        shear.reinitialize();
        HField& divergence = divergenceAccessor.getHostRW();
        divergence.reinitialize();
        HField& sigmaI = sigmaIAccessor.getHostRW();
        sigmaI.reinitialize();
        HField& sigmaII = sigmaIIAccessor.getHostRW();
        sigmaII.reinitialize();
    }

    virtual void update(const TimestepTime& tst) = 0;

    /*!
     * Returns whether the dynamics implementation uses the damage field.
     */
    virtual bool usesDamage() const { return m_usesDamage; }

    virtual void advectField(double timestep, ModelArray& field,
        double lowerLimit = -std::numeric_limits<double>::infinity(),
        double upperLimit = std::numeric_limits<double>::infinity())
    {
    }

#ifdef USE_KOKKOS
    virtual void advectField(double timestep, const DeviceViewMA& field,
        double lowerLimit = -std::numeric_limits<double>::infinity(),
        double upperLimit = std::numeric_limits<double>::infinity())
    {
    }
#endif

    virtual void prepareAdvection() = 0;

protected:
    // Shared ice velocity arrays
    ModelArrayAccessor<Protected::ICE_U, RW> uiceAccessor;
    ModelArrayAccessor<Protected::ICE_V, RW> viceAccessor;
    // Ice-ocean stress (for the coupler, mostly)
    ModelArrayAccessor<Protected::IO_STRESS_X, RW> tauxAccessor;
    ModelArrayAccessor<Protected::IO_STRESS_Y, RW> tauyAccessor;
    // Diagnostic outputs of shear, divergence and the stress invariants
    ModelArrayAccessor<Protected::SHEAR, RW> shearAccessor;
    ModelArrayAccessor<Protected::DIV, RW> divergenceAccessor;
    ModelArrayAccessor<Protected::SIGMAI, RW> sigmaIAccessor;
    ModelArrayAccessor<Protected::SIGMAII, RW> sigmaIIAccessor;
    // References to the DG0 finite volume data arrays
    ModelArrayAccessor<Shared::DAMAGE, RW> damageAccessor;

    // References to the forcing velocity arrays
    ModelArrayAccessor<Protected::WIND_U> uwindAccessor;
    ModelArrayAccessor<Protected::WIND_V> vwindAccessor;
    ModelArrayAccessor<Protected::OCEAN_U> uoceanAccessor;
    ModelArrayAccessor<Protected::OCEAN_V> voceanAccessor;
    ModelArrayAccessor<Protected::SSH> sshAccessor;

    // Does this implementation of the dynamics use damage?
    bool m_usesDamage;

    // Store the h_ice and c_ice DG fields here, rather than in the kernel.
    ModelArrayAccessor<Shared::H_ICE_DG, RW> hiceDGAccessor;
    ModelArrayAccessor<Shared::C_ICE_DG, RW> ciceDGAccessor;
    ModelArrayAccessor<Shared::H_SNOW_DG, RW> hsnowDGAccessor;

    /*
     * Checks and returns if the provided data map is spherical
     */
    static bool checkSpherical(const ModelState::DataMap& ms)
    {
        // Decide between Cartesian (x & y) and spherical (longitude & latitude)
        if (ms.count(longitudeName) > 0 && ms.count(latitudeName) > 0) {
            return true;
        } else if (ms.count(xName) > 0 && ms.count(yName) > 0) {
            return false;
        } else {
            // Throw a runtime_error exception which can either be handled or not
            throw std::runtime_error("Input data must contain either Cartesian (" + xName + ", "
                + yName + ") or spherical (" + longitudeName + ", " + latitudeName
                + ") coordinates.");
        }
    }
};
}

#endif /* IDYNAMICS_HPP */
