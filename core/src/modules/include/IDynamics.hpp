/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IDYNAMICS_HPP
#define IDYNAMICS_HPP

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
        : uice(ModelArray::Type::H)
        , vice(ModelArray::Type::H)
        , taux(ModelArray::Type::H)
        , tauy(ModelArray::Type::H)
        , shear(ModelArray::Type::H)
        , divergence(ModelArray::Type::H)
        , sigmaI(ModelArray::Type::H)
        , sigmaII(ModelArray::Type::H)
        , stress11(ModelArray::Type::DGSTRESS)
        , stress12(ModelArray::Type::DGSTRESS)
        , stress22(ModelArray::Type::DGSTRESS)
        , damage(getStore())
        , uwind(getStore())
        , vwind(getStore())
        , uocean(getStore())
        , vocean(getStore())
        , ssh(getStore())
        , m_usesDamage(usesDamageIn)
        , hiceDG(getStore())
        , ciceDG(getStore())
        , hsnowDG(getStore())
    {
        getStore().registerArray(Protected::DIV, &divergence, RO);
        getStore().registerArray(Protected::ICE_U, &uice, RO);
        getStore().registerArray(Protected::ICE_V, &vice, RO);
        getStore().registerArray(Protected::IO_STRESS_X, &taux, RO);
        getStore().registerArray(Protected::IO_STRESS_Y, &tauy, RO);
        getStore().registerArray(Protected::SHEAR, &shear, RO);
        getStore().registerArray(Protected::SIGMAI, &sigmaI, RO);
        getStore().registerArray(Protected::SIGMAII, &sigmaII, RO);

#ifdef USE_XIOS
        // Set XIOS field types
        Xios& xiosHandler = Xios::getInstance();
        xiosHandler.setPrognosticFieldType(coordsName, ModelArray::Type::VERTEX);
        xiosHandler.setPrognosticFieldType(hiceName, ModelArray::AdvectionType);
        xiosHandler.setPrognosticFieldType(ciceName, ModelArray::AdvectionType);
        xiosHandler.setPrognosticFieldType(damageName, ModelArray::AdvectionType);
        xiosHandler.setPrognosticFieldType(hsnowName, ModelArray::AdvectionType);
#endif
    }
    virtual ~IDynamics() = default;

    ModelState getStatePrognostic() const override
    {
        ModelState state = { {
                                 { uName, uice },
                                 { vName, vice },
                             },
            getConfiguration() };

        if (m_usesDamage) {
            ModelState::DataMap damageState = { { damageName, damage } };
            state.merge(damageState);
        }

        return state;
    }

    ModelState getStateDiagnostic() const override
    {
        ModelState state = { {
                                 { uIOStressName, taux },
                                 { vIOStressName, tauy },
                                 { uName, uice },
                                 { vName, vice },
                                 { shearName, shear },
                                 { divergenceName, divergence },
                                 { sigmaIName, sigmaI },
                                 { sigmaIIName, sigmaII },
                                 { stress11Name, stress11 },
                                 { stress11Name, stress12 },
                                 { stress11Name, stress22 },
                             },
            {} };
        return state.merge(getStatePrognostic());
    }

    std::string getName() const override { return "IDynamics"; }
    virtual void setData(const ModelState::DataMap& ms) override
    {
        uice.resize();
        vice.resize();
        shear.resize();
        divergence.resize();
        sigmaI.resize();
        sigmaII.resize();
        stress11.resize();
        stress12.resize();
        stress22.resize();
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

    virtual void prepareAdvection() = 0;

protected:
    // Shared ice velocity arrays
    HField uice;
    HField vice;
    // Ice-ocean stress (for the coupler, mostly)
    HField taux;
    HField tauy;
    // Diagnostic outputs of shear, divergence and the stress invariants
    HField shear;
    HField divergence;
    HField sigmaI;
    HField sigmaII;
    // Stress components. Diagnostic for some dynamics, prognostic for brittle and others.
    DGSField stress11;
    DGSField stress12;
    DGSField stress22;
    // References to the DG0 finite volume data arrays
    ModelArrayRef<Shared::DAMAGE, RW> damage;

    // References to the forcing velocity arrays
    ModelArrayRef<Protected::WIND_U> uwind;
    ModelArrayRef<Protected::WIND_V> vwind;
    ModelArrayRef<Protected::OCEAN_U> uocean;
    ModelArrayRef<Protected::OCEAN_V> vocean;
    ModelArrayRef<Protected::SSH> ssh;

    // Does this implementation of the dynamics use damage?
    bool m_usesDamage;

    // Store the h_ice and c_ice DG fields here, rather than in the kernel.
    ModelArrayRef<Shared::H_ICE_DG, RW> hiceDG;
    ModelArrayRef<Shared::C_ICE_DG, RW> ciceDG;
    ModelArrayRef<Shared::H_SNOW_DG, RW> hsnowDG;

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
