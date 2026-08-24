/*!
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef FINITEELEMENTFLUXES_HPP
#define FINITEELEMENTFLUXES_HPP

#include "include/Configured.hpp"
#include "include/IFluxCalculation.hpp"
#include "include/IIceAlbedo.hpp"
#include "include/IIceOceanHeatFlux.hpp"
#include "include/ModelComponent.hpp"

namespace Nextsim {

//! A class that implements ice fluxes and open water fluxes according finiteelement.cpp.
class FiniteElementFluxes : public IFluxCalculation, public Configured<FiniteElementFluxes> {
public:
    FiniteElementFluxes()
        : iIceAlbedoImpl(nullptr)
        , Q_lh_owAccessor(getStore(), RO, ModelArray::Type::H)
        , Q_sh_owAccessor(getStore(), RO, ModelArray::Type::H)
        , Q_lw_owAccessor(getStore(), RO, ModelArray::Type::H)
        , Q_lh_iaAccessor(getStore(), RO, ModelArray::Type::H)
        , Q_sh_iaAccessor(getStore(), RO, ModelArray::Type::H)
        , Q_sw_iaAccessor(getStore(), RO, ModelArray::Type::H)
        , Q_lw_iaAccessor(getStore(), RO, ModelArray::Type::H)
        , rho_airAccessor(getStore(), RO, ModelArray::Type::H)
        , cp_airAccessor(getStore(), RO, ModelArray::Type::H)
        , sh_airAccessor(getStore(), RO, ModelArray::Type::H)
        , sh_waterAccessor(getStore(), RO, ModelArray::Type::H)
        , sh_iceAccessor(getStore(), RO, ModelArray::Type::H)
        , dshice_dTAccessor(getStore(), RO, ModelArray::Type::H)
        , sstAccessor(getStore())
        , sssAccessor(getStore())
        , t_airAccessor(getStore())
        , t_dew2Accessor(getStore())
        , p_airAccessor(getStore())
        , windSpeedAccessor(getStore())
        , u_airAccessor(getStore())
        , v_airAccessor(getStore())
        , hsnowAccessor(getStore())
        , ciceAccessor(getStore())
        , tsurfAccessor(getStore())
        , sw_inAccessor(getStore())
        , lw_inAccessor(getStore())
        , iceAlbedoAccessor(getStore())
        , icePenSWAccessor(getStore())
    {
    }
    ~FiniteElementFluxes() = default;

    enum {
        DRAGOCEANQ_KEY,
        DRAGOCEANT_KEY,
        DRAGICET_KEY,
        OCEANALBEDO_KEY,
        I0_KEY,
    };
    void configure() override;
    ConfigMap getConfiguration() const override;

    void setData(const ModelState::DataMap&) override;

    ModelState getStateDiagnostic() const override;

    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    std::string getName() const override { return "FiniteElementFluxes"; }

    void update(const TimestepTime& tst) override;

    //! Updates the fluxes over open water.
    void updateOW(const TimestepTime& tst);

    //! Updates the fluxes over ice.
    void updateIce(const TimestepTime& tst);

    //! Updates the atmospheric fluxes.
    void updateAtmosphere(const TimestepTime& tst);

private:
    // local namespace to prevent conflicts with other implementations
    struct Private {
        static inline constexpr TextTag Q_LH_OW = "Q_LH_OW";
        static inline constexpr TextTag Q_SH_OW = "Q_SH_OW";
        static inline constexpr TextTag Q_LW_OW = "Q_LW_OW";
        static inline constexpr TextTag Q_LH_IA = "Q_LH_IA";
        static inline constexpr TextTag Q_SH_IA = "Q_SH_IA";
        static inline constexpr TextTag Q_SW_IA = "Q_SW_IA";
        static inline constexpr TextTag Q_LW_IA = "Q_LW_IA";
        static inline constexpr TextTag RHO_AIR = "RHO_AIR";
        static inline constexpr TextTag CP_AIR = "CP_AIR";
        static inline constexpr TextTag SH_AIR = "SH_AIR";
        static inline constexpr TextTag SH_WATER = "SH_WATER";
        static inline constexpr TextTag SH_ICE = "SH_ICE";
        static inline constexpr TextTag DSHICE_DT = "DSHICE_DT";
    };
    // Owned diagnostic fields
    ModelArrayAccessor<Private::Q_LH_OW, RW> Q_lh_owAccessor; // Open water latent heat flux [W m⁻²]
    ModelArrayAccessor<Private::Q_SH_OW, RW>
        Q_sh_owAccessor; // Open water sensible heat flux [W m⁻²]
    ModelArrayAccessor<Private::Q_LW_OW, RW>
        Q_lw_owAccessor; // Open water net longwave radiative flux [W m⁻²]
    ModelArrayAccessor<Private::Q_LH_IA, RW> Q_lh_iaAccessor; // Ice latent heat flux [W m⁻²]
    ModelArrayAccessor<Private::Q_SH_IA, RW> Q_sh_iaAccessor; // Ice sensible heat flux [W m⁻²]
    ModelArrayAccessor<Private::Q_SW_IA, RW>
        Q_sw_iaAccessor; // Ice incident shortwave radiative flux [W m⁻²]
    ModelArrayAccessor<Private::Q_LW_IA, RW>
        Q_lw_iaAccessor; // Ice net longwave radiative flux [W m⁻²]
    // Derived air properties
    ModelArrayAccessor<Private::RHO_AIR, RW> rho_airAccessor;
    ModelArrayAccessor<Private::CP_AIR, RW>
        cp_airAccessor; // Specific heat capacity of the wet air [J kg⁻¹ K⁻¹]
    // Specific humidity and T derivative
    ModelArrayAccessor<Private::SH_AIR, RW> sh_airAccessor;
    ModelArrayAccessor<Private::SH_WATER, RW> sh_waterAccessor;
    ModelArrayAccessor<Private::SH_ICE, RW> sh_iceAccessor;
    ModelArrayAccessor<Private::DSHICE_DT, RW> dshice_dTAccessor;

    // Input fields
    ModelArrayAccessor<Protected::SST> sstAccessor;
    ModelArrayAccessor<Protected::SSS> sssAccessor;
    ModelArrayAccessor<Protected::T_AIR> t_airAccessor;
    ModelArrayAccessor<Protected::DEW_2M> t_dew2Accessor;
    ModelArrayAccessor<Protected::P_AIR> p_airAccessor;
    ModelArrayAccessor<Protected::WIND_SPEED> windSpeedAccessor;
    ModelArrayAccessor<Protected::WIND_U> u_airAccessor;
    ModelArrayAccessor<Protected::WIND_V> v_airAccessor;
    ModelArrayAccessor<Shared::H_SNOW_DG> hsnowAccessor; // cell-averaged value
    ModelArrayAccessor<Shared::C_ICE_DG> ciceAccessor;
    ModelArrayAccessor<Protected::T_SURF> tsurfAccessor;
    ModelArrayAccessor<Protected::SW_IN> sw_inAccessor;
    ModelArrayAccessor<Protected::LW_IN> lw_inAccessor;

    ModelArrayAccessor<Protected::ICE_ALBEDO> iceAlbedoAccessor;
    ModelArrayAccessor<Protected::ICE_PEN_SW> icePenSWAccessor;

    static FloatType dragOcean_q;
    static FloatType dragOcean_t;
    static FloatType dragIce_t;

    static FloatType m_oceanAlbedo;

    IIceAlbedo* iIceAlbedoImpl;
};

} /* namespace Nextsim */

#endif /* FINITEELEMENTFLUXES_HPP */
