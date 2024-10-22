/*!
 * @file IOceanBoundary.hpp
 *
 * @date Sep 12, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IOCEANBOUNDARY_HPP
#define IOCEANBOUNDARY_HPP

#include "include/ModelComponent.hpp"

namespace Nextsim {

namespace CouplingFields {
    constexpr TextTag SST = "SST"; // sea surface temperature ˚C
    constexpr TextTag SSS = "SSS"; // sea surface salinity PSU
    constexpr TextTag MLD = "MLD"; // Mixed layer or slab ocean depth m
    constexpr TextTag OCEAN_U = "U"; // x(east)-ward ocean current m s⁻¹
    constexpr TextTag OCEAN_V = "V"; // y(north)-ward ocean current m s⁻¹
}
//! An interface class for the oceanic inputs into the ice physics.
class IOceanBoundary : public ModelComponent {
public:
    IOceanBoundary()
    {
        m_couplingArrays.registerArray(CouplingFields::SST, &sst, RW);
        m_couplingArrays.registerArray(CouplingFields::SSS, &sss, RW);
        m_couplingArrays.registerArray(CouplingFields::OCEAN_U, &u, RW);
        m_couplingArrays.registerArray(CouplingFields::OCEAN_V, &v, RW);

        getStore().registerArray(Shared::Q_IO, &qio, RW);
        getStore().registerArray(Protected::SST, &sst, RO);
        getStore().registerArray(Protected::SSS, &sss, RO);
        getStore().registerArray(Protected::MLD, &mld, RO);
        getStore().registerArray(Protected::ML_BULK_CP, &cpml, RO);
        getStore().registerArray(Protected::TF, &tf, RO);
        getStore().registerArray(Protected::OCEAN_U, &u, RO);
        getStore().registerArray(Protected::OCEAN_V, &v, RO);
    }
    virtual ~IOceanBoundary() = default;

    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel&) const override { return getState(); }

    std::string getName() const override { return "IOceanBoundary"; }
    void setData(const ModelState::DataMap& ms) override
    {
        qio.resize();
        sst.resize();
        sss.resize();
        mld.resize();
        cpml.resize();
        tf.resize();
        u.resize();
        v.resize();

        if (ms.count("sst")) {
            sst = ms.at("sst");
        }
        if (ms.count("sss")) {
            sss = ms.at("sss");
        }
    }

    /*!
     * Performs the implementation specific updates before the physics calculations.
     *
     * @param tst The timestep start and duration .
     */
    virtual void updateBefore(const TimestepTime& tst) = 0;
    /*!
     *  Performs the implementation specific updates after the physics calculations.
     *
     * @param tst the timestep start and duration. Note that in some sense this
     *            update occurs near the end of the timestep at time tst.start + tst.duration
     */
    virtual void updateAfter(const TimestepTime& tst) = 0;

protected:
    HField qio; // Ice-ocean heat flux, W m⁻²
    HField sst; // Coupled or slab ocean sea surface temperature, ˚C
    HField sss; // Coupled or slab ocean sea surface salinity, PSU
    HField mld; // Mixed layer or slab ocean depth m
    HField tf; // Freezing point of the mixed layer, ˚C
    HField cpml; // Heat capacity of the mixed layer, J K⁻¹ m²
    UField u; // x(east)-ward ocean current, m s⁻¹
    VField v; // y(north)-ward ocean current, m s⁻¹

    ModelArrayReferenceStore m_couplingArrays;
};
} /* namespace Nextsim */

#endif /* IOCEANBOUNDARY_HPP */
