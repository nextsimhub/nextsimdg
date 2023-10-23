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

//! An interface class for the oceanic inputs into the ice physics.
class IOceanBoundary : public ModelComponent {
public:
    IOceanBoundary()
    {
        m_couplingArrays.resize(static_cast<size_t>(CouplingFields::COUNT));
        m_couplingArrays[static_cast<size_t>(CouplingFields::SST)] = &sst;
        m_couplingArrays[static_cast<size_t>(CouplingFields::SSS)] = &sss;
        m_couplingArrays[static_cast<size_t>(CouplingFields::OCEAN_U)] = &u;
        m_couplingArrays[static_cast<size_t>(CouplingFields::OCEAN_V)] = &v;

        registerSharedArray(SharedArray::Q_IO, &qio);
        registerProtectedArray(ProtectedArray::SST, &sst);
        registerProtectedArray(ProtectedArray::SSS, &sss);
        registerProtectedArray(ProtectedArray::MLD, &mld);
        registerProtectedArray(ProtectedArray::ML_BULK_CP, &cpml);
        registerProtectedArray(ProtectedArray::TF, &tf);
        registerProtectedArray(ProtectedArray::OCEAN_U, &u);
        registerProtectedArray(ProtectedArray::OCEAN_V, &v);
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
    enum class CouplingFields {
        SST, // sea surface temperature ˚C
        SSS, // sea surface salinity PSU
        MLD, // Mixed layer or slab ocean depth m
        OCEAN_U, // x(east)-ward ocean current m s⁻¹
        OCEAN_V, // y(north)-ward ocean current m s⁻¹
        COUNT
    };
    HField qio; // Ice-ocean heat flux, W m⁻²
    HField sst; // Coupled or slab ocean sea surface temperature, ˚C
    HField sss; // Coupled or slab ocean sea surface salinity, PSU
    HField mld; // Mixed layer or slab ocean depth m
    HField tf; // Freezing point of the mixed layer, ˚C
    HField cpml; // Heat capacity of the mixed layer, J K⁻¹ m²
    UField u; // x(east)-ward ocean current, m s⁻¹
    VField v; // y(north)-ward ocean current, m s⁻¹

    MARBackingStore m_couplingArrays;
};
} /* namespace Nextsim */

#endif /* IOCEANBOUNDARY_HPP */
