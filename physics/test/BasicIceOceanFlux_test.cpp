/*!
 * @file BasicIceOceanFlux_test.cpp
 *
 * @date Sep 29, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

// Does this class need testing? Not really, but it got removed from
// FiniteElementFluxes_test and I thought the tests should continue to exist
// somewhere

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "include/BasicIceOceanHeatFlux.hpp"

#include "include/IOceanBoundary.hpp"
#include "include/ModelArray.hpp"
#include "include/Time.hpp"
#include "include/UnescoFreezing.hpp"
#include "include/constants.hpp"

namespace Nextsim {

TEST_CASE("Melting conditions", "[BasicIceOceanHeatFlux]")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 1 });

    class OceanData : public IOceanBoundary {
    public:
        OceanData()
            : IOceanBoundary()
        {
        }
        void setData(const ModelState::DataMap& state) override
        {
            IOceanBoundary::setData(state);
            UnescoFreezing uf;
            sst = -1.;
            sss = 32.;
            mld = 10.25;
            tf = uf(sss[0]);
            cpml = Water::cp * Water::rho * mld[0];
            u = 0;
            v = 0;
        }
        void updateBefore(const TimestepTime& tst) override { }
        void updateAfter(const TimestepTime& tst) override { }
    } ocnBdy;
    ocnBdy.setData(ModelState().data);

    class ProgData : public ModelComponent {
    public:
        ProgData()
        {
            registerProtectedArray(ProtectedArray::H_ICE, &hice);
            registerProtectedArray(ProtectedArray::C_ICE, &cice);
            registerProtectedArray(ProtectedArray::H_SNOW, &hsnow);
            registerProtectedArray(ProtectedArray::T_ICE, &tice0);
            registerProtectedArray(ProtectedArray::HTRUE_ICE, &hice0);
            registerProtectedArray(ProtectedArray::HTRUE_SNOW, &hsnow0);
        }
        std::string getName() const override { return "ProgData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            cice[0] = 0.5;
            hice[0] = 0.1; // Here we are using the cell-averaged thicknesses
            hsnow[0] = 0.01;
            tice0[0] = -1.;

            hice0[0] = hice[0] / cice[0];
            hsnow0[0] = hsnow[0] / cice[0];
        }

        HField hice;
        HField cice;
        HField hsnow;
        HField tice0;
        HField hice0; // ice averaged ice thickness
        HField hsnow0; // ice averaged snow thickness
        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } iceState;
    iceState.setData(ModelState().data);

    HField qio;
    qio.resize();
    ModelComponent::registerExternalSharedArray(ModelComponent::SharedArray::Q_IO, &qio);

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };
    BasicIceOceanHeatFlux biohf;
    biohf.update(tst);

    double prec = 1e-5;
    REQUIRE(qio[0] == Approx(53717.8).epsilon(prec));
}

TEST_CASE("Freezing conditions", "[BasicIceOceanHeatFlux]")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 1 });

    class OceanData : public IOceanBoundary {
    public:
        OceanData()
            : IOceanBoundary()
        {
        }
        void setData(const ModelState::DataMap& state) override
        {
            IOceanBoundary::setData(state);
            UnescoFreezing uf;
            sst = -1.75;
            sss = 32.;
            mld = 10.25;
            tf = uf(sss[0]);
            cpml = Water::cp * Water::rho * mld[0];
            u = 0;
            v = 0;
        }
        void updateBefore(const TimestepTime& tst) override { }
        void updateAfter(const TimestepTime& tst) override { }
    } ocnBdy;
    ocnBdy.setData(ModelState().data);

    class ProgData : public ModelComponent {
    public:
        ProgData()
        {
            registerProtectedArray(ProtectedArray::H_ICE, &hice);
            registerProtectedArray(ProtectedArray::C_ICE, &cice);
            registerProtectedArray(ProtectedArray::H_SNOW, &hsnow);
            registerProtectedArray(ProtectedArray::T_ICE, &tice0);
            registerProtectedArray(ProtectedArray::HTRUE_ICE, &hice0);
            registerProtectedArray(ProtectedArray::HTRUE_SNOW, &hsnow0);
        }
        std::string getName() const override { return "ProgData"; }

        void setData(const ModelState::DataMap&) override
        {
            noLandMask();
            cice[0] = 0.5;
            hice[0] = 0.1; // Here we are using the cell-averaged thicknesses
            hsnow[0] = 0.01;
            tice0[0] = -9.;

            hice0[0] = hice[0] / cice[0];
            hsnow0[0] = hsnow[0] / cice[0];
        }

        HField hice;
        HField cice;
        HField hsnow;
        HField tice0;
        HField hice0; // ice averaged ice thickness
        HField hsnow0; // ice averaged snow thickness

        ModelState getState() const override { return ModelState(); }
        ModelState getState(const OutputLevel&) const override { return getState(); }
    } iceState;
    iceState.setData(ModelState().data);

    HField qio;
    qio.resize();
    ModelComponent::registerExternalSharedArray(ModelComponent::SharedArray::Q_IO, &qio);

    TimestepTime tst = { TimePoint("2000-001"), Duration("P0-0T0:10:0") };
    BasicIceOceanHeatFlux biohf;
    biohf.update(tst);

    double prec = 1e-5;
    REQUIRE(qio[0] == Approx(73.9465).epsilon(prec));
}
} // namespace Nextsim
