/*!
 * @file ConfigOutput_test.cpp
 *
 * @date 11 May 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/ConfigOutput.hpp"

#include "include/IStructure.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/ModelMetadata.hpp"
#include "include/ModelState.hpp"
#include "include/Module.hpp"
#include "include/NZLevels.hpp"

#include <sstream>

namespace Nextsim {

const std::string hiceName = "hice";
const std::string ciceName = "cice";
const std::string hsnowName = "hsnow";
const std::string ticeName = "tice";

TEST_SUITE_BEGIN("ConfigOutput");
TEST_CASE("Test periodic output")
{
    size_t nx = 2;
    size_t ny = 5;
    size_t nz = 3;
    NZLevels::set(nz);

    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
    ModelArray::setDimension(ModelArray::Dimension::Z, NZLevels::get());

    std::stringstream config;
    config << "[Modules]" << std::endl;
    config << "Nextsim::IDiagnosticOutput = Nextsim::ConfigOutput" << std::endl;
    config << std::endl;
    config << "[ConfigOutput]" << std::endl;
    config << "period = 3600" << std::endl; // Output every hour
    config << "start = 2020-01-11T00:00:00Z" << std::endl; // start after 10 days
    config << "field_names = " << hiceName << "," << ciceName << "," << ticeName << std::endl;
    config << "filename = diag%m%d.nc" << std::endl;
    config << "file_period = 86400" << std::endl; // Files every day

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    ConfiguredModule::parseConfigurator();

    Module::setImplementation<IStructure>("ParametricGrid");

    HField hice(ModelArray::Type::H);
    HField cice(ModelArray::Type::H);
    HField hsnow(ModelArray::Type::H);
    ZField tice(ModelArray::Type::Z);

    hice.resize();
    cice.resize();
    hsnow.resize();
    tice.resize();

    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::H_ICE, &hice);
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::C_ICE, &cice);
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::H_SNOW, &hsnow);
    ModelComponent::registerExternalProtectedArray(ModelComponent::ProtectedArray::T_ICE, &tice);

    ModelMetadata meta;
    meta.setTime(TimePoint("2020-01-01T00:00:00Z"));

    IDiagnosticOutput& ido = Module::getImplementation<IDiagnosticOutput>();
    tryConfigure(ido);

    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                tice(i, j, k) = 0.1 * k + 0.4 + 0.01 * (j * nx + i);
            }
        }
    }
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            hice(i, j) = 0 + 0.01 * (j * nx + i);
            cice(i, j) = 0.1 + 0.01 * (j * nx + i);
            hsnow(i, j) = 0.2 + 0.01 * (j * nx + i);
        }
    }
    for (size_t day = 1; day <= 20; ++day) {
        double dayIncr = 100.;
        hice += dayIncr;
        cice += dayIncr;
        hsnow += dayIncr;
        tice += dayIncr;
        for (size_t hour = 0; hour < 24; ++hour) {
            double hourIncr = 1;
            hice += hourIncr;
            cice += hourIncr;
            hsnow += hourIncr;
            ModelState state;

            ido.outputState(meta);
            meta.incrementTime(Duration(3600.));
        }
    }
}
TEST_SUITE_END();
}
