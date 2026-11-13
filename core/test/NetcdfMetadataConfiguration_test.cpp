/*!
 * @author  M.A. Kowalski <mak60@cam.ac.uk>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/NetcdfMetadataConfiguration.hpp"

#include <ncDim.h>
#include <ncFile.h>
#include <ncVar.h>

#include <cstdio>
#include <sstream>
#include <string>

namespace Nextsim {

struct RAIITempFile {
    std::string filename;
    RAIITempFile()
        : filename(std::tmpnam(nullptr))
    {
    }
    ~RAIITempFile() { std::remove(filename.c_str()); }
};

TEST_SUITE_BEGIN("NetcdfMetadataConfiguration");

TEST_CASE("Write and read a netCDF configuration file")
{
    const RAIITempFile tempFile;

    // Test values for each type
    const double doubleValue = 3.14159;
    const float floatValue = 2.71828f;
    const int intValue = 42;
    const unsigned uintValue = 65535u;
    const std::string stringValue = "test string";

    // Write the configuration variables to a temporary netCDF file
    {

        netCDF::NcFile ncFile(tempFile.filename, netCDF::NcFile::replace);

        // Data that will be read
        netCDF::NcVar doubleVar = ncFile.addVar("double", netCDF::ncDouble);
        doubleVar.putVar(&doubleValue);

        netCDF::NcVar floatVar = ncFile.addVar("float", netCDF::ncFloat);
        floatVar.putVar(&floatValue);

        netCDF::NcVar intVar = ncFile.addVar("int", netCDF::ncInt);
        intVar.putVar(&intValue);

        netCDF::NcVar uintVar = ncFile.addVar("uint", netCDF::ncUint);
        uintVar.putVar(&uintValue);

        netCDF::NcVar stringVar = ncFile.addVar("string", netCDF::ncString);
        const char* cString = stringValue.c_str();
        stringVar.putVar(&cString);

        // Data that will not be read as configuration option
        {
            // Rank 1 variable
            // We make is have decent size to ensure trigger segfaults if read
            // accidentally.
            constexpr int size = 100;
            std::array<double, size> data;
            data.fill(17.0);
            netCDF::NcDim xDim = ncFile.addDim("x_array", size);
            netCDF::NcVar arrayVar = ncFile.addVar("array", netCDF::ncDouble, xDim);
            arrayVar.putVar(data.data());
        }
        {
            // Rank 2 variable
            constexpr int size_0 = 10;
            constexpr int size_1 = 20;
            std::array<double, size_0 * size_1> data;
            data.fill(42.0);
            netCDF::NcDim xDim = ncFile.addDim("x_array2d", size_0);
            netCDF::NcDim yDim = ncFile.addDim("y_array2d", size_1);
            netCDF::NcVar arrayVar = ncFile.addVar("array2d", netCDF::ncDouble, { xDim, yDim });
            arrayVar.putVar(data.data());
        }
    }

    // We explicitly don't want to use the static parameter in Model::restartOptionName
    // We want to test that what we provide is respected
    // Point the configurator to our temporary file
    const std::string localRestartOptionName = "config_file";
    {
        Configurator::clear();
        std::stringstream configStream;
        configStream << localRestartOptionName << " = " << tempFile.filename << std::endl;
        Configurator::addSStream(configStream);
    }

    // Get the configuration stored in the restart file
    Nextsim::NetcdfMetadataConfiguration ncdfMC;
    Nextsim::Configurator::setAdditionalConfiguration(&ncdfMC);
    Nextsim::Configurator::getAdditionalConfiguration(localRestartOptionName);

    // Read the configuration

    boost::program_options::options_description opt;
    opt.add_options()("double", boost::program_options::value<double>(), "")("float",
        boost::program_options::value<float>(), "")("int", boost::program_options::value<int>(),
        "")("uint", boost::program_options::value<unsigned>(), "")(
        "string", boost::program_options::value<std::string>(), "");
    auto vm = Configurator::parse(opt);
    REQUIRE(vm["double"].as<double>() == doctest::Approx(doubleValue));
    REQUIRE(vm["float"].as<float>() == doctest::Approx(floatValue));
    REQUIRE(vm["int"].as<int>() == intValue);
    REQUIRE(vm["uint"].as<unsigned>() == uintValue);
    REQUIRE(vm["string"].as<std::string>() == stringValue);
}

TEST_SUITE_END();

} /* namespace Nextsim */
