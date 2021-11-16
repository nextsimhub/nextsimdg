/*!
 * @file Configurator_test.cpp
 *
 * @date Oct 8, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "ArgV.hpp"
#include "Configurator.hpp"
#include "Configured.hpp"

#include <iostream>
#include <memory>
#include <sstream>
#include <string>

#include <boost/program_options.hpp>

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

class Config1 {
public:
    Config1()
        : value(0)
    {
    }
    void configure()
    {
        const std::string valueKey = "config.value";

        boost::program_options::options_description opt("Options");
        opt.add_options()(valueKey.c_str(), boost::program_options::value<int>()->default_value(-1),
            "Specify a value");
        boost::program_options::variables_map vm = Nextsim::Configurator::parse(opt);

        value = vm[valueKey].as<int>();
    }

    bool checkValue(int target) { return value == target; }

private:
    int value;
};

class Config2 : public Nextsim::Configured<Config2> {
public:
    enum {
        VALUE_KEY,
        NAME_KEY,
    };

    Config2();
    int getValue() { return value; }
    std::string getName() { return name; }

    void configure() override;

private:
    int value;
    std::string name;
};

template<>
const std::map<int, std::string> Nextsim::Configured<Config2>::keyMap = {
        {Config2::VALUE_KEY, "config.value"},
        {Config2::NAME_KEY, "config.name"},
};

Config2::Config2()
    : value(0)
{
addOption<int>(keyMap.at(VALUE_KEY), -1);
addOption<std::string>(keyMap.at(NAME_KEY), "");
}

void Config2::configure()
{
    value = retrieveValue<int>(keyMap.at(VALUE_KEY));
    name = retrieveValue<std::string>(keyMap.at(NAME_KEY));
}


class Config3 : public Nextsim::Configured<Config3> {
public:
    enum {
        VALUE_KEY,
        WEIGHT_KEY,
    };
    Config3()
    : value(0)
    , weight(0.)
{
}
    int getValue() { return value; }
    double getWeight() { return weight; }

    void configure() override;

private:
    int value;
    double weight;
};

template<>
const std::map<int, std::string> Nextsim::Configured<Config3>::keyMap = {
        {Config3::VALUE_KEY, "config.value"},
        {Config3::WEIGHT_KEY, "data.weight"},
};

void Config3::configure()
{
    value = Configured::getConfiguration(keyMap.at(VALUE_KEY), -1);
    weight = Configured::getConfiguration(keyMap.at(WEIGHT_KEY), 1.);
}

namespace Nextsim {

TEST_CASE("Parse one config stream using the raw configurator", "[Configurator]")
{
    // Since tests are not different execution environments, clear the streams from other tests.
    Configurator::clearStreams();
    Config1 config;

    int target = 42;
    std::stringstream text;
    text << "[config]" << std::endl << "value = " << target << std::endl;

    // Check for the initialized value
    REQUIRE(config.checkValue(0));

    config.configure();
    // Check for the default initialized value
    REQUIRE(config.checkValue(-1));

    std::unique_ptr<std::istream> pcstream(new std::stringstream(text.str()));

    Configurator::addStream(std::move(pcstream));

    config.configure();

    // Check for the configured value
    REQUIRE(config.checkValue(target));
}

TEST_CASE("Parse one config stream using the pointer configuration function", "[Configurator]")
{
    // Since tests are not different execution environments, clear the streams from other tests.
    Configurator::clearStreams();
    Config2::clearConfigurationMap();
    Config2 config;

    int target = 69105;
    std::string targetName = "Zork";
    std::stringstream text;
    text << "[config]" << std::endl
         << "value = " << target << std::endl
         << "name = " << targetName << std::endl;

    Configurator::addStream(std::unique_ptr<std::istream>(new std::stringstream(text.str())));

    tryConfigure(&config);

    REQUIRE(config.getValue() == target);
    REQUIRE(config.getName() == targetName);
}

TEST_CASE("Parse two config streams for one class, reference helper function", "[Configurator]")
{
    // Since tests are not different execution environments, clear the streams from other tests.
    Configurator::clearStreams();
    Config2::clearConfigurationMap();
    Config2 config;

    int target = 69105;
    std::string targetName = "Zork";
    std::stringstream text1;
    text1 << "[config]" << std::endl << "value = " << target << std::endl;
    Configurator::addStream(std::unique_ptr<std::istream>(new std::stringstream(text1.str())));

    std::stringstream text2;
    text2 << "[config]" << std::endl << "name = " << targetName << std::endl;
    Configurator::addStream(std::unique_ptr<std::istream>(new std::stringstream(text2.str())));

    Config2& ref = config;
    tryConfigure(ref);

    REQUIRE(config.getValue() == target);
    REQUIRE(config.getName() == targetName);
}

TEST_CASE("Parse config streams for two overlapping class, try", "[Configurator]")
{
    // Since tests are not different execution environments, clear the streams from other tests.
    Configurator::clearStreams();
    Config2::clearConfigurationMap();
    Config3::clearConfigurationMap();
    Config2 config;
    Config3 confih;


    int target = 69105;
    std::string targetName = "Zork II";
    std::stringstream text1;
    text1 << "[config]" << std::endl
          << "value = " << target << std::endl
          << "name = " << targetName << std::endl;
    Configurator::addStream(std::unique_ptr<std::istream>(new std::stringstream(text1.str())));

    double targetWeight = 0.467836;
    std::stringstream text2;
    text2 << "[data]" << std::endl << "weight = " << targetWeight << std::endl;
    Configurator::addStream(std::unique_ptr<std::istream>(new std::stringstream(text2.str())));

    // Test that the target values are not already present before configuring
    REQUIRE(config.getValue() != target);
    REQUIRE(config.getName() != targetName);

    REQUIRE(confih.getValue() != target);
    REQUIRE(confih.getWeight() != targetWeight);

    tryConfigure(config);
    tryConfigure(confih);

    REQUIRE(config.getValue() == target);
    REQUIRE(config.getName() == targetName);

    REQUIRE(confih.getValue() == target);
    REQUIRE(confih.getWeight() == targetWeight);
}

} /* namespace Nextsim */
