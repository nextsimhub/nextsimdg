/*!
 * @file Configurator_test.cpp
 *
 * @date Oct 8, 2021
 * @author Tim Spain
 */

#include "Configurator.hpp"
#include "Configured.hpp"
#include "ArgV.hpp"

#include <string>
#include <iostream>
#include <sstream>
#include <memory>

#include <boost/program_options.hpp>

#define CATCH_CONFIG_MAIN
//#include <catch2/catch.hpp>
#include "/opt/home/include/catch2/catch.hpp"

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
        opt.add_options()
                (valueKey.c_str(), boost::program_options::value<int>()->default_value(1), "Specify a value")
                ;
        boost::program_options::variables_map vm = Nextsim::Configurator::parse(opt);
        std::cerr << "Read configuration" << std::endl;
        std::cerr << vm[valueKey].as<int>();
        std::cerr << "Assigned value" << std::endl;
    }

    bool checkValue(int target)
    {
        return value == target;
    }
private:
    int value;
};

class Config2 : public Nextsim::Configured {
public:
    Config2()
        : value(0)
    {}
private:
    int value;
    std::string name;
};

class Config3 : public Nextsim::Configured {
public:
    Config3()
        : value(0)
        , weight(0.)
    {}
private:
    int value;
    double weight;
};

namespace Nextsim {

TEST_CASE("Parse one config stream using the raw configurator", "[Configurator]")
{
    Config1 config;

    int target = 42;
    std::stringstream text;
    text << "[config]" << std::endl
            << "value = " << target << std::endl;

    // Check for the initialized value
    REQUIRE(config.checkValue(0));

    config.configure();
    // Check for the default initialized value
    REQUIRE(config.checkValue(-1));

    return;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(text.str()));

    Configurator::addStream(std::move(pcstream));

    config.configure();

    // Check for the configured value
    REQUIRE(config.checkValue(target));
}
} /* namespace Nextsim */
