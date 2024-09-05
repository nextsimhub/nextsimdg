/*!
 * @file HiblerSpread.cpp
 *
 * @date Apr 5, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/HiblerSpread.hpp"

namespace Nextsim {

double HiblerSpread::h0 = 0;
double HiblerSpread::phiM = 0;

static const double h0Default = 0.25;
static const double phimDefault = 0.5;

static const std::map<int, std::string> localKeyMap = {
    { HiblerSpread::H0_KEY, "Hibler.h0" },
    { HiblerSpread::PHIM_KEY, "Hibler.phiM" },
};

void HiblerSpread::configure()
{
    h0 = Configured::getConfiguration(localKeyMap.at(H0_KEY), h0Default);
    phiM = Configured::getConfiguration(localKeyMap.at(PHIM_KEY), phimDefault);
}

ModelState HiblerSpread::getStateRecursive(const OutputSpec& os) const
{
    return { {},
        {
            { localKeyMap.at(H0_KEY), h0 },
            { localKeyMap.at(PHIM_KEY), phiM },
        } };
}

HiblerSpread::HelpMap& HiblerSpread::getHelpText(HelpMap& map, bool getAll)
{
    map["HiblerSpread"] = {
        { localKeyMap.at(H0_KEY), ConfigType::NUMERIC, { "0", "∞" }, std::to_string(h0Default), "m",
            "The thickness of newly frozen ice." },
        { localKeyMap.at(PHIM_KEY), ConfigType::NUMERIC, { "0", "∞" }, std::to_string(phimDefault), "",
            "Power-law exponent for melting ice." },
    };
    return map;
}
HiblerSpread::HelpMap& HiblerSpread::getHelpRecursive(HelpMap& map, bool getAll)
{
    return getHelpText(map, getAll);
}

void HiblerSpread::freeze(const TimestepTime& tstep, double hice, double hsnow, double deltaHi,
    double newIce, double& cice, double& qow, double& deltaCfreeze)
{
    static const double ooh0 = 1. / h0;
    deltaCfreeze = newIce * ooh0;
}

void HiblerSpread::melt(const TimestepTime& tstep, double hice, double hsnow, double deltaHi,
    double& cice, double& qow, double& deltaCmelt)
{
    if (cice < 1) {
        deltaCmelt = deltaHi * cice * phiM / hice;
    }
}

}
