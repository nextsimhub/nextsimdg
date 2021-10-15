/*!
 * @file Configured.cpp
 *
 * @date Oct 8, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/Configured.hpp"
#include "include/Configurator.hpp"

#include <set>
#include <memory>

namespace Nextsim {

std::set<Configured*> Configured::configuredObjects;

Configured::Configured()
{
    configuredObjects.insert(this);
}

void Configured::parse()
{
    vm = Configurator::parse(opt);
}

void Configured::addConfiguredObject(Configured* pcfg)
{
    configuredObjects.insert(pcfg);
}

void Configured::configureAll()
{
    for (Configured* pConfigured: configuredObjects) {
        pConfigured->parse();
    }
    configuredObjects.clear();
}
}
