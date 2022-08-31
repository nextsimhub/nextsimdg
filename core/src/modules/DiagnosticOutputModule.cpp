/*!
 * @file DiagnosticOutputModule.cpp
 *
 * @date May 25, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/DiagnosticOutputModule.hpp"

#include "include/SimpleOutput.hpp"
#include "include/ConfigOutput.hpp"

#include <string>

namespace Module {
const std::string SIMPLEOUTPUT = "Nextsim::SimpleOutput";
const std::string CONFIGOUTPUT = "Nextsim::ConfigOutput";

template <>
Module<Nextsim::IDiagnosticOutput>::map Module<Nextsim::IDiagnosticOutput>::functionMap = {
    {SIMPLEOUTPUT, newImpl<Nextsim::IDiagnosticOutput, Nextsim::SimpleOutput>},
    {CONFIGOUTPUT, newImpl<Nextsim::IDiagnosticOutput, Nextsim::ConfigOutput>},
};

template <>
Module<Nextsim::IDiagnosticOutput>::fn Module<Nextsim::IDiagnosticOutput>::spf = functionMap.at(SIMPLEOUTPUT);
template <>
std::unique_ptr<Nextsim::IDiagnosticOutput> Module<Nextsim::IDiagnosticOutput>::staticInstance
= std::move(newImpl<Nextsim::IDiagnosticOutput, Nextsim::SimpleOutput>());

template <>
std::string Module<Nextsim::IDiagnosticOutput>::moduleName(){    return "Nextsim::IDiagnosticOutput";}

template <>
Nextsim::IDiagnosticOutput& getImplementation<Nextsim::IDiagnosticOutput>()
{
    return getImplTemplate<Nextsim::IDiagnosticOutput, DiagnosticOutputModule>();
}
template <>
void setImplementation<Nextsim::IDiagnosticOutput>(const std::string& implName)
{
    setImplTemplate<DiagnosticOutputModule>(implName);
}
template <>
std::unique_ptr<Nextsim::IDiagnosticOutput> getInstance()
{
    return getInstTemplate<Nextsim::IDiagnosticOutput, DiagnosticOutputModule>();
}
DiagnosticOutputModule::Constructor DiagnosticOutputModule::ctor;
DiagnosticOutputModule::Constructor::Constructor()
{
    addToConfiguredModules<Nextsim::IDiagnosticOutput, DiagnosticOutputModule>();
}

} /* namespace Module */
