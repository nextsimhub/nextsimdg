/*!
 * @file DiagnosticOutputModule.hpp
 *
 * @date May 25, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef DIAGNOSTICOUTPUTMODULE_HPP
#define DIAGNOSTICOUTPUTMODULE_HPP

#include "include/Module.hpp"

#include "include/IDiagnosticOutput.hpp"

namespace Module {

template <> Module<Nextsim::IDiagnosticOutput>::map Module<Nextsim::IDiagnosticOutput>::functionMap;
class DiagnosticOutputModule : public Module<Nextsim::IDiagnosticOutput> {
    struct Constructor {
        Constructor();
    };
    static Constructor ctor;
};

} /* namespace Module */

#endif /* DIAGNOSTICOUTPUTMODULE_HPP */
