/*!
 * @file Finalizer.cpp
 *
 * @date Sep 11, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/Finalizer.hpp"

namespace Nextsim {

void Finalizer::registerFunc(const FinalFn& fn)
{
    auto& fns = functions();
    fns.push_back(fn);
}

bool Finalizer::registerUnique(const FinalFn& fn)
{
    if (contains(fn))
        return false;
    registerFunc(fn);
    return true;
}

bool Finalizer::contains(const FinalFn& fn)
{
    auto& fns = functions();
    for (const auto& stored : fns) {
        // Compare function addresses
        if (*stored.target<FnType*>() == *fn.target<FnType*>())
            return true;
    }
    return false;
}

void Finalizer::finalize()
{
    while (!functions().empty()) {
        try {
            // Execute the last function in the list
            functions().back()();
        } catch (const std::exception& e) {
            // Remove the function from the finalization list even if an exception was thrown.
            functions().pop_back();
            throw;
        }
        // Remove the function from the finalization list.
        functions().pop_back();
    }
}

size_t Finalizer::count() { return functions().size(); }

Finalizer::FinalFnContainer& Finalizer::functions()
{
    static FinalFnContainer set;
    return set;
}

}
