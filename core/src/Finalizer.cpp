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
            if (stored == fn)
                return true;
        }
        return false;
    }

void Finalizer::finalize()
    {
        while (!functions().empty()) {
            auto fn = functions().back();
            // Remove the function from the finalization list before executing it.
            functions().pop_back();
            fn();
        }
    }

size_t Finalizer::count() { return functions().size(); }

Finalizer::FinalFnContainer& Finalizer::functions()
    {
        static FinalFnContainer set;
        return set;
    }

}
