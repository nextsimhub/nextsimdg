/*!
 * @file ModelFinalizer.hpp
 *
 * @date Sep 3, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef FINALIZER_HPP
#define FINALIZER_HPP

#include <functional>
#include <list>
#include <set>
#include <vector>

#include <iostream> // FIXME remove me

namespace Nextsim {

class Finalizer {
public:
    using FinalFn = std::function<void()>;

    /*!
     * Adds a function to be called at finalization. Functions are ordered last-in, first out.
     * @param fn The function to be added. Must have void() signature.
     */
    inline static void registerFunc(const FinalFn& fn)
    {
        auto& fns = functions();
        fns.push_back(fn);
    }

    /*!
     * @brief Adds a function to be called at finalization unless it will already be called at
     * finalization.
     *
     * @details Adds a function in the same way as atfinal(), unless the function already exists in
     * the list of functions to be called at finalization. Identity of functions is equality of the
     * values of their target<void()>() member functions. Functions are ordered last-in, first out.
     * @param fn The function to be added. Must have void() signature.
     * @return true if the function was added, false if it already existed and was not added.
     */
    inline static bool registerUnique(const FinalFn& fn)
    {
        if (contains(fn))
            return false;
        registerFunc(fn);
        return true;
    }

    /*!
     * @brief Returns whether a function will already be called at finalization.
     *
     * @details If the supplied function will already be called at finalization, return true, else
     * return false. Identity of functions is equality of the values of their target<void()>()
     * member functions.
     * @param fn The function to be queried. Must have void() signature.
     * @return true is the function is contained in those to be called at finalization, false if it
     * is not.
     *
     */
    inline static bool contains(const FinalFn& fn)
    {
        auto& fns = functions();
        for (const auto& stored : fns) {
            if (stored.target<void()>() == fn.target<void()>())
                return true;
        }
        return false;
    }

    /*!
     * @brief Runs the finalize functions and erases them.
     *
     * @details In last-in, first-out order, the finalization functions are executed. If an
     * exception is thrown, the throwing function will have been removed from the finalization
     * functions, should execution be able to recover. If no exception is thrown, then there will
     * be no finalization functions remaining. Otherwise, re-running finalize will begin execution
     * with the function assigned previous to the one that threw the exception.
     */
    inline static void finalize()
    {
        while (!functions().empty()) {
            auto fn = functions().back();
            // Remove the function from the finalization list before executing it.
            functions().pop_back();
            fn();
        }
    }

    /*!
     * Returns the number of finalize functions that will be run.
     * @return the number of finalize functions that will be run.
     */
    inline static size_t count() { return functions().size(); }

    /*!
     * Eliminates all functions to be run at finalization.
     */
//    inline static void clear() { functions().clear(); }

private:
    typedef std::list<FinalFn> FinalFnContainer;
    inline static FinalFnContainer& functions()
    {
        static FinalFnContainer set;
        return set;
    }
};

} /* namespace Nextsim */

#endif /* FINALIZER_HPP */
