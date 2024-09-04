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
    typedef std::function<void()> FinalFn;

    /*!
     * Adds a function to be called at finalization. Functions are ordered last-in, first out.
     * @param fn The function to be added. Must have void() signature.
     */
    inline static void atfinal(const FinalFn& fn)
    {
        auto& fnSet = functions();
        fnSet.push_front(fn);
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
    inline static bool atfinalUnique(const FinalFn& fn)
    {
        if (contains(fn)) return false;
        atfinal(fn);
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
        auto& fnSet = functions();
        for (const auto& stored : fnSet) {
            if (stored.target<void()>() == fn.target<void()>()) return true;
        }
        return false;
    }

    /*!
     * @brief Runs the finalize functions without erasing all the functions.
     */
    inline static void run()
    {
        auto& fns = functions();
        for (auto& fn : fns) {
            fn();
        }
    }

    /*!
     * Runs the finalize functions and erases them.
     */
    inline static void finalize()
    {
        run();
        clear();
    }

    /*!
     * Returns the number of finalize functions that will be run.
     * @return the number of finalize functions that will be run.
     */
    inline static size_t count() { return functions().size(); }

    /*!
     * Eliminates all functions to be run at finalization.
     */
    inline static void clear() { functions().clear(); }
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
