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

namespace Nextsim {

class Finalizer {
public:
    using FinalFn = std::function<void()>;
    /*!
     * Adds a function to be called at finalization. Functions are ordered last-in, first out.
     * @param fn The function to be added. Must have void() signature.
     */
    static void registerFunc(const FinalFn& fn);

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
    static bool registerUnique(const FinalFn& fn);

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
    static bool contains(const FinalFn& fn);

    /*!
     * @brief Runs the finalize functions and erases them.
     *
     * @details In last-in, first-out order, the finalization functions are executed. If an
     * exception is thrown, the throwing function will have been removed from the finalization
     * functions, should execution be able to recover. If no exception is thrown, then there will
     * be no finalization functions remaining. Otherwise, re-running finalize will begin execution
     * with the function assigned previous to the one that threw the exception.
     */
    static void finalize();

    /*!
     * Returns the number of finalize functions that will be run.
     * @return the number of finalize functions that will be run.
     */
    static size_t count();

private:
    typedef std::list<FinalFn> FinalFnContainer;
    /*!
     * Returns a reference to a static container of the finalizer functions.
     */
    static FinalFnContainer& functions();
};

} /* namespace Nextsim */

#endif /* FINALIZER_HPP */
