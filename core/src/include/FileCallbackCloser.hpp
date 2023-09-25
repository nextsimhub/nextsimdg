/*!
 * @file FileCallbackCloser.hpp
 *
 * @date 15 May 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef FILECALLBACKCLOSER_HPP
#define FILECALLBACKCLOSER_HPP

#include <functional>
#include <list>
#include <string>

namespace Nextsim {

/*!
 * A class to provide callback-based file closure, to avoid passing filenames
 * though classes where they don't belong.
 */
class FileCallbackCloser {

public:
    typedef std::function<void(const std::string&)> ClosingFn;

    /*!
     * Adds a callback function which will be executed when the close function is called.
     * @param fn the function to be called on file closure
     */
    static void onClose(ClosingFn fn);

    /*!
     * Closes the given filename be calling all onClose callbacks.
     * @param filename the filename to be closed.
     */
    static void close(const std::string& fileName);

    /*!
     * Clears all onClose functions
     */
    static void clearAllClose();

private:
    static std::list<ClosingFn> closingFns;
};
} /* namespace Nextsim */

#endif /* FILECALLBACKCLOSER_HPP */
