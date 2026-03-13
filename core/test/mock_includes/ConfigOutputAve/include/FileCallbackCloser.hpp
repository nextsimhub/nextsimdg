/*!
 * @file FileCallbackCloser.hpp
 *
 * @date Mar 13, 2026
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MOCK_INCLUDES_FILECALLBACKCLOSER_HPP
#define MOCK_INCLUDES_FILECALLBACKCLOSER_HPP

#include <string>

class FileCallbackCloser {
public:
    static void close(const std::string& s) { };
};

#endif /* MOCK_INCLUDES_FILECALLBACKCLOSER_HPP */
