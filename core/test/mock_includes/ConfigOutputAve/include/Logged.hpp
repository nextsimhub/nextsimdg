/*!
 * @file Logged.hpp
 *
 * @date Mar 13, 2026
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MOCK_INCLUDES_LOGGED_HPP
#define MOCK_INCLUDES_LOGGED_HPP

#include <string>

class Logged {
public:
    static void info(const std::string& s) { }
    static void warning(const std::string& s) { }

};

#endif /* MOCK_INCLUDES_LOGGED_HPP */
