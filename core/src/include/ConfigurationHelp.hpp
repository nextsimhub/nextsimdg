/*!
 * @file ConfigurationHelp.hpp
 *
 * @date Aug 12, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CONFIGURATIONHELP_HPP
#define CONFIGURATIONHELP_HPP

#include <string>
#include <vector>
#include <ostream>

namespace Nextsim {

class ConfigurationHelp {
public:

    enum class ConfigType {
        STRING,
        NUMERIC,
        INTEGER,
        MODULE
    };

    std::ostream& print(std::ostream& os) const;

    std::string name;
    ConfigType type;
    std::vector<std::string> range;
    std::string defaultValue;
    std::string units;
    std::string text;
protected:
    std::ostream& printString(std::ostream& os) const;
    std::ostream& printNumeric(std::ostream& os) const;
    std::ostream& printInteger(std::ostream& os) const;
    std::ostream& printModule(std::ostream& os) const;
private:
    // Common printing functions
};

inline std::ostream& operator<<(std::ostream& os, const ConfigurationHelp& ch) { return ch.print(os); }
} /* namespace Nextsim */

#endif /* CONFIGURATIONHELP_HPP */
