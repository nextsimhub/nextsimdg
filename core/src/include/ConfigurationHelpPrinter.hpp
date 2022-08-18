/*!
 * @file ConfigurationHelpPrinter.hpp
 *
 * @date Aug 18, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CONFIGURATIONHELPPRINTER_HPP
#define CONFIGURATIONHELPPRINTER_HPP

#include "include/ConfigurationHelp.hpp"

#include <ostream>
#include <string>

static std::string b;

namespace Nextsim {

class ConfigurationHelpPrinter {
public:
    typedef ConfigurationHelp::ConfigType ConfigType;

    enum class Output {
        ANSI,
        MARKDOWN,
    };

    static void setOutput(Output out);

    static std::ostream& print(std::ostream& os, const ConfigurationHelp::HelpMap& map);
    static std::ostream& print(std::ostream& os, const ConfigurationHelp& help);

protected:
    static std::ostream& printString(std::ostream& os, const ConfigurationHelp& help);
    static std::ostream& printNumeric(std::ostream& os, const ConfigurationHelp& help);
    static std::ostream& printInteger(std::ostream& os, const ConfigurationHelp& help);
    static std::ostream& printModule(std::ostream& os, const ConfigurationHelp& help);

};

inline std::ostream& operator<<(std::ostream& os, const ConfigurationHelp& help)
{
    return ConfigurationHelpPrinter::print(os, help);
}
} /* namespace Nextsim */

#endif /* CONFIGURATIONHELPPRINTER_HPP */
