/*!
 * @file ConfigurationHelp.hpp
 *
 * @date 20 Nov 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CONFIGURATIONHELP_HPP
#define CONFIGURATIONHELP_HPP

#include <iomanip>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace Nextsim {

class ConfigurationHelp {
public:
    typedef std::list<ConfigurationHelp> OptionList;
    typedef std::map<std::string, OptionList> HelpMap;

    enum class ConfigType { STRING, NUMERIC, INTEGER, MODULE, BOOLEAN };

    std::string name;
    ConfigType type;
    std::vector<std::string> range;
    std::string defaultValue;
    std::string units;
    std::string text;

    // Ever so slightly better formatting than std::to_string
    template <typename T> static std::string toString(const T input, const int n = 6)
    {
        std::ostringstream output;
        output << std::setprecision(n) << input;
        return output.str();
    }
};

} /* namespace Nextsim */

#endif /* CONFIGURATIONHELP_HPP */
