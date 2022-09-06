/*!
 * @file ConfigurationHelp.hpp
 *
 * @date Aug 12, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CONFIGURATIONHELP_HPP
#define CONFIGURATIONHELP_HPP

#include <list>
#include <map>
#include <string>
#include <vector>

namespace Nextsim {

class ConfigurationHelp {
public:
    typedef std::list<ConfigurationHelp> OptionList;
    typedef std::map<std::string, OptionList> HelpMap;

    enum class ConfigType { STRING, NUMERIC, INTEGER, MODULE };

    std::string name;
    ConfigType type;
    std::vector<std::string> range;
    std::string defaultValue;
    std::string units;
    std::string text;
};

} /* namespace Nextsim */

#endif /* CONFIGURATIONHELP_HPP */
