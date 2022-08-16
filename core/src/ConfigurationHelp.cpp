/*!
 * @file ConfigurationHelp.cpp
 *
 * @date Aug 12, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ConfigurationHelp.hpp"

namespace Nextsim {

std::ostream& ConfigurationHelp::print(std::ostream& os) const
{
    switch (type) {
    case (ConfigType::STRING):
        return printString(os);
        break;
    case (ConfigType::NUMERIC):
        return printNumeric(os);
        break;
    case (ConfigType::INTEGER):
        return printInteger(os);
        break;
    case (ConfigType::MODULE):
        return printModule(os);
        break;
    default:
        return os;
    }
}

std::ostream& ConfigurationHelp::printString(std::ostream& os) const
{
    os << "\033[4m" << name << "\033[m" << std::endl;
    os << "\033[3mstring\033[m";
    if (!defaultValue.empty()) {
        os << "      default = " << defaultValue;
    }
    os << std::endl;
    os << text << std::endl;
    return os;
}

std::ostream& ConfigurationHelp::printNumeric(std::ostream& os) const
{
    os << "\033[4m" << name << "\033[m" << std::endl;
    os << "numeric:    " << range[0] << "—" << range[1] << " " << units;
    if (!defaultValue.empty()) {
        os << "      default = " << defaultValue << " " << units;
    }
    os << std::endl;
    os << text << std::endl;
    return os;
}

std::ostream& ConfigurationHelp::printInteger(std::ostream& os) const
{
    os << name << std::endl;
    os << "integer:    " << range[0] << "—" << range[1] << " " << units;
    if (!defaultValue.empty()) {
        os << "      default = " << defaultValue << " " << units;
    }
    os << std::endl;
    os << text << std::endl;
    return os;
}

std::ostream& ConfigurationHelp::printModule(std::ostream& os) const
{
    os << name << std::endl;
    os << "module: [";
    for (auto impl : range) {
        os << impl << ",";
    }
    os << "]";
    os << text << std::endl;
    return os;
}

} /* namespace Nextsim */
