/*!
 * @file ConfigurationHelpPrinter.cpp
 *
 * @date Aug 18, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ConfigurationHelpPrinter.hpp"

namespace Nextsim {

static std::string h1 = "";
static std::string h1End = "";
static std::string h2 = "";
static std::string h2End = "";
static std::string b = "";
static std::string u = "";
static std::string i = "";
static std::string bEnd = "";
static std::string uEnd = "";
static std::string iEnd = "";

static std::string section(const std::string& title)
{
    return h1 + title + h1End;
}

static std::string option(const std::string& name)
{
    return h2 + name + h2End;
}

static std::string defaultVal(const std::string& name)
{
    return b + name + bEnd;
}

static std::string type(const std::string& desc)
{
    return i + desc + iEnd;
}

std::ostream& ConfigurationHelpPrinter::print(std::ostream& os, const ConfigurationHelp::HelpMap& map)
{
    for (auto configEntry : map) {
        os << section(configEntry.first) << std::endl << std::endl;
        for (auto optionEntry : configEntry.second) {
            os << optionEntry << std::endl;
        }
        os << std::endl;
    }
    return os;
}


std::ostream& ConfigurationHelpPrinter::print(std::ostream& os, const ConfigurationHelp& help)
{
    switch (help.type) {
    case (ConfigType::STRING):
        return printString(os, help);
        break;
    case (ConfigType::NUMERIC):
        return printNumeric(os, help);
        break;
    case (ConfigType::INTEGER):
        return printInteger(os, help);
        break;
    case (ConfigType::MODULE):
        return printModule(os, help);
        break;
    default:
        return os;
    }
}

std::ostream& ConfigurationHelpPrinter::printString(std::ostream& os, const ConfigurationHelp& help)
{
    os << option(help.name) << std::endl;
    os << type("string");
    if (!help.defaultValue.empty()) {
        os << "           (default = " << help.defaultValue << ")";
    }
    os << std::endl;
    os << help.text << std::endl;
    return os;
}

std::ostream& ConfigurationHelpPrinter::printNumeric(std::ostream& os, const ConfigurationHelp& help)
{
    os << option(help.name) << std::endl;
    os << type("numeric") + "    range: " << help.range[0] << "—" << help.range[1] << " " << help.units;
    if (!help.defaultValue.empty()) {
        os << " (default = " << help.defaultValue << ")";
    }
    os << std::endl;
    os << help.text << std::endl;
    return os;
}

std::ostream& ConfigurationHelpPrinter::printInteger(std::ostream& os, const ConfigurationHelp& help)
{
    os << option(help.name) << std::endl;
    os << type("integer") + "     range: " << help.range[0] << "—" << help.range[1] << " " << help.units;
    if (!help.defaultValue.empty()) {
        os << " default = " << help.defaultValue << ")";
    }
    os << std::endl;
    os << help.text << std::endl;
    return os;
}

std::ostream& ConfigurationHelpPrinter::printModule(std::ostream& os, const ConfigurationHelp& help)
{
    os << option(help.name) << std::endl;
    os << type("module") + " [";

    os << defaultVal(help.defaultValue);
    for (auto impl : help.range) {
        if (impl != help.defaultValue) {
            os << ", " << impl;
        }
    }
    os << "]" << std::endl;
    os << help.text << std::endl;
    return os;
}

static std::string ansiMode(std::string mode)
{
    return "\033[" + mode + "m";
}

void ConfigurationHelpPrinter::setOutput(Output out)
{
    if (out == Output::ANSI) {
        std::string esc = "\033[";
        std::string cls = "m";
        b = ansiMode("1");
        u = ansiMode("4");
        i = ansiMode("3");
        h1 = u + b;
        h2 = u;
        bEnd = uEnd = iEnd = h1End = h2End = ansiMode("");
    } else if (out == Output::MARKDOWN) {
        h1 = "# ";
        h1End = "";
        h2 = "## ";
        h2End = "";
        b = "**";
        bEnd = b;
        i = "_";
        iEnd = i;
        // No pure underline in GitHub Markdown
    }
}

} /* namespace Nextsim */
