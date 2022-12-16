/*!
 * @file ConfigurationHelpPrinter.hpp
 *
 * @brief The header file for the ConfigurationHelpPrinter class.
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

/*!
 * A class to handle the printing of the ConfigurationHelp class to screen.
 */
class ConfigurationHelpPrinter {
public:
    typedef ConfigurationHelp::ConfigType ConfigType;

    enum class Output {
        ANSI,
        MARKDOWN,
    };

    static const std::string allStr;
    static const std::string availStr;

    /*!
     * @brief Sets the output style.
     *
     * @details As the details of the depend on whether the help information
     * should be formatted for text or screen, this function allow selection
     * between the two, with Output::ANSI utilising ANSI formatting codes for
     * terminal help and Output::MARKDOWN utilising markdown formatting for text.
     *
     * @param out The Output enum corresponding to the desired formatting.
     */
    static void setOutput(Output out);

    /*!
     * @brief Prints one or more items of help to the passed stream.
     *
     * @details Given a target stream, a ConfigurationMap::HelpMap and a target
     * name, print the relevant help to the stream in the format previously
     * requested by setting setOutput().
     *
     * @param os An output character stream.
     * @param map A ConfigurationMap::HelpMap containing the configuration help to be processed.
     * @param target A std::string naming the requested help item. In addition
     *              to a valid name, "all" and "avail" are also valid values,
     *              which print all configured and all available values, respectively.
     *
     * @returns A reference to the updated output stream.
     */
    static std::ostream& print(
        std::ostream& os, const ConfigurationHelp::HelpMap& map, const std::string& target);
    /*!
     * @brief Prints a specific help item to the passed stream.
     *
     * @details Given a target stream and a ConfigurationHelp item, print the
     * details of the item to the stream, using the format previous set using
     * setFormat().
     *
     * @param os An output character stream.
     * @param help The help item to be printed.
     *
     * @returns A reference to the updated output stream.
     */
    static std::ostream& print(std::ostream& os, const ConfigurationHelp& help);

protected:
    // Individual functions to print each possible type of configuration item.
    static std::ostream& printString(std::ostream& os, const ConfigurationHelp& help);
    static std::ostream& printNumeric(std::ostream& os, const ConfigurationHelp& help);
    static std::ostream& printInteger(std::ostream& os, const ConfigurationHelp& help);
    static std::ostream& printModule(std::ostream& os, const ConfigurationHelp& help);
};

/*!
 * Overloads the << output operator to work for ConfigurationHelp objects on
 * std::ostream output streams.
 *
 * @param os An output character stream to be written to.
 * @param help The help item to be printed.
 *
 * @returns A reference to the updated output stream.
 *
 */
inline std::ostream& operator<<(std::ostream& os, const ConfigurationHelp& help)
{
    return ConfigurationHelpPrinter::print(os, help);
}
} /* namespace Nextsim */

#endif /* CONFIGURATIONHELPPRINTER_HPP */
