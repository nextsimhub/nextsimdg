/*!
 * @file CommandLineParser.hpp
 *
 * @date Oct 8, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_COMMANDLINEPARSER_HPP
#define SRC_INCLUDE_COMMANDLINEPARSER_HPP

#include <vector>
#include <string>
#include <boost/program_options.hpp>

namespace Nextsim {

class CommandLineParser {
public:
    CommandLineParser(int argc, char* argv[]);
    virtual ~CommandLineParser() = default;

    std::vector<std::string> getConfigFileNames();
private:
    CommandLineParser() = default;

    boost::program_options::variables_map m_arguments;
    std::vector<std::string> m_configFilenames;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_COMMANDLINEPARSER_HPP */
