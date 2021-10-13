/*!
 * @file Configurator.hpp
 *
 * @date Oct 8, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_CONFIGURATOR_HPP
#define SRC_INCLUDE_CONFIGURATOR_HPP

#include <vector>
#include <memory>
#include <istream>
#include <string>
#include <boost/program_options.hpp>

namespace Nextsim {

/*!
 * A class to handle the sources of configuration, both files and the command line.
 */
class Configurator {
public:
    Configurator() = default;
    virtual ~Configurator() = default;

    /*!
     *  Add a config file to the configuration sources
     *
     *  @param filename the name of the file to be read.
     */

    inline static void addFile(const std::string& filename)
    {
        addStream(std::unique_ptr<std::istream>(new std::fstream(filename)));
    }

    /*!
     * Add several config files to the configuration sources
     *
     * Takes a container of the names of files to be used as configuration
     * sources. The individual filenames should be stored as std::strings.
     *
     * @param container an iterable container holding std::string filenames.
     */
    template<typename C>
    static void addFiles(const C& container)
    {
        for (auto& filename: container)
            addFile(filename);
    }
    /*!
     * Add a istream source of configuration data.
     *
     * @param pis a std::unique_ptr to a std::istream containing the config data.
     */
    inline static void addStream(std::unique_ptr<std::istream> pis)
    {
        sources.push_back(std::move(pis));
    }
    /*!
     * Add several istream sources of configuration data.
     *
     * The container should hold std::unique_ptrs to std::istreams holding the
     * data.
     *
     * @param container an iterable container of std::unique_ptrs to
     * std::istream data sources.
     */
    template<typename C>
    static void addStreams(const C& container)
    {
        for(auto& stream: container) {
            addStream(stream);
        }
    }

    /*!
     * Remove all previously assigned data sources, both files and istreams.
     */
    inline static void clearStreams()
    {
        sources.clear();
    }
    /*!
     * Set the command line data to be parsed.
     *
     * The data is formatted as the C standard argc and argv values. Any values
     * defined here will override the corresponding values that might be found
     * in the config files.
     *
     * @param argc the number of arguments to be parsed
     * @param argv an array of zero terminated character arrays making up the
     * command line arguments, with an addition null at argv[argc].
     */
    inline static void setCommandLine(int argc, char* argv[])
    {
        m_argc = argc;
        m_argv = argv;
    }
    static boost::program_options::variables_map parse(const boost::program_options::options_description& opt);

private:
    static std::vector<std::unique_ptr<std::istream>> sources;

    static int m_argc;
    static char** m_argv;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_CONFIGURATOR_HPP */
