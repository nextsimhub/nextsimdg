/*!
 * @file Configurator.hpp
 *
 * @brief The header file for the Configurator class.
 *
 * @date Oct 8, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CONFIGURATOR_HPP
#define CONFIGURATOR_HPP

#include "include/ConfigMap.hpp"

#include <boost/program_options.hpp>
#include <fstream>
#include <istream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace Nextsim {

class NoAdditionalConfiguration;

/*!
 * @brief A class to handle the sources of configuration, both files and the command line.
 *
 * @details If an option is configured twice the value of the option will not
 * be updated. Whatever is parsed first sets the value of that option. This
 * allows the command line to override values in config files, as it is always
 * parsed first
 */
class Configurator {
public:
    Configurator() = default;
    virtual ~Configurator() = default;

    /*!
     *  Adds a config file to the configuration sources
     *
     *  @param filename the name of the file to be read.
     */

    inline static void addFile(const std::string& filename)
    {
        addStream(std::unique_ptr<std::istream>(new std::fstream(filename)));
    }

    /*!
     * @brief Adds several config files to the configuration sources
     *
     * @details Takes a container of the names of files to be used as
     * configuration sources. The individual filenames should be stored as
     * std::strings.
     *
     * @param container an iterable container holding std::string filenames.
     */
    template <typename C> static void addFiles(const C& container)
    {
        for (auto& filename : container)
            addFile(filename);
    }
    /*!
     * @brief Adds a istream source of configuration data.
     *
     * @param pis a std::unique_ptr to a std::istream containing the config data.
     */
    inline static void addStream(std::unique_ptr<std::istream> pis)
    {
        sources.push_back(std::move(pis));
    }
    /*!
     * @brief Adds several istream sources of configuration data.
     *
     * @details The container should hold std::unique_ptrs to std::istreams holding the
     * data.
     *
     * @param container an iterable container of std::unique_ptrs to
     * std::istream data sources.
     */
    template <typename C> static void addStreams(const C& container)
    {
        for (auto& stream : container) {
            addStream(stream);
        }
    }

    /*!
     * @brief Adds a std::stringstream configuration stream, wrapping all the pointer mechanics.
     * @param sstream The std::stringstream containiong the configuration text to be added.
     */
    static void addSStream(const std::stringstream& sstream);

    /*!
     * Removes previously assigned stream data sources, both files and istreams.
     */
    inline static void clearStreams() { sources.clear(); }

    /*!
     * Removes all data sources, both streams and command line.
     */
    inline static void clear()
    {
        clearStreams();
        setCommandLine(0, nullptr);
    }
    /*!
     * @brief Sets the command line data to be parsed.
     *
     * @details The data is formatted as the C standard argc and argv values.
     * Any values defined here will override the corresponding values that
     * might be found in the config files.
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

    /*!
     * @brief Gets the additional configuration according to the supplied implementation.
     *
     * @param source the (implementation-dependent) source of additional configuration.
     */
    static void getAdditionalConfiguration(const std::string& source);

    /*!
     * @brief Parses all configuration sources.
     *
     * @details Parses all the the stored configuration sources for the
     * configuration options specified in the options description. The command
     * line options are parsed first. Subsequent matching options will not
     * update the value of the option, so whatever is parsed first sets the
     * value of that option. This allows the command line to override values in
     * config files.
     *
     * @param opt An instance of boost::program_options describing the options
     * to be configured.
     */
    static boost::program_options::variables_map parse(
        const boost::program_options::options_description& opt);

private:
    static std::vector<std::unique_ptr<std::istream>> sources;

    static int m_argc;
    static char** m_argv;

public:
    class AdditionalConfiguration {
    public:
        virtual ~AdditionalConfiguration() = default;
        /*!
         * @brief Reads the additional configuration from the provided source,
         * which is interpreted in an implementation defined way.
         * @param source std::string containing the source specification.
         */
        virtual std::stringstream read(const std::string& source) = 0;
    };
    /*!
     * @brief Sets the source of any additional configuration data, such as a
     * netCDF restart file.
     *
     * @param pAC a pointer to the object that will provide the read() function to
     * read any additional configuration.
     */
    static void setAdditionalConfiguration(AdditionalConfiguration* pAC);

private:
    static AdditionalConfiguration* p_addConf;
};

//! A default implementation of Configurator::AdditionalConfiguration
class NoAdditionalConfiguration : public Configurator::AdditionalConfiguration {
    std::stringstream read(const std::string& source) override { return std::stringstream(); }
};

} /* namespace Nextsim */

#endif /* CONFIGURATOR_HPP */
