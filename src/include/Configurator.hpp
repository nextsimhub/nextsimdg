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

class Configurator {
public:
    Configurator() = default;
    virtual ~Configurator() = default;

    static inline void addStream(std::unique_ptr<std::istream> pis)
    {
        sources.push_back(std::move(pis));
    }
    template<typename C>
    static void addStreams(const C& container)
    {
        for(auto& stream: container) {
            addStream(stream);
        }
    }

    static inline void clearStreams()
    {
        sources.clear();
    }
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
