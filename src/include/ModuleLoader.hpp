/*
 * @file ModuleLoader.hpp
 *
 * @date Sep 23, 2021
 * @author Tim Spain
 */

#ifndef SRC_INCLUDE_MODULELOADER_HPP
#define SRC_INCLUDE_MODULELOADER_HPP

#include <boost/program_options.hpp>
#include <map>
#include <memory>
#include <set>
#include <string>

class ModuleLoader {
public:
    static ModuleLoader& getLoader()
    {
        static ModuleLoader instance; // C++11 magic static
        instance.init();
        return instance;
    }

    // typedef boost::program_options::variables_map VariablesMap;
    typedef std::map<std::string, std::string> VariablesMap;

    void init();
    void init(const VariablesMap&);
    inline const std::set<std::string>& listModules() const { return m_modules; }
    inline const std::set<std::string>& listImplementations(const std::string& module) const
    {
        return m_availableImplementationNames.at(module);
    }
    template <class T> std::unique_ptr<T> getInstance() const;

    template <class T> T& getImplementation();

    void setImplementation(const std::string& module, const std::string& impl);
    // Singleton function definitions
private:
    ModuleLoader() {};

public:
    ModuleLoader(const ModuleLoader&) = delete;
    void operator=(const ModuleLoader&) = delete;

private:
    // Make sure initialization happens
    bool isInit = false;
    // One module could have many names (but probably shouldn't)
    std::set<std::string> m_modules;
    // Names of available implementations
    std::map<std::string, std::set<std::string>> m_availableImplementationNames;
};

#endif /* SRC_INCLUDE_MODULELOADER_HPP */
