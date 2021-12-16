/*!
 * @file ModuleLoader.hpp
 *
 * @date Sep 23, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_MODULELOADER_HPP
#define SRC_INCLUDE_MODULELOADER_HPP

#include <boost/program_options.hpp>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <string>

//! A class to manage run-time polymorphism within the model.
class ModuleLoader {
public:
    //! Returns the default loader instance
    static ModuleLoader& getLoader()
    {
        static ModuleLoader instance; // C++11 magic static
        instance.init();
        return instance;
    }

    typedef std::map<std::string, std::string> VariablesMap;

    //! Initializes the loader with no modules loaded.
    void init();
    /*!
     * @brief Initializes the model with an initial set of module
     * implementations selected.
     *
     * @param map The mapping between interface class names to the requested
     * implementation class names.
     */
    void init(const VariablesMap& map);

    //! Returns the set of modules that can be implemented as a set of strings.
    inline const std::set<std::string>& listModules() const { return m_modules; }
    /*!
     * Lists the available implementation names for a given module.
     *
     * @param module The name of the module of interest.
     */
    inline const std::list<std::string>& listImplementations(const std::string& module) const
    {
        return m_availableImplementationNames.at(module);
    }

    /*!
     * @brief Returns a newly created instance of the implementing class.
     *
     * @details The module of interest is specified as a template argument. A
     * new instance of the already-selected implementing class is returned
     * using a unique_ptr to the instance.
     */
    template <class T> std::unique_ptr<T> getInstance() const;

    /*!
     * @brief Returns a reference to a static instance of the implementing class.
     *
     * @details The class holds one instance of every implementing class. Once
     * the implementing class is set, this function will return a reference to
     * that instance.
     */
    template <class T> T& getImplementation();

    /*!
     * @brief Sets the implementation class of a module.
     *
     * @details Given a module, specified by its name, set the implementing
     * class by name. The names should match the name given in the module
     * specification file.
     *
     * @param module The fully qualified name of the module to be implemented.
     * @param impl The fully qualified name of the implementing class.
     */
    void setImplementation(const std::string& module, const std::string& impl);

    /*!
     * @brief Sets the default implementation for a module.
     *
     * @param module The module for which the default implementation is to be
     * loaded.
     */
    void setDefault(const std::string& module);

    //! Sets the default implementation for all modules.
    void setAllDefaults();

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
    std::map<std::string, std::list<std::string>> m_availableImplementationNames;
};

#endif /* SRC_INCLUDE_MODULELOADER_HPP */
