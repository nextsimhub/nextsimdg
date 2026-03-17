/*!
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MOCK_INCLUDES_IDIAGNOSTICOUTPUT_HPP
#define MOCK_INCLUDES_IDIAGNOSTICOUTPUT_HPP

#include <map>

#include "include/ConfigurationHelp.hpp"
#include "include/ModelState.hpp"
#include "include/Time.hpp"

class IDiagnosticOutput {
public:
    virtual ~IDiagnosticOutput() = default;

    virtual void setFilenamePrefix(const std::string& filePrefix) = 0;
    virtual void outputState(const Nextsim::ModelState& state) = 0;
    virtual void setData(const Nextsim::ModelState::DataMap& state) { }
    virtual void setModelStart(const Nextsim::TimePoint& ModelStart) { }

    // functions usually declared in ModelComponent
    virtual std::string getName() const = 0;
    const static std::map<std::string, std::string> externalNames;
};
typedef Nextsim::ConfigurationHelp::HelpMap HelpMap;
typedef Nextsim::ConfigurationHelp::ConfigType ConfigType;

class ModelMetadata {
public:
    Nextsim::TimePoint& time()
    {
        static Nextsim::TimePoint m_time;
        return m_time;
    };
    void affixCoordinates(Nextsim::ModelState& ms) { };
    static ModelMetadata& getInstance()
    {
        static ModelMetadata meta;
        return meta;
    };
};

class ModelComponent {
    class Store {
    public:
        using StoreMap = std::map<std::string, Nextsim::ModelArray*>;
        StoreMap data;
        StoreMap getAllData() { return data; }
    };
public:
    static Store& getStore()
    {
        static Store store;
        return store;
    }
    void setData(const std::string& name, Nextsim::ModelArray* p_data)
    {
        getStore().data[name] = p_data;
    }
};

#endif /* MOCK_INCLUDES_IDIAGNOSTICOUTPUT_HPP */
