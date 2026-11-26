/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODELSTATE_HPP
#define MODELSTATE_HPP

#include "ModelArray.hpp"

#include "include/ConfigMap.hpp"

#include <functional>
#include <map>
#include <string>
#include <utility>

namespace Nextsim {

struct ModelState {
private:
public:
    typedef std::map<std::string, ModelArray> DataMap;

    DataMap data;
    ConfigMap config;

    inline ModelState& merge(ModelState& source)
    {
        data.merge(source.data);
        config.merge(source.config);
        return *this;
    }

    inline ModelState& merge(ModelState&& source)
    {
        data.merge(std::move(source.data));
        config.merge(std::move(source.config));
        return *this;
    }

    inline ModelState& merge(ConfigMap& configSrc)
    {
        config.merge(configSrc);
        return *this;
    }

    inline ModelState& merge(ConfigMap&& configSrc)
    {
        config.merge(std::move(configSrc));
        return *this;
    }

    inline ModelState& merge(DataMap& dataSrc)
    {
        data.merge(dataSrc);
        return *this;
    }

    inline ModelState& merge(DataMap&& dataSrc)
    {
        data.merge(std::move(dataSrc));
        return *this;
    }
};

struct ModelStateView {
    using DataMap = std::map<std::string, std::reference_wrapper<const ModelArray>>;

    ModelStateView() = default;
    explicit ModelStateView(const ModelState& state)
        : config(state.config)
    {
        for (const auto& [name, array] : state.data) {
            addReference(name, array);
        }
    }

    ModelStateView(const ModelStateView&) = delete;
    ModelStateView& operator=(const ModelStateView&) = delete;

    void addReference(const std::string& name, const ModelArray& array)
    {
        data.emplace(name, std::cref(array));
    }

    void merge(ModelState&& state)
    {
        config.merge(std::move(state.config));
        for (auto& [name, array] : state.data) {
            auto [it, inserted] = ownedData.emplace(name, std::move(array));
            if (inserted) {
                addReference(it->first, it->second);
            }
        }
    }

    DataMap data;
    ConfigMap config;

private:
    std::map<std::string, ModelArray> ownedData;
};

} /* namespace Nextsim */

#endif /* MODELSTATE_HPP */
