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
    typedef std::map<std::string, ModelArray> DataMap;
    typedef std::function<const ModelArray&()> DataReference;
    typedef std::map<std::string, DataReference> ReferenceMap;

    DataMap data;
    ConfigMap config;
    ReferenceMap references;

    void addReference(const std::string& name, const ModelArray& array)
    {
        references.emplace(name, [&array]() -> const ModelArray& { return array; });
    }

    template <typename Accessor>
    void addReference(const std::string& name, const Accessor& accessor)
    {
        references.emplace(
            name, [&accessor]() -> const ModelArray& { return accessor.getHostRO(); });
    }

    template <typename Fn> void forEachDataReference(Fn&& fn) const
    {
        for (const auto& [name, array] : data) {
            const ModelArray* arrayPtr = &array;
            fn(name, [arrayPtr]() -> const ModelArray& { return *arrayPtr; });
        }
        for (const auto& [name, array] : references) {
            if (data.count(name) == 0) {
                fn(name, array);
            }
        }
    }

    template <typename Fn> void forEachData(Fn&& fn) const
    {
        forEachDataReference(
            [&fn](const std::string& name, const DataReference& array) { fn(name, array()); });
    }

    void materializeReferences()
    {
        for (const auto& [name, array] : references) {
            if (data.count(name) == 0) {
                data.emplace(name, array());
            }
        }
        references.clear();
    }

    inline ModelState& merge(ModelState& source)
    {
        data.merge(source.data);
        references.merge(source.references);
        config.merge(source.config);
        return *this;
    }

    inline ModelState& merge(ModelState&& source)
    {
        data.merge(std::move(source.data));
        references.merge(std::move(source.references));
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

} /* namespace Nextsim */

#endif /* MODELSTATE_HPP */
