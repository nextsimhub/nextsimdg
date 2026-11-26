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
    typedef std::map<std::string, std::reference_wrapper<const ModelArray>> ReferenceMap;

    DataMap data;
    ConfigMap config;
    ReferenceMap references;

    void addReference(const std::string& name, const ModelArray& array)
    {
        references.emplace(name, std::cref(array));
    }

    template <typename Fn> void forEachData(Fn&& fn) const
    {
        for (const auto& [name, array] : data) {
            fn(name, array);
        }
        for (const auto& [name, array] : references) {
            if (data.count(name) == 0) {
                fn(name, array.get());
            }
        }
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
