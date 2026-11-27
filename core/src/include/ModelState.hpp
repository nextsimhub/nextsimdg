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
    //! Callback returning a read-only array; invoking it may synchronize host data.
    typedef std::function<const ModelArray&()> DataReference;
    typedef std::map<std::string, DataReference> ReferenceMap;

    DataMap data;
    ConfigMap config;
    ReferenceMap references;

    /*!
     * @brief Adds a deferred reference to an existing array.
     *
     * The array is not copied. It must outlive this ModelState's use of the
     * reference. The reference is evaluated only when data is requested.
     *
     * @param name Name used to identify the array in this state.
     * @param array Array to reference.
     */
    void addReference(const std::string& name, const ModelArray& array)
    {
        references.emplace(name, [&array]() -> const ModelArray& { return array; });
    }

    /*!
     * @brief Adds a deferred reference through a model-array accessor.
     *
     * The accessor is not queried until the reference is evaluated, avoiding
     * an early host-device synchronization. The accessor must outlive this
     * ModelState's use of the reference.
     *
     * @tparam Accessor Type providing a const getHostRO() method.
     * @param name Name used to identify the array in this state.
     * @param accessor Accessor that supplies the array on demand.
     */
    template <typename Accessor>
    void addReference(const std::string& name, const Accessor& accessor)
    {
        references.emplace(
            name, [&accessor]() -> const ModelArray& { return accessor.getHostRO(); });
    }

    /*!
     * @brief Iterates over owned data and deferred data references.
     *
     * Owned data takes precedence when a name exists in both maps. The
     * callback receives the name and a deferred getter; invoking the getter
     * resolves the corresponding array and may synchronize it to the host.
     *
     * @param fn Callable accepting (const std::string&, const DataReference&).
     */
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

    /*!
     * @brief Iterates over all arrays, resolving deferred references.
     *
     * Owned data takes precedence when a name exists in both maps. Deferred
     * references are resolved immediately before their array is passed to the
     * callback.
     *
     * @param fn Callable accepting (const std::string&, const ModelArray&).
     */
    template <typename Fn> void forEachData(Fn&& fn) const
    {
        forEachDataReference(
            [&fn](const std::string& name, const DataReference& array) { fn(name, array()); });
    }

    /*!
     * @brief Copies deferred arrays into owned data and removes the references.
     *
     * Existing owned data takes precedence for duplicate names. Every
     * remaining reference is evaluated, which may synchronize device data to
     * the host.
     */
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
