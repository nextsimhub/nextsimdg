/*!
 * @file ModelState.hpp
 *
 * @date Feb 28, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODELSTATE_HPP
#define MODELSTATE_HPP

#include "ModelArray.hpp"

#include "include/ConfigMap.hpp"

#include <map>
#include <string>

namespace Nextsim {

struct ModelState {
private:
public:
    typedef std::map<std::string, ModelArray> DataMap;

    DataMap data;
    ConfigMap config;

    inline void merge(ModelState& source) {
        data.merge(source.data);
        config.merge(source.config);
    }

    inline void merge(ModelState&& source) {
        data.merge(source.data);
        config.merge(source.config);
    }
};

} /* namespace Nextsim */

#endif /* MODELSTATE_HPP */
