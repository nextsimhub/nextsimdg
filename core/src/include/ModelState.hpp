/*!
 * @file ModelState.hpp
 *
 * @date Feb 28, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODELSTATE_HPP
#define MODELSTATE_HPP

#include "ModelArray.hpp"

#include <boost/program_options.hpp>
#include <map>
#include <string>

namespace Nextsim {

//typedef std::map<std::string, ModelArray> ModelState;
struct ModelState {
public:
    typedef std::map<std::string, ModelArray> DataMap;
    typedef boost::program_options::variables_map ConfigMap;

    DataMap data;
    ConfigMap config;

};

} /* namespace Nextsim */

#endif /* MODELSTATE_HPP */
