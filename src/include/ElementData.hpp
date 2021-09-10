/*!
 * @file ElementData.hpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_ELEMENTDATA_HPP_
#define SRC_INCLUDE_ELEMENTDATA_HPP_

#include <PhysicsData.hpp>
#include <PrognosticData.hpp>

namespace Nextsim {

class ElementData: public PrognosticData, public PhysicsData {
public:
    ElementData();
    virtual ~ElementData();
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_ELEMENTDATA_HPP_ */
