/*!
 * @file ElementData.hpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_ELEMENTDATA_HPP
#define SRC_INCLUDE_ELEMENTDATA_HPP

#include "PhysicsData.hpp"
#include "PrognosticData.hpp"
#include "ExternalData.hpp"

namespace Nextsim {

class ElementData: public PrognosticData,
                   public PhysicsData,
                   public ExternalData {
public:
    ElementData();
    virtual ~ElementData();

    using PrognosticData::operator=;
    using PhysicsData::operator=;
    using ExternalData::operator=;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_ELEMENTDATA_HPP */
