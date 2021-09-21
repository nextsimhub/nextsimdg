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

class UnusedData: public BaseElementData {

};

template<class Phys>
class ElementData: public PrognosticData,
                   public PhysicsData,
                   public ExternalData,
                   public Phys,
                   public UnusedData {
public:
    ElementData() = default;
    ~ElementData() = default;

    using PrognosticData::operator=;
    using PhysicsData::operator=;
    using ExternalData::operator=;
    using Phys::operator=;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_ELEMENTDATA_HPP */
