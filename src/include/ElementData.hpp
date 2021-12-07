/*!
 * @file ElementData.hpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_ELEMENTDATA_HPP
#define SRC_INCLUDE_ELEMENTDATA_HPP

#include "ExternalData.hpp"
#include "PrognosticData.hpp"
#include "include/PhysicsData.hpp"

namespace Nextsim {

class UnusedData : public BaseElementData {
};

template <class Phys>
class ElementData : public PrognosticData,
                    public PhysicsData,
                    public ExternalData,
                    public Phys,
                    public UnusedData,
                    public Configured<ElementData<Phys>> {
public:
    ElementData() = default;
    ~ElementData() = default;

    using PrognosticData::operator=;
    using PhysicsData::operator=;
    using ExternalData::operator=;
    using Phys::operator=;

    void configure() override
    {
        PrognosticData::configure();
        Phys::configure();
    }
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_ELEMENTDATA_HPP */
