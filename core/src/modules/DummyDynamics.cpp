/*!
 * @file DummyDynamics.cpp
 *
 * @date 16 Feb 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/DummyDynamics.hpp"

#include "include/gridNames.hpp"

#include <string>
#include <vector>


namespace Nextsim {

static const std::vector<std::string> namedFields = { hiceName, ciceName, uName, vName };
DummyDynamics::DummyDynamics()
    : IDynamics()
{
}

void DummyDynamics::setData(const ModelState::DataMap& ms)
{
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    IDynamics::setData(ms);
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;

    // Set the data in the kernel arrays.
    for (const auto& fieldName : namedFields) {
        std::cerr << __FILE__ << ":" << __LINE__ << " for " << fieldName << std::endl;
        kernel.setData(fieldName, ms.at(fieldName));
    }
}

void DummyDynamics::update(const TimestepTime& tst)
{
    kernel.setData(hiceName, hice.data());
    kernel.setData(ciceName, cice.data());
    kernel.setData(uName, uice);
    kernel.setData(vName, vice);

    kernel.update(tst);

    hice.data() = kernel.getDG0Data(hiceName);
    cice.data() = kernel.getDG0Data(ciceName);
    uice = kernel.getDG0Data(uName);
    vice = kernel.getDG0Data(vName);
}

}
