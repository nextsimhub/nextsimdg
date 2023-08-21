/*!
 * @file Dynamics.cpp
 *
 * @date 27 Mar 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#include "include/Dynamics.hpp"

#include "include/gridNames.hpp"

#include <string>
#include <vector>


namespace Nextsim {

static const std::vector<std::string> namedFields = { hiceName, ciceName, uName, vName };
Dynamics::Dynamics()
    : IDynamics()
{
    registerProtectedArray(ProtectedArray::ICE_U, &uice);
    registerProtectedArray(ProtectedArray::ICE_V, &vice);
}

void Dynamics::setData(const ModelState::DataMap& ms)
{
    IDynamics::setData(ms);

    kernel.initialisation();

    uice = ms.at(uName);
    vice = ms.at(vName);

    // Set the data in the kernel arrays.
    for (const auto& fieldName : namedFields) {
        kernel.setData(fieldName, ms.at(fieldName));
    }
}



void Dynamics::update(const TimestepTime& tst)
{   std::cout << tst.start << std::endl;

    kernel.setData(hiceName, hice.data());
    kernel.setData(ciceName, cice.data());

    //kernel.setData(uName, uice);
    //kernel.setData(vName, vice);
    
    kernel.update(tst);

    hice.data() = kernel.getDG0Data(hiceName);
    cice.data() = kernel.getDG0Data(ciceName);
    
    //uice = kernel.getDG0Data(uName);
    //vice = kernel.getDG0Data(vName);
}

}
