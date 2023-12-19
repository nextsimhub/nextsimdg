/*!
 * @file Dynamics.cpp
 *
 * @date 7 Sep 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#include "include/Dynamics.hpp"

#include "include/gridNames.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace Nextsim {

static const std::vector<std::string> namedFields = { hiceName, ciceName, uName, vName };
Dynamics::Dynamics()
        : IDynamics()
{
    getStore().registerArray(Protected::ICE_U, &uice, RO);
    getStore().registerArray(Protected::ICE_V, &vice, RO);
}

void Dynamics::setData(const ModelState::DataMap& ms)
{
    IDynamics::setData(ms);

    bool isSpherical;
    // Decide between Cartesian (x & y) and spherical (longitude & latitude)
    if (ms.count(longitudeName) > 0 && ms.count(latitudeName) > 0) {
        isSpherical = true;
    } else if (ms.count(xName) > 0 && ms.count(yName) > 0) {
        isSpherical = false;
    } else {
        // Throw a runtime_error exception which can either be handled or not
        throw std::runtime_error(
                "Input data must contain either Cartesian (" + xName + ", " + yName
                        + ") or spherical (" + longitudeName + ", " + latitudeName
                        + ") coordinates.");
    }

    // TODO: Remove this when spherical coordinates are fully implemented
    if (isSpherical) std::cout << "Spherical coordinates are not yet implemented. Reverting to Cartesian." << std::endl;
    isSpherical = false;
    ModelArray fake25kmCoords(ModelArray::Type::VERTEX);
    double d = 25000; // 25 km in metres
    // Fill the fake coordinate array
    for (size_t j = 0; j < ModelArray::size(ModelArray::Dimension::YVERTEX); ++j) {
        for (size_t i = 0; i < ModelArray::size(ModelArray::Dimension::XVERTEX); ++i) {
            fake25kmCoords.components({i, j})[0] = d * i;
            fake25kmCoords.components({i, j})[1] = d * j;
        }
    }
    ModelArray& coords = fake25kmCoords;
    // End of code to be removed

    std::cerr << "dynamics kernel initialization" << std::endl;
    std::cerr << "mask.size() = " << ms.at(maskName).trueSize() << std::endl;
    // ModelArray& coords = ms.at(coordsName);
    // TODO: Some encoding of the periodic edge boundary conditions
    kernel.initialisation(coords, isSpherical, ms.at(maskName));

    uice = ms.at(uName);
    vice = ms.at(vName);

    // Set the data in the kernel arrays.
    for (const auto &fieldName : namedFields) {
        kernel.setData(fieldName, ms.at(fieldName));
    }
}

void Dynamics::update(const TimestepTime& tst)
{
    std::cout << tst.start << std::endl;

    // set the updated ice thickness and concentration
    kernel.setData(hiceName, hice.data());
    kernel.setData(ciceName, cice.data());

    // set the forcing velocities
    kernel.setData(uWindName, uwind.data());
    kernel.setData(vWindName, vwind.data());
    kernel.setData(uOceanName, uocean.data());
    kernel.setData(vOceanName, vocean.data());

    //kernel.setData(uName, uice);
    //kernel.setData(vName, vice);

    kernel.update(tst);

    hice.data() = kernel.getDG0Data(hiceName);
    cice.data() = kernel.getDG0Data(ciceName);

    uice = kernel.getDG0Data(uName);
    vice = kernel.getDG0Data(vName);
}

}
