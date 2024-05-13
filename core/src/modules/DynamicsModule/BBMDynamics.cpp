/*!
 * @file BBMDynamics.cpp
 *
 * @date Jan 5, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/BBMDynamics.hpp"

#include "include/gridNames.hpp"

namespace Nextsim {

static const std::vector<std::string> namedFields = { hiceName, ciceName, uName, vName };
BBMDynamics::BBMDynamics()
    : IDynamics(true)
    , kernel(params)
{
    getStore().registerArray(Protected::ICE_U, &uice, RO);
    getStore().registerArray(Protected::ICE_V, &vice, RO);
}

void BBMDynamics::setData(const ModelState::DataMap& ms)
{
    IDynamics::setData(ms);

    bool isSpherical = checkSpherical(ms);

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

    // ModelArray& coords = ms.at(coordsName);
    // TODO: Some encoding of the periodic edge boundary conditions
    kernel.initialise(coords, isSpherical, ms.at(maskName));

    uice = ms.at(uName);
    vice = ms.at(vName);

    // Set the data in the kernel arrays.
    for (const auto& fieldName : namedFields) {
        kernel.setData(fieldName, ms.at(fieldName));
    }
}

void BBMDynamics::update(const TimestepTime& tst)
{
    std::cout << tst.start << std::endl;

    // Fill the updated damage array with the initial value
    damage = damage0.data();

    // set the updated ice thickness, concentration and damage
    kernel.setData(hiceName, hice.data());
    kernel.setData(ciceName, cice.data());
    kernel.setData(damageName, damage);

    // set the forcing velocities
    kernel.setData(uWindName, uwind.data());
    kernel.setData(vWindName, vwind.data());
    kernel.setData(uOceanName, uocean.data());
    kernel.setData(vOceanName, vocean.data());

    /*
     * Ice velocity components are stored in the dynamics, and not changed by the model outside the
     * dynamics kernel. Hence they are not set at this point.
     */

    kernel.update(tst);

    hice.data() = kernel.getDG0Data(hiceName);
    cice.data() = kernel.getDG0Data(ciceName);
    damage = kernel.getDG0Data(damageName);

    uice = kernel.getDG0Data(uName);
    vice = kernel.getDG0Data(vName);
}

} /* namespace Nextsim */