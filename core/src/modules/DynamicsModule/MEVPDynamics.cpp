/*!
 * @file MEVPDynamics.cpp
 *
 * @date 18 Jul 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#include "include/MEVPDynamics.hpp"

#include "include/gridNames.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace Nextsim {

// Degrees to radians as a hex float
static const double radians = 0x1.1df46a2529d39p-6;

void MEVPDynamics::configure()
{
    Module::Module<Nextsim::IDamageHealing>::setImplementation("Nextsim::NoHealing");
}

static const std::vector<std::string> namedFields = { hiceName, ciceName, hsnowName, uName, vName };
MEVPDynamics::MEVPDynamics()
    : IDynamics()
    , kernel(params)
{
    getStore().registerArray(Protected::ICE_U, &uice, RO);
    getStore().registerArray(Protected::ICE_V, &vice, RO);
}

void MEVPDynamics::setData(const ModelState::DataMap& ms)
{
    IDynamics::setData(ms);

    bool isSpherical = checkSpherical(ms);

    ModelArray coords = ms.at(coordsName);
    if (isSpherical) {
        coords *= radians;
    }
    // TODO: Some encoding of the periodic edge boundary conditions
    kernel.initialise(coords, isSpherical, ms.at(maskName));

    uice = ms.at(uName);
    vice = ms.at(vName);

    // Set the data in the kernel arrays.
    for (const auto& fieldName : namedFields) {
        kernel.setData(fieldName, ms.at(fieldName));
    }
}

void MEVPDynamics::update(const TimestepTime& tst)
{
    std::cout << tst.start << std::endl;

    // set the updated ice thickness and concentration
    kernel.setData(hiceName, hice.data());
    kernel.setData(ciceName, cice.data());
    kernel.setData(hsnowName, hsnow.data());

    // set the forcing velocities
    kernel.setData(uWindName, uwind.data());
    kernel.setData(vWindName, vwind.data());
    kernel.setData(uOceanName, uocean.data());
    kernel.setData(vOceanName, vocean.data());

    // kernel.setData(uName, uice);
    // kernel.setData(vName, vice);

    kernel.update(tst);

    hice.data() = kernel.getDG0Data(hiceName);
    cice.data() = kernel.getDG0Data(ciceName);
    hsnow.data() = kernel.getDG0Data(hsnowName);

    uice = kernel.getDG0Data(uName);
    vice = kernel.getDG0Data(vName);
}

// All data for prognostic output
ModelState MEVPDynamics::getState() const
{
    // Get the velocities from IDynamics
    ModelState state(IDynamics::getState());

    // Kernel prognostic fields
    state.merge({
        { hiceName, kernel.getDGData(hiceName) },
        { ciceName, kernel.getDGData(ciceName) },
        { hsnowName, kernel.getDGData(hsnowName) },
    });

    return state;
}

ModelState MEVPDynamics::getStateRecursive(const OutputSpec& os) const
{
    // Base class state
    ModelState state(IDynamics::getStateRecursive(os));

    if (os.allComponents()) {
        state.merge({
            { hiceName, kernel.getDGData(hiceName) },
            { ciceName, kernel.getDGData(ciceName) },
            { hsnowName, kernel.getDGData(hsnowName) },
        });
    }
    return state;
}

}
