/*!
 * @file BBMDynamics.cpp
 *
 * @date 07 Oct 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#include "include/BBMDynamics.hpp"

#include "include/gridNames.hpp"

namespace Nextsim {

static const std::vector<std::string> namedFields = { hiceName, ciceName, uName, vName };
static const std::map<std::string, std::pair<ModelArray::Type, double>> defaultFields = {
    { damageName, { ModelArray::Type::H, 1.0 } },
};
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

    ModelArray coords = ms.at(coordsName);
    if (isSpherical) {
        coords *= radians;
    }
    // TODO: Some encoding of the periodic edge boundary conditions
    kernel.initialise(coords, isSpherical, ms.at(maskName));

    uice = ms.at(uName);
    vice = ms.at(vName);

    // Set the data in the kernel arrays.
    // Required data
    for (const auto& fieldName : namedFields) {
        kernel.setData(fieldName, ms.at(fieldName));
    }
    // Data that can have a default value
    for (const auto entry : defaultFields) {
        // Directly add data that is supplied
        const std::string& fieldName = entry.first;
        if (ms.count(fieldName) > 0) {
            kernel.setData(fieldName, ms.at(fieldName));
        } else {
            // Fill data that is not supplied, masking if the mask is available
            ModelArray data(entry.second.first);
            data.resize();
            // Fill the default value
            data = entry.second.second;
            // Mask the default data
            kernel.setData(fieldName, mask(data));
        }
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
    kernel.setData(sshName, ssh.data());

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

    taux = kernel.getDG0Data(uIOStressName);
    tauy = kernel.getDG0Data(vIOStressName);
}

// All data for prognostic output
ModelState BBMDynamics::getState() const
{
    // Get the velocities from IDynamics
    ModelState state(IDynamics::getState());

    // Kernel prognostic fields
    state.merge({
        { hiceName, kernel.getDGData(hiceName) },
        { ciceName, kernel.getDGData(ciceName) },
        { damageName, kernel.getDGData(damageName) },
    });

    return state;
}

ModelState BBMDynamics::getStateRecursive(const OutputSpec& os) const
{
    // Base class state
    ModelState state(IDynamics::getStateRecursive(os));

    if (os.allComponents()) {
        state.merge({
            { hiceName, kernel.getDGData(hiceName) },
            { ciceName, kernel.getDGData(ciceName) },
            { damageName, kernel.getDGData(damageName) },
        });
    }
    return state;
}

} /* namespace Nextsim */
