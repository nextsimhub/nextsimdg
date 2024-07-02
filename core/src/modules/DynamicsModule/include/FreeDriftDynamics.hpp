/*!
 * @file FreeDriftDynamics.hpp
 *
 * @date 27 Mar 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef FREEDRIFTDYNAMICS_HPP
#define FREEDRIFTDYNAMICS_HPP

#include "include/FreeDriftDynamicsKernel.hpp"
#include "include/IDynamics.hpp"

#include "include/ModelArray.hpp"
#include "include/ModelComponent.hpp"
#include "include/gridNames.hpp"

#ifndef DGCOMP
#define DGCOMP 3 // define to make red lines go away in the IDE
#error "Number of DG components (DGCOMP) not defined" // But throw an error anyway
#endif

namespace Nextsim {
static const std::vector<std::string> namedFields = { hiceName, ciceName, uName, vName };


class FreeDriftDynamics : public IDynamics {
public:
    FreeDriftDynamics()
        : IDynamics()
        , kernel(params)
    {
        getStore().registerArray(Protected::ICE_U, &uice, RO);
        getStore().registerArray(Protected::ICE_V, &vice, RO);
    }

    std::string getName() const override { return "FreeDriftDynamics"; }
    void update(const TimestepTime& tst) override
    {
        std::cout << tst.start << std::endl;

        // set the updated ice thickness and concentration
        kernel.setData(hiceName, hice.data());
        kernel.setData(ciceName, cice.data());

        // set the forcing velocities
        kernel.setData(uOceanName, uocean.data());
        kernel.setData(vOceanName, vocean.data());

        kernel.update(tst);

        hice.data() = kernel.getDG0Data(hiceName);
        cice.data() = kernel.getDG0Data(ciceName);

        uice = kernel.getDG0Data(uName);
        vice = kernel.getDG0Data(vName);
    }

    void setData(const ModelState::DataMap& ms) override
    {
        // Degrees to radians as a hex float
        static const double radians = 0x1.1df46a2529d39p-6;

        IDynamics::setData(ms);

        bool isSpherical = checkSpherical(ms);

        ModelArray coords = ms.at(coordsName);
        if (isSpherical) {
            coords *= radians;
        }
        // TODO: Some encoding of the periodic edge boundary conditions
        kernel.initialise(coords, false, ms.at(maskName));

        // Set the data in the kernel arrays.
        for (const auto& fieldName : namedFields) {
            kernel.setData(fieldName, ms.at(fieldName));
        }
    }

private:
    FreeDriftDynamicsKernel<DGCOMP> kernel;
    DynamicsParameters params;
};
}

#endif /* FREEDRIFTDYNAMICS_HPP */
