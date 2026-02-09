/*!
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef FREEDRIFTDYNAMICS_HPP
#define FREEDRIFTDYNAMICS_HPP

#include "include/FreeDriftDynamicsKernel.hpp"
#include "include/IDynamics.hpp"
#include "include/dgVectorHolder.hpp"

#include "include/ModelArray.hpp"
#include "include/ModelComponent.hpp"
#include "include/gridNames.hpp"

#ifndef DGCOMP
#define DGCOMP 3 // Define to prevent errors from static analysis tools
#error "Number of DG components (DGCOMP) not defined" // But throw an error anyway
#endif

namespace Nextsim {
static const std::vector<std::string> namedFields = { uName, vName };

// TODO: This class should be configurable and configure the densities, drag coefficients, Coriolis
// parameter, and turning angle.
class FreeDriftDynamics : public IDynamics {
public:
    FreeDriftDynamics()
        : IDynamics()
        , kernel(params)
    {
    }

    std::string getName() const override { return "FreeDriftDynamics"; }
    void update(const TimestepTime& tst) override
    {
        std::cout << tst.start << std::endl;

        // set the forcing velocities
        kernel.setData(uOceanName, uoceanAccessor.getHostRO());
        kernel.setData(vOceanName, voceanAccessor.getHostRO());

        kernel.update(tst);

        uiceAccessor.getHostRW() = kernel.getDG0Data(uName);
        viceAccessor.getHostRW() = kernel.getDG0Data(vName);
    }

    void setData(const ModelState::DataMap& ms) override
    {
        IDynamics::setData(ms);

        // Set the DG field data. Needs to be done before initialise() because the Kokkos kernel
        // requires the actual vectors to init its views.
        kernel.setDGArray(hiceName, hiceDGAccessor.getHostRW());
        kernel.setDGArray(ciceName, ciceDGAccessor.getHostRW());
        kernel.setDGArray(hsnowName, hsnowDGAccessor.getHostRW());

        bool isSpherical = checkSpherical(ms);

        ModelArray coords = ms.at(coordsName);
        if (isSpherical) {
            coords *= PhysicalConstants::deg2rad;
        }
        // TODO: Some encoding of the periodic edge boundary conditions
        kernel.initialise(coords, false, ms.at(maskName));

        // Set the data in the kernel arrays.
        for (const auto& fieldName : namedFields) {
            kernel.setData(fieldName, ms.at(fieldName));
        }
    }

    void advectField(double timestep, ModelArray& field,
        double lowerLimit = -std::numeric_limits<double>::infinity(),
        double upperLimit = std::numeric_limits<double>::infinity()) override
    {
        kernel.advectField(timestep, field, lowerLimit, upperLimit);
    }

    void prepareAdvection() override { kernel.prepareAdvection(); }

private:
    FreeDriftDynamicsKernel<DGCOMP> kernel;
    DynamicsParameters params;
};
}

#endif /* FREEDRIFTDYNAMICS_HPP */
