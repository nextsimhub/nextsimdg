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

namespace Nextsim {
static const std::vector<std::string> namedFields = { hiceName, ciceName, uName, vName };
class FreeDriftDynamics : public IDynamics {
public:
    FreeDriftDynamics()
        : IDynamics()
        , kernel()
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
        IDynamics::setData(ms);

        bool isSpherical = checkSpherical(ms);

        // TODO: Remove this when spherical coordinates are fully implemented
        ModelArray coords;
        if (isSpherical) {
            std::cout << "Spherical coordinates are not yet implemented. Reverting to Cartesian." << std::endl;
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
            coords = fake25kmCoords;
            // End of code to be removed
        } else {
            coords = ms.at(coordsName);
        }
        kernel.initialise(coords, false, ms.at(maskName));

        // Set the data in the kernel arrays.
        for (const auto& fieldName : namedFields) {
            kernel.setData(fieldName, ms.at(fieldName));
        }
    }

private:
    // TODO: How to get the template parameters here?
    FreeDriftDynamicsKernel<6> kernel;
};
}

#endif /* FREEDRIFTDYNAMICS_HPP */
