/*!
 * @file IStructureModule.cpp
 *
 * @date Feb 18, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/IStructureModule.hpp"

#include "include/DevGrid.hpp"
#include "include/ParametricGrid.hpp"
#include "include/RectangularGrid.hpp"

namespace Module {

const std::string RECT_GRID = "RectangularGrid";
const std::string PARAMETRICGRID = "ParametricGrid";

template <>
Module<Nextsim::IStructure>::map Module<Nextsim::IStructure>::functionMap = {
    { DEV_GRID, newImpl<Nextsim::IStructure, Nextsim::DevGrid> },
    { RECT_GRID, newImpl<Nextsim::IStructure, Nextsim::RectangularGrid> },
    { PARAMETRICGRID, newImpl<Nextsim::IStructure, Nextsim::ParametricGrid> },
};
template <>
Module<Nextsim::IStructure>::fn Module<Nextsim::IStructure>::spf = functionMap.at(DEV_GRID);
template <>
std::unique_ptr<Nextsim::IStructure> Module<Nextsim::IStructure>::staticInstance
    = std::move(newImpl<Nextsim::IStructure, Nextsim::DevGrid>());

template <> std::string Module<Nextsim::IStructure>::moduleName() { return "IStructure"; }

template <> Nextsim::IStructure& getImplementation<Nextsim::IStructure>()
{
    return getImplTemplate<Nextsim::IStructure, IStructureModule>();
}

template <> void setImplementation<Nextsim::IStructure>(const std::string& implName)
{
    setImplTemplate<IStructureModule>(implName);
}

template <> std::unique_ptr<Nextsim::IStructure> getInstance()
{
    return getInstTemplate<Nextsim::IStructure, IStructureModule>();
}
} /* namespace Module */
