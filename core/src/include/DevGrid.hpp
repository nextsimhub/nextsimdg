/*!
 * @file DevGrid.hpp
 *
 * @date Dec 20, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_DEVGRID_HPP
#define CORE_SRC_INCLUDE_DEVGRID_HPP

#include "IStructure.hpp"

namespace Nextsim {

class DevGrid : public IStructure {
public:
    DevGrid();
    virtual ~DevGrid();

    void init(netCDF::NcGroup& grp) override;
    void dumpMeta(netCDF::NcGroup& metaGroup) const override;
    void dumpData(netCDF::NcGroup& dataGroup) const override;

private:
    const static std::string ourStructureName;
};

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_DEVGRID_HPP */
