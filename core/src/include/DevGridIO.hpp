/*!
 * @file DevGridIO.hpp
 *
 * @date Jan 14, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_DEVGRIDIO_HPP
#define CORE_SRC_INCLUDE_DEVGRIDIO_HPP

#include "include/ElementData.hpp"
#include "include/IDevGridIO.hpp"

#include <vector>

namespace Nextsim {

class DevGrid;

class DevGridIO : public IDevGridIO {
public:
    DevGridIO(DevGrid& grid)
        : IDevGridIO(grid)
    {
    }
    virtual ~DevGridIO() = default;

    void init(std::vector<ElementData>& data, const std::string& filePath) const override;
    void dump(const std::vector<ElementData>& data, const std::string& filePath) const override;
    ModelState getModelState(const std::string& filePath) const override;
    void dumpModelState(const ModelState& state, const std::string& filePath) const override;

private:
    DevGrid* grid;
};

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_DEVGRIDIO_HPP */
