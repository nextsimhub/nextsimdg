/*!
 * @file DevGridIO.hpp
 *
 * @date Jan 14, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef DEVGRIDIO_HPP
#define DEVGRIDIO_HPP

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

    ModelState getModelState(const std::string& filePath) const override;
    void dumpModelState(const ModelState& state, const std::string& filePath) const override;

private:
    DevGrid* grid;
};

} /* namespace Nextsim */

#endif /* DEVGRIDIO_HPP */
