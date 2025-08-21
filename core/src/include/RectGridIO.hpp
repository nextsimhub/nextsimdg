/*!
 *
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Kacper Kornet <kk562@cam.ac.uk>
*/

#ifndef RECTGRIDIO_HPP
#define RECTGRIDIO_HPP

#include "StructureModule/include/RectangularGrid.hpp"

namespace Nextsim {

class RectGridIO : public RectangularGrid::IRectGridIO {
public:
    RectGridIO(RectangularGrid& grid)
        : IRectGridIO(grid)
    {
    }
    virtual ~RectGridIO() = default;

    typedef RectangularGrid::GridDimensions GridDimensions;

    ModelState getModelState(const std::string& filePath) override;

    void dumpModelState(
        const ModelState& state, const std::string& filePath, bool isRestart) const override;

private:
    RectGridIO() = default;
};

} /* namespace Nextsim */

#endif /* RECTGRIDIO_HPP */
