/*!
 * @file RectGridIO.hpp
 *
 * @date Feb 8, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_RECTGRIDIO_HPP
#define CORE_SRC_INCLUDE_RECTGRIDIO_HPP

#include "include/RectangularGrid.hpp"

namespace Nextsim {

class RectGridIO : public RectangularGrid::IRectGridIO {
public:
    RectGridIO(RectangularGrid& grid)
        : IRectGridIO(grid)
    {
    }
    virtual ~RectGridIO() = default;

    typedef RectangularGrid::GridDimensions GridDimensions;

    void init(
        std::vector<ElementData>& dg, const std::string& filePath, GridDimensions& dims) override;

    ModelState getModelState(const std::string& filePath) override;

    void dump(const std::vector<ElementData>& dg, const std::string& filePath,
        const GridDimensions& dims) const override;

    void dumpModelState(const ModelState& state, const std::string& filePath) const override;

private:
    RectGridIO() = default;
};

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_RECTGRIDIO_HPP */
