/*!
 * @file RectGridIO.hpp
 *
 * @date Feb 8, 2022
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

#ifdef USE_MPI
    ModelState getModelState(const std::string& filePath, const std::string& partitionFile,
        ModelMetadata& metadata) override;
#else
    ModelState getModelState(const std::string& filePath) override;
#endif

    void dumpModelState(const ModelState& state, const ModelMetadata& metadata,
        const std::string& filePath, bool isRestart) const override;

private:
    RectGridIO() = default;

    void readPartitionData(const std::string& partitionFile, ModelMetadata& metadata);
};

} /* namespace Nextsim */

#endif /* RECTGRIDIO_HPP */
