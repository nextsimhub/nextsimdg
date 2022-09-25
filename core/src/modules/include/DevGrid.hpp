/*!
 * @file DevGrid.hpp
 *
 * @date Dec 20, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef DEVGRID_HPP
#define DEVGRID_HPP

#include "include/IStructure.hpp"

#include "include/IDevGridIO.hpp"
#include "include/ModelState.hpp"

#include <map>

namespace Nextsim {

class DevGridIO;

//! A class to hold a grid of ElementData instances in a fixed sized square grid.
class DevGrid : public IStructure {
public:
    DevGrid()
        : pio(nullptr)
    {
    }

#ifdef USE_MPI
    DevGrid(MPI_Comm comm)
        : IStructure(comm)
        , pio(nullptr)
    {
    }
#endif // USE_MPI

    //! Destructor. The lifetime of pio should be the lifetime of the instance.
    virtual ~DevGrid()
    {
        if (pio) {
            delete pio;
        }
    }

    const static int nx;
#ifdef USE_MPI
    static int idx_x; // Global coordinate in "x" dimension of upper left corner
    static int idx_y; // Global coordinate in "y" dimension of upper left corner
    static int local_nx; // Local extent in "x" dimension
    static int local_ny; // Local extent in "y" dimension
#endif // USE_MPI
    const static std::string structureName;

    // Read/write override functions
#ifdef USE_MPI
    ModelState getModelState(
        const std::string& restartFilePath, const std::string& partitionFilePath) override
    {
        return pio ? pio->getModelState(restartFilePath, partitionFilePath) : ModelState();
    }
#else
    ModelState getModelState(const std::string& filePath) override
    {
        return pio ? pio->getModelState(filePath) : ModelState();
    }
#endif // USE_MPI

    void dumpModelState(const ModelState& state, const ModelMetadata& metadata,
        const std::string& filePath, bool isRestart = false) const override
    {
        if (pio)
            pio->dumpModelState(state, metadata, filePath, isRestart);
    }
    const std::string& structureType() const override { return structureName; };

    int nIceLayers() const override { return 1; };

    //! Sets the pointer to the class that will perform the IO. Should be an instance of DevGridIO
    void setIO(IDevGridIO* p) { pio = p; }

    const static std::string xDimName;
    const static std::string yDimName;
    const static std::string nIceLayersName;

private:
    IDevGridIO* pio;

    friend DevGridIO;
};

} /* namespace Nextsim */

#endif /* DEVGRID_HPP */
