/*!
 * @file RectangularGrid.hpp
 *
 * @date Feb 7, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef RECTANGULARGRID_HPP
#define RECTANGULARGRID_HPP

#include "IStructure.hpp"

#include "include/ModelMetadata.hpp"
#include "include/ModelState.hpp"

namespace Nextsim {

class RectGridIO;

class RectangularGrid : public IStructure {
public:
    struct GridDimensions {
        int nx;
        int ny;
        int nz;
    };

    RectangularGrid()
        : pio(nullptr)
        , nx(0)
        , ny(0)
        , nz(0)
    {
    }

#ifdef USE_MPI
    RectangularGrid(MPI_Comm comm)
        : IStructure(comm)
        , pio(nullptr)
        , nx(0)
        , ny(0)
        , nz(0)
    {
    }
#endif // USE_MPI

    RectangularGrid(const GridDimensions& dims)
        : pio(nullptr)
    {
        setDimensions(dims);
    }

#ifdef USE_MPI
    RectangularGrid(MPI_Comm comm, const GridDimensions& dims)
        : IStructure(comm)
        , pio(nullptr)
    {
        setDimensions(dims);
    }
#endif // USE_MPI

    virtual ~RectangularGrid()
    {
        if (pio) {
            delete pio;
        }
    }

    // Read/write override functions
#ifdef USE_MPI
    ModelState getModelState(
        const std::string& modelFilePath, const std::string& partitionFilePath) override
    {
        return pio ? pio->getModelState(modelFilePath, partitionFilePath) : ModelState();
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

    int nIceLayers() const override { return nz; };

    void setDimensions(const GridDimensions& dims)
    {
        nx = dims.nx;
        ny = dims.ny;
        nz = dims.nz;
    }

    class IRectGridIO {
    public:
        IRectGridIO(RectangularGrid& grid)
            : grid(&grid)
        {
        }
        virtual ~IRectGridIO() = default;

#ifdef USE_MPI
        /*!
         * @brief Generates the ModelState based on the data in the given file for an MPI process.
         *
         * @param modelFilePath the name of the model file to be read.
         * @param partitionFilePath the name of the partitioining file to be read.
         * @param rank the MPI rank of the process.
         * @param comm the MPI communicator.
         */
        virtual ModelState getModelState(
            const std::string& modelFilePath, const std::string& partitionFilePath) const
            = 0;
#else
        virtual ModelState getModelState(const std::string& filePath) = 0;
#endif // USE_MPI

        /*!
         * @brief Dumps the given ModelState to the given file path.
         *
         * @param state The ModelState data
         * @param filePath The path to attempt to write the data to.
         */
        virtual void dumpModelState(const ModelState& state, const ModelMetadata& metadata,
            const std::string& filePath, bool isRestart) const
            = 0;

    protected:
        IRectGridIO() = default;
        RectangularGrid* grid;
    };

    //! Sets the pointer to the class that will perform the IO. Should be an instance of DevGridIO
    void setIO(IRectGridIO* p) { pio = p; }

    const static std::string structureName;

private:
    int nx;
    int ny;
    int nz; // Number of ice layers
#ifdef USE_MPI
    int global_nx; // Global extent in "x" dimension
    int global_ny; // Global extent in "y" dimension
    int idx_x; // Global coordinate in "x" dimension of upper left corner
    int idx_y; // Global coordinate in "y" dimension of upper left corner
#endif // USE_MPI

    const static std::string xDimName;
    const static std::string yDimName;
    const static std::string nIceLayersName;

    IRectGridIO* pio;

    friend RectGridIO;
};

} /* namespace Nextsim */

#endif /* RECTANGULARGRID_HPP */
