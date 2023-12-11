/*!
 * @file RectangularGrid.hpp
 *
 * @date Feb 7, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Kacper Kornet <kk562@cam.ac.uk>
 */

#ifndef RECTANGULARGRID_HPP
#define RECTANGULARGRID_HPP

#include "include/IStructure.hpp"

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

    RectangularGrid(const GridDimensions& dims)
        : pio(nullptr)
    {
        setDimensions(dims);
    }

    virtual ~RectangularGrid()
    {
        if (pio) {
            delete pio;
        }
    }

    // Read/write override functions
#ifdef USE_MPI
    ModelState getModelState(const std::string& filePath, const std::string& partitionFile,
        ModelMetadata& metadata) override
    {
        return pio ? pio->getModelState(filePath, partitionFile, metadata) : ModelState();
    }
#else
    ModelState getModelState(const std::string& filePath) override
    {
        return pio ? pio->getModelState(filePath) : ModelState();
    }
#endif

    void dumpModelState(
        const ModelState& state, const ModelMetadata& metadata, const std::string& filePath, bool isRestart = false) const override
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
        virtual ModelState getModelState(
            const std::string& filePath, const std::string& partitionFile, ModelMetadata& metadata)
            = 0;
#else
        virtual ModelState getModelState(const std::string& filePath) = 0;
#endif

        /*!
         * @brief Dumps the given ModelState to the given file path.
         *
         * @param state The ModelState data
         * @param filePath The path to attempt to write the data to.
         */
        virtual void dumpModelState(
            const ModelState& state, const ModelMetadata& metadata, const std::string& filePath, bool isRestart) const = 0;

    protected:
        IRectGridIO() = default;

    private:
        RectangularGrid* grid;
    };

    //! Sets the pointer to the class that will perform the IO. Should be an instance of IRectGridIO
    void setIO(IRectGridIO* p) { pio = p; }

    const static std::string structureName;

private:
    int nx;
    int ny;
    int nz; // Number of ice layers

    const static std::string xDimName;
    const static std::string yDimName;
    const static std::string nIceLayersName;

    IRectGridIO* pio;

    friend RectGridIO;
};

} /* namespace Nextsim */

#endif /* RECTANGULARGRID_HPP */
