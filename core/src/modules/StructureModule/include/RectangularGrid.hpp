/*!
 *
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
    };

    RectangularGrid()
        : pio(nullptr)
        , nx(0)
        , ny(0)
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
    ModelState getModelState(const std::string& filePath) override
    {
        return pio ? pio->getModelState(filePath) : ModelState();
    }

    void dumpModelState(
        const ModelState& state, const std::string& filePath, bool isRestart = false) const override
    {
        if (pio)
            pio->dumpModelState(state, filePath, isRestart);
    }
    const std::string& structureType() const override { return structureName; };

    void setDimensions(const GridDimensions& dims)
    {
        nx = dims.nx;
        ny = dims.ny;
    }

    class IRectGridIO {
    public:
        IRectGridIO(RectangularGrid& grid)
            : grid(&grid)
        {
        }
        virtual ~IRectGridIO() = default;

        virtual ModelState getModelState(const std::string& filePath) = 0;

        /*!
         * @brief Dumps the given ModelState to the given file path.
         *
         * @param state The ModelState data
         * @param filePath The path to attempt to write the data to.
         */
        virtual void dumpModelState(
            const ModelState& state, const std::string& filePath, bool isRestart) const
            = 0;

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

    const static std::string xDimName;
    const static std::string yDimName;

    IRectGridIO* pio;

    friend RectGridIO;
};

} /* namespace Nextsim */

#endif /* RECTANGULARGRID_HPP */
