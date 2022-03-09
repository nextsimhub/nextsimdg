/*!
 * @file RectangularGrid.hpp
 *
 * @date Feb 7, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_MODULES_INCLUDE_RECTANGULARGRID_HPP_
#define CORE_SRC_MODULES_INCLUDE_RECTANGULARGRID_HPP_

#include "IStructure.hpp"

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
        data.resize(nx * ny);
    }

    virtual ~RectangularGrid()
    {
        if (pio) {
            delete pio;
        }
    }

    // Read/write override functions
    void init(const std::string& filePath) override;

    ModelState getModelState(const std::string& filePath) override
    {
        return pio ? pio->getModelState(filePath) : ModelState();
    }

    void dump(const std::string& filePath) const override;

    void dumpModelState(const ModelState& state, const std::string& filePath) const override
    {
        if (pio)
            pio->dumpModelState(state, filePath);
    }
    std::string structureType() const override { return structureName; };

    int nIceLayers() const override { return nz; };

    // Cursor manipulation override functions
    int resetCursor() override;
    bool validCursor() const override;
    ElementData& cursorData() override;
    const ElementData& cursorData() const override;
    void incrCursor() override;

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

        virtual void init(
            std::vector<ElementData>& dg, const std::string& filePath, GridDimensions& dims)
            = 0;
        /*!
         * @brief Writes data from the vector of data elements into the file location.
         *
         * @param dg The vector of ElementData instances containing the data.
         * @param filePath The location of the NetCDF restart file to be written.
         */
        virtual void dump(const std::vector<ElementData>& dg, const std::string& filePath,
            const GridDimensions& dims) const = 0;

        virtual ModelState getModelState(const std::string& filePath) = 0;

        /*!
         * @brief Dumps the given ModelState to the given file path.
         *
         * @param state The ModelState data
         * @param filePath The path to attempt to write the data to.
         */
        virtual void dumpModelState(const ModelState& state, const std::string& filePath) const = 0;
    protected:
        IRectGridIO() = default;

    private:
        RectangularGrid* grid;
    };

    //! Sets the pointer to the class that will perform the IO. Should be an instance of DevGridIO
    void setIO(IRectGridIO* p) { pio = p; }

    const static std::string structureName;

private:
    int nx;
    int ny;
    int nz; // Number of ice layers

    const static std::string xDimName;
    const static std::string yDimName;
    const static std::string nIceLayersName;

    std::vector<ElementData> data;

    std::vector<ElementData>::iterator iCursor;

    IRectGridIO* pio;

    friend RectGridIO;
};

} /* namespace Nextsim */

#endif /* CORE_SRC_MODULES_INCLUDE_RECTANGULARGRID_HPP_ */
