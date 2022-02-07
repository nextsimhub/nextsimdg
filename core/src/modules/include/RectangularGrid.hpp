/*!
 * @file RectangularGrid.hpp
 *
 * @date Feb 7, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_MODULES_INCLUDE_RECTANGULARGRID_HPP_
#define CORE_SRC_MODULES_INCLUDE_RECTANGULARGRID_HPP_

#include "IStructure.hpp"

namespace Nextsim {

class RectangularGrid : public IStructure {
public:
    RectangularGrid()
        : pio (nullptr)
    {
    }

    virtual ~RectangularGrid()
    {
        if (pio) {
            delete pio;
        }
    }

    // Read/write override functions
    void init(const std::string& filePath) override;

    void dump(const std::string& filePath) const override;


    std::string structureType() const override { return structureName; };

    int nIceLayers() const override { return nz; };

    // Cursor manipulation override functions
    int resetCursor() override;
    bool validCursor() const override;
    ElementData& cursorData() override;
    const ElementData& cursorData() const override;
    void incrCursor() override;

    class IRectGridIO {
    public:
        IRectGridIO(RectangularGrid& grid)
            : grid(&grid)
        {
        }
        virtual ~IRectGridIO() = default;

        virtual void init(std::vector<ElementData>& dg, const std::string& filePath) const = 0;
        /*!
         * @brief Writes data from the vector of data elements into the file location.
         *
         * @param dg The vector of ElementData instances containing the data.
         * @param filePath The location of the NetCDF restart file to be written.
         */
        virtual void dump(const std::vector<ElementData>& dg, const std::string& fielPath) const = 0;

    private:
        RectangularGrid* grid;
};

private:
    const static std::string structureName;
    // Only one size of RectangularGrid at a time
    static int nx;
    static int ny;
    static int nz; // Number of ice layers

    const static std::string xDimName;
    const static std::string yDimName;
    const static std::string nIceLayersName;

    std::vector<ElementData> data;

    std::vector<ElementData>::iterator iCursor;

    IRectGridIO* pio;
};

} /* namespace Nextsim */

#endif /* CORE_SRC_MODULES_INCLUDE_RECTANGULARGRID_HPP_ */
