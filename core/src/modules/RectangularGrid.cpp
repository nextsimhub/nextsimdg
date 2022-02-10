/*!
 * @file RectangularGrid.cpp
 *
 * @date Feb 7, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/RectangularGrid.hpp"

namespace Nextsim {

const std::string RectangularGrid::structureName = "simple_rectangular";
const std::string RectangularGrid::xDimName = "x";
const std::string RectangularGrid::yDimName = "y";
const std::string RectangularGrid::nIceLayersName = "nLayers";

void RectangularGrid::init(const std::string& filePath)
{
    GridDimensions dims;
    data = std::vector<ElementData>();
    if (pio && !filePath.empty()) {
        pio->init(data, filePath, dims);
    }
    setDimensions(dims);
}

void RectangularGrid::dump(const std::string& filePath) const
{
    GridDimensions dims = { nx, ny, nz };
    if (pio && !filePath.empty()) {
        pio->dump(data, filePath, dims);
    }
}

// Cursor manipulation override functions
int RectangularGrid::resetCursor()
{
    iCursor = data.begin();
    return IStructure::resetCursor();
}
bool RectangularGrid::validCursor() const { return iCursor != data.end(); }
ElementData& RectangularGrid::cursorData() { return *iCursor; }
const ElementData& RectangularGrid::cursorData() const { return *iCursor; }
void RectangularGrid::incrCursor() { ++iCursor; }

} /* namespace Nextsim */
