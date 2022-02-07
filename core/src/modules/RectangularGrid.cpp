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

int RectangularGrid::nx;
int RectangularGrid::ny;
int RectangularGrid::nz;

void RectangularGrid::init(const std::string& filePath)
{
    data.resize(nx * ny);
    if (pio && !filePath.empty()) {
        pio->init(data, filePath);
    }
}

void RectangularGrid::dump(const std::string& filePath) const
{
    if (pio && !filePath.empty()) {
        pio->dump(data, filePath);
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
