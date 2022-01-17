/*!
 * @file DevGrid.cpp
 *
 * @date Dec 20, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/DevGrid.hpp"

#include <cstddef>
#include <vector>

namespace Nextsim {

const std::string DevGrid::ourStructureName = "devgrid";
const std::string DevGrid::xDimName = "x";
const std::string DevGrid::yDimName = "y";
const int DevGrid::nx = 10;

void DevGrid::init(const std::string& filePath)
{
    data.resize(nx * nx);
    if (pio && !filePath.empty()) {
        pio->init(data, filePath);
    }
};

void DevGrid::dump(const std::string& filePath) const
{
    if (pio  && !filePath.empty()) {
        pio->dump(data, filePath);
    }
};

// Cursor manipulation override functions
int DevGrid::resetCursor()
{
    iCursor = data.begin();
    return IStructure::resetCursor();
}
bool DevGrid::validCursor() const { return iCursor != data.end(); }
ElementData& DevGrid::cursorData() { return *iCursor; }
const ElementData& DevGrid::cursorData() const { return *iCursor; }
void DevGrid::incrCursor() { ++iCursor; }

IStructure& DevGrid::Cursor::operator=(const int i) const
{
    if (0 == i)
        owner.resetCursor();
    return owner;
}

DevGrid::Cursor::operator bool() const { return owner.validCursor(); }

ElementData& DevGrid::Cursor::operator*() const { return owner.cursorData(); }

ElementData* DevGrid::Cursor::operator->() const { return &owner.cursorData(); }

IStructure& DevGrid::Cursor::operator++() const
{
    owner.incrCursor();
    return owner;
}

} /* namespace Nextsim */
