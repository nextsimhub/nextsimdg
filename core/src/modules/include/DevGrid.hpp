/*!
 * @file DevGrid.hpp
 *
 * @date Dec 20, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_DEVGRID_HPP
#define CORE_SRC_INCLUDE_DEVGRID_HPP

#include "include/IStructure.hpp"

#include "include/ElementData.hpp"
#include "include/PrognosticData.hpp"

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

    //! Destructor. The lifetime of pio should be the lifetime of the instance.
    virtual ~DevGrid()
    {
        if (pio) {
            delete pio;
        }
    }

    const static int nx;

    // Read/write override functions
    void init(const std::string& filePath) override;

    void dump(const std::string& filePath) const override;

    std::string structureType() const override { return ourStructureName; };

    // Cursor manipulation override functions
    int resetCursor() override;
    bool validCursor() const override;
    ElementData& cursorData() override;
    const ElementData& cursorData() const override;
    void incrCursor() override;

    /*!
     * @brief A class that deals with all the netCDF related parts of DevGrid.
     *
     * @details DevGrid is a module, and currently merely mentioning it in the
     * code requires the cpp file to be linked. Since this would pull in the
     * NetCDF libraries, the DevGrid::IDevGridIO interface was created to
     * separate out all the code that actually uses the NetCDF libraries. See
     * DevGridIO for the implementing class.
     */
    class IDevGridIO {
    public:
        IDevGridIO(DevGrid& grid)
            : grid(&grid)
        {
        }
        virtual ~IDevGridIO() = default;
        /*!
         * @brief Reads data from the file location into the vector of data elements.
         *
         * @param dg The vector of ElementData instances to be filled.
         * @param filePath The location of the NetCDF restart file to be read.
         */
        virtual void init(std::vector<ElementData>& dg, const std::string& filePath) const = 0;
        /*!
         * @brief Writes data from the vector of data elements into the file location.
         *
         * @param dg The vector of ElementData instances containing the data.
         * @param filePath The location of the NetCDF restart file to be written.
         */
        virtual void dump(
            const std::vector<ElementData>& dg, const std::string& fielPath) const = 0;

    protected:
        DevGrid* grid;
    };

    //! Sets the pointer to the class that will perform the IO. Should be an instance of DevGridIO
    void setIO(IDevGridIO* p) { pio = p; }

private:
    const static std::string ourStructureName;
    const static std::string xDimName;
    const static std::string yDimName;

    std::vector<ElementData> data;

    std::vector<ElementData>::iterator iCursor;

    IDevGridIO* pio;

    friend DevGridIO;
};

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_DEVGRID_HPP */
