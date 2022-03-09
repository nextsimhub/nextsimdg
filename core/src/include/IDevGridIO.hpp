/*!
 * @file IDevGridIO.hpp
 *
 * @date Jan 20, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_IDEVGRIDIO_HPP_
#define CORE_SRC_INCLUDE_IDEVGRIDIO_HPP_

#include "include/ElementData.hpp"
#include "include/ModelState.hpp"

namespace Nextsim {

class DevGrid;
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
     * @brief Generates the ModelState based on the data in the given file.
     *
     * @param filePath The location of the NetCDF restart file to be read.
     */
    virtual ModelState getModelState(const std::string& filePath) const = 0;

    /*!
     * @brief Writes data from the vector of data elements into the file location.
     *
     * @param dg The vector of ElementData instances containing the data.
     * @param filePath The location of the NetCDF restart file to be written.
     */
    virtual void dump(const std::vector<ElementData>& dg, const std::string& fielPath) const = 0;

    /*!
     * @brief Dumps the given ModelState to the given file path.
     *
     * @param state The ModelState data
     * @param filePath The path to attempt to write the data to.
     */
    virtual void dumpModelState(const ModelState& state, const std::string& filePath) const = 0;

protected:
    DevGrid* grid;
};

}
#endif /* CORE_SRC_INCLUDE_IDEVGRIDIO_HPP_ */
