/*!
 * @file IDevGridIO.hpp
 *
 * @date Jan 20, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IDEVGRIDIO_HPP
#define IDEVGRIDIO_HPP

#include "include/ModelMetadata.hpp"
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
     * @brief Generates the ModelState based on the data in the given file.
     *
     * @param filePath The location of the NetCDF restart file to be read.
     */
    virtual ModelState getModelState(const std::string& filePath) const = 0;

    /*!
     * @brief Dumps the given ModelState to the given file path.
     *
     * @param state The ModelState data
     * @param filePath The path to attempt to write the data to.
     * @param isRestart Should this file be written as a restart file or a
     *          diagnostic dump?
     */
    virtual void dumpModelState(const ModelState& state, const ModelMetadata& metadata,
        const std::string& filePath, bool isRestart) const
        = 0;

protected:
    DevGrid* grid;
};

}
#endif /* IDEVGRIDIO_HPP */
