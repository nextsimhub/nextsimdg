/*!
 * @file DGNetcdf.hpp
 *
 * @date Oct 4, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef DGNETCDF_HPP
#define DGNETCDF_HPP

#include "include/ModelArray.hpp"
#include <string>

namespace netCDF {
class NcFile;
class NcGroup;
}

namespace Nextsim {

/*!
 * A class to read and write ModelArrays with components from and to NetCDF
 * files. These are the DG and DGSTRESS Types.
 */
class DGNetcdf {
public:
    DGNetcdf();
    enum FileMode {
        READ,
        WRITE,
        REPLACE,
        NEWFILE,
    };
    DGNetcdf(const std::string& filePath, FileMode fm);
    ~DGNetcdf();

    // A user defined dtor invokes The Rule of Five. Here, delete them all
    DGNetcdf(const DGNetcdf&) = delete;
    DGNetcdf(const DGNetcdf&&) = delete;
    DGNetcdf& operator=(const DGNetcdf&) = delete;
    DGNetcdf& operator=(DGNetcdf&&) = delete;

    /*!
     * @brief Writes a ModelArray DGField to the referenced netCDF::NcGroup with the given name.
     *
     * @param group The netCDF group to write the data to.
     * @param dgField The ModelArray holding the data with DG components.
     * @param fieldName The name of the field to contain the data within the group.
     *
     * @return A reference to the netCDF group.
     */
    netCDF::NcGroup& write(netCDF::NcGroup& group, const ModelArray& dgField, const std::string& fieldName);

    /*!
     * @brief Reads a ModelArray DGField from the referenced netCDF::NcGroup with the given name.
     *
     * @param group The netCDF group to read the data from.
     * @param dgField The ModelArray to hold the data with DG components.
     * @param fieldName The name of the field containing the data within the group.
     *
     * @return A reference to the newly filled ModelArray.
     */
    ModelArray& read(netCDF::NcGroup& group, ModelArray& dgField, const std::string& fieldName);
private:
    netCDF::NcFile* p_file;
};

} /* namespace Nextsim */

#endif /* DGNETCDF_HPP */
