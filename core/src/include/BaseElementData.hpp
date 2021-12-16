/*!
 * @file BaseElementData.hpp
 * @date Sep 8, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_BASEELEMENTDATA_HPP
#define SRC_INCLUDE_BASEELEMENTDATA_HPP

namespace Nextsim {
/*!
 * @brief The base class for the per element data classes.
 *
 * @details This base class handles the common features of the per-element
 * data classes. These include handling output.
 */
class BaseElementData {
public:
    BaseElementData() = default;
    // TODO: implement output handling
    ~BaseElementData() = default;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_BASEELEMENTDATA_HPP */
