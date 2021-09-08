/*!
 * @file BaseElementData.hpp
 * @date Sep 8, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_BASEELEMENTDATA_HPP
#define SRC_INCLUDE_BASEELEMENTDATA_HPP

namespace Nextsim {
/*!
 * @brief The base class to the per element data classes
 *
 * @detailed This base class will handle the common features of the per-element
 * data classes. These include handling output.
 */
class BaseElementData {
public:
    BaseElementData();
    // TODO: implement output handling
    virtual ~BaseElementData();
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_BASEELEMENTDATA_HPP */
