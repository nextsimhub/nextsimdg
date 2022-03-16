/*!
 * @file IVerticalIceGrowth.hpp
 *
 * @date Mar 16, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IVERTICALICEGROWTH_HPP
#define IVERTICALICEGROWTH_HPP

#include "include/ModelModule.hpp"

namespace Nextsim {
class IVerticalIceGrowth : public ModelModule {
public:
    ~IVerticalIceGrowth() = default;


protected:
    IVerticalIceGrowth() = default;
};

} /* namespace Nextsim */

#endif /* IVERTICALICEGROWTH_HPP */
