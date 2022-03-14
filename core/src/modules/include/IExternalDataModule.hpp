/*!
 * @file IExternalDataModule.hpp
 *
 * @date Mar 14, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_MODULES_INCLUDE_IEXTERNALDATAMODULE_HPP
#define CORE_SRC_MODULES_INCLUDE_IEXTERNALDATAMODULE_HPP

#include "include/Module.hpp"

#include "include/IExternalData.hpp"

namespace Module {

template<> Module <Nextsim::IExternalData>::map Module<Nextsim::IExternalData>::functionMap;
class IExternalDataModule : public Module<Nextsim::IExternalData> {
public:
    struct Constructor {
        Constructor();
    };
    static Constructor ctor;
};

} /* namespace Module */

#endif /* CORE_SRC_MODULES_INCLUDE_IEXTERNALDATAMODULE_HPP */
