/*!
 * @file SlabOcean.cpp
 *
 * @date 27 Jan 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/SlabOcean.hpp"

namespace Nextsim {

const std::string SlabOcean::sstSlabName = "sst_slab";
const std::string SlabOcean::sssSlabName = "sss_slab";

void SlabOcean::setData(const ModelState::DataMap& ms)
{
    sstSlab.resize();
    sssSlab.resize();
    qdw.resize();
    fdw.resize();
}

ModelState SlabOcean::getState() const
{
    return { {
                 { sstSlabName, sstSlab },
                 { sssSlabName, sssSlab },
             },
        {} };
}
ModelState SlabOcean::getState(const OutputLevel&) const { return getState(); }

std::unordered_set<std::string> SlabOcean::hFields() const { return { sstSlabName, sssSlabName }; }

void SlabOcean::update(const TimestepTime&)
{

}

} /* namespace Nextsim */
