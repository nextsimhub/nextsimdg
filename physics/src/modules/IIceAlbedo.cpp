#include "include/IIceAlbedo.hpp"

namespace Nextsim {

IIceAlbedo::IIceAlbedo()
    : iceAlbedoAccessor(getStore(), RO, ModelArray::Type::H)
    , icePenSWAccessor(getStore(), RO, ModelArray::Type::H)
    , hsnowAccessor(getStore())
    , ciceAccessor(getStore())
    , tsurfAccessor(getStore())
{
}

void IIceAlbedo::setData(const ModelState::DataMap&)
{
    iceAlbedoAccessor.getHostRW().reinitialize();
    icePenSWAccessor.getHostRW().reinitialize();
}
}
