/*!
 * @file Logged.cpp
 * @date 12 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "Logged.hpp"

namespace Nextsim {

Logged::Logged()
{
    // TODO Auto-generated constructor stub
}

// TODO Empty implementations of the declared functions
void Logged::log(const std::string& message, Logged::level lvl)
{
    switch (lvl) {
    case (INFO):
        info(message);
        break;
    case (DEBUG):
    case (NOTICE):
    case (WARNING):
    case (ERROR):
    case (CRITICAL):
    case (ALERT):
    // TODO implement these levels
        break;
    case (EMERGENCY):
        emergency(message);
        break;
    default:
        break;
    }
}

void Logged::info(const std::string& message)
{
    // TODO Replace empty implementation
}

void Logged::emergency(const std::string& message)
{
    // TODO Replace empty implementation
}

} /* namespace Nextsim */
