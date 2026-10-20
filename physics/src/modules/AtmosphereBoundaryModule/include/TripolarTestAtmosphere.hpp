/*!
 * @author  Einar Ólason <einar.olason@nersc.no>
 */

#ifndef TRIPOLARTESTATMOSPHERE_HPP
#define TRIPOLARTESTATMOSPHERE_HPP

#include "include/IAtmosphereBoundary.hpp"

namespace Nextsim {

class TripolarTestAtmosphere : public IAtmosphereBoundary {
public:
    TripolarTestAtmosphere();

    std::string getName() const override { return "TripolarTestAtmosphere"; }
    void setData(const ModelState::DataMap&) override;
    void update(const TimestepTime& tst) override;
};

} /* namespace Nextsim */

#endif // TRIPOLARTESTATMOSPHERE_HPP
