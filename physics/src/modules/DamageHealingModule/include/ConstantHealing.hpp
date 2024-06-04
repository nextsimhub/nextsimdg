/*!
 * @file ConstantHealing.hpp
 *
 * @date Jun 3, 2024
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef CONSTANTHEALING_HPP
#define CONSTANTHEALING_HPP

#include "include/Configured.hpp"
#include "include/IDamageHealing.hpp"

namespace Nextsim {

//! A class implementing constant healing of damage (e.g. Bouillon and Rampal, 2015)
class ConstantHealing : public IDamageHealing, public Configured<ConstantHealing> {
public:
    ConstantHealing()
        : IDamageHealing()
    {
    }
    virtual ~ConstantHealing() = default;

    void configure() override;
    enum { TD_KEY };

    ModelState getStateRecursive(const OutputSpec& os) const override;

    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap&, bool getAll);

    void update(const TimestepTime& tstep) override;

private:
    static double tD;
};

}

#endif // CONSTANTHEALING_HPP
