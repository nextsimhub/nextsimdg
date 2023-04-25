/*!
 * @file MinimumIce.hpp
 *
 * @date Nov 17, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MINIMUMICE_HPP
#define MINIMUMICE_HPP

#include "include/Configured.hpp"

namespace Nextsim {

class MinimumIce : public Configured<MinimumIce> {
public:
    MinimumIce() = default;

    enum {
        MINC_KEY,
        MINH_KEY,
    };
    void configure() override;
    ConfigMap getConfiguration() const override;

    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    static inline double thickness() { return minh; }
    static inline double concentration() { return minc; }

private:
    static double minc; // Minimum sea ice concentration
    static double minh; // Minimum sea ice thickness
};

} /* namespace Nextsim */

#endif /* MINIMUMICE_HPP */
