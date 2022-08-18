/*!
 * @file MissingData.hpp
 *
 * @date Jun 14, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MISSINGDATA_HPP
#define MISSINGDATA_HPP

#include "include/Configured.hpp"

namespace Nextsim {

class MissingData : public Configured<MissingData> {
public:
    inline static double value() { return m_value; }
    void configure() override;
    enum {
        MISSINGVALUE_KEY,
    };

    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll) { return getHelpText(map, getAll); }
    static HelpMap& getHelpText(HelpMap& map, bool getAll);
private:
    static double m_value;
};

} /* namespace Nextsim */

#endif /* MISSINGDATA_HPP */
