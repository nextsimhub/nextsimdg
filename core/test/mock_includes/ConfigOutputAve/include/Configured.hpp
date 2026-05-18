/*!
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MOCK_INCLUDES_CONFIGURED_HPP
#define MOCK_INCLUDES_CONFIGURED_HPP

#include "include/ConfigMap.hpp"

template <typename C>
class Configured {
public:
    virtual ~Configured() = default;
    virtual void configure() = 0;
    static std::string getConfiguration(const std::string& name, const std::string& defaultValue) {
        const std::string pfx = "ConfigOutput";
        if (name == pfx+".period") {
            return std::string("3600");
        } else if (name == pfx+".start") {
            return std::string("2000-01-01T00:00:00Z");
        } else if (name == pfx+".field_names") {
            return std::string("cice");
        } else if (name == pfx+".filename") {
            return "";
        } else if (name == pfx+".file_period") {
            return std::string("10000000000");
        } else {
            return std::string();
        }
    }
};

#endif /* MOCK_INCLUDES_CONFIGURED_HPP */
