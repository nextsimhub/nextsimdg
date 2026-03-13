/*!
 * @file Configured.hpp
 *
 * @date Mar 13, 2026
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
    template <typename T>
    static inline T getConfiguration(const std::string& name, const T& defaultValue) {
        return T();
    }
};

#endif /* MOCK_INCLUDES_CONFIGURED_HPP */
