/*!
 * @file BenchmarkOcean.hpp
 *
 * @date 19 Apr 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef BENCHMARKOCEAN_HPP
#define BENCHMARKOCEAN_HPP

#include "IOceanBoundary.hpp"

namespace Nextsim {

class BenchmarkOcean : public IOceanBoundary {
public:
    BenchmarkOcean() = default;
    ~BenchmarkOcean() = default;

    void setData(const ModelState::DataMap&) override;
    std::string getName() const override { return "BenchmarkOcean"; }

    void updateBefore(const TimestepTime& tst) override;
    void updateAfter(const TimestepTime& tst) override { }
};

} /* namespace Nextsim */

#endif /* BENCHMARKOCEAN_HPP */
