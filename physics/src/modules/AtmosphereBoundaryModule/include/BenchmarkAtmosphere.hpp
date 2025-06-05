/*!
 * @file BenchmarkAtmosphere.hpp
 *
 * @date 24 Sep 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef BENCHMARKATMOSPHERE_HPP
#define BENCHMARKATMOSPHERE_HPP

#include "include/IAtmosphereBoundary.hpp"

namespace Nextsim {

class BenchmarkAtmosphere : public IAtmosphereBoundary {
public:
    BenchmarkAtmosphere()
        : IAtmosphereBoundary()
    {
    }
    ~BenchmarkAtmosphere() = default;

    void setData(const ModelState::DataMap&) override;
    std::string getName() const override { return "BenchmarkAtmosphere"; }

    void update(const TimestepTime& tst) override;

private:
    TimePoint t0;
    bool t0Set = false;
};

} /* namespace Nextsim */

#endif /* BENCHMARKATMOSPHERE_HPP */
