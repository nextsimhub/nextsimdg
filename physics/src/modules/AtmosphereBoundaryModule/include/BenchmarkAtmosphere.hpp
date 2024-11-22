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

class BenchmarkAtmosphere : public IAtmosphereBoundary, public Configured<IAtmosphereBoundary> {
public:
    BenchmarkAtmosphere()
        : IAtmosphereBoundary()
        , t0Set(false)
    {
    }
    ~BenchmarkAtmosphere() = default;

    enum { SPINUP_KEY };

    void setData(const ModelState::DataMap&) override;
    std::string getName() const override { return "BenchmarkAtmosphere"; }

    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    void configure() override;

    void update(const TimestepTime& tst) override;

private:
    TimePoint t0;
    bool t0Set;

    static double spinupTime;
};

} /* namespace Nextsim */

#endif /* BENCHMARKATMOSPHERE_HPP */
