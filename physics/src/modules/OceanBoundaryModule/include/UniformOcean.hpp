/*!
 * @file UniformOcean.hpp
 *
 * @date 30 Mar 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef UNIFORMOCEAN_HPP
#define UNIFORMOCEAN_HPP

#include "include/IOceanBoundary.hpp"

namespace Nextsim {

//* Ocean boundary values that are constant with space and time.
class UniformOcean : public IOceanBoundary {
public:
    UniformOcean()
        // The same defaults as ConstantOceanBoundary
        : UniformOcean(-1.5, 32., 10, 0., 0.)
    {
    }
    UniformOcean(double sstIn, double sssIn, double mldIn, double uIn = 0., double vIn = 0.,
        double qioIn = 0.)
        : sst0(sstIn)
        , sss0(sssIn)
        , mld0(mldIn)
        , u0(uIn)
        , v0(vIn)
        , qio0(qioIn)
    {
    }

    std::string getName() const override { return "UniformOcean"; }
    void setData(const ModelState::DataMap&) override;
    void updateBefore(const TimestepTime&) override { }
    void updateAfter(const TimestepTime&) override { }
    // TODO ^add the SlabOcean when it becomes available

    UniformOcean& setSST(double);
    UniformOcean& setSSS(double);
    UniformOcean& setMLD(double);
    UniformOcean& setQio(double);
    UniformOcean& setU(double);
    UniformOcean& setV(double);

private:
    double sst0;
    double sss0;
    double mld0;
    double u0;
    double v0;
    double qio0;
};

} /* namespace Nextsim */

#endif /* UNIFORMOCEAN_HPP */
