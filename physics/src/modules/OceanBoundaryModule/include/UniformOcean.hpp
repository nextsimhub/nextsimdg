/*!
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
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
    UniformOcean(FloatType sstIn, FloatType sssIn, FloatType mldIn, FloatType uIn = 0.,
        FloatType vIn = 0., FloatType qioIn = 0.)
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

    UniformOcean& setSST(FloatType);
    UniformOcean& setSSS(FloatType);
    UniformOcean& setMLD(FloatType);
    UniformOcean& setQio(FloatType);
    UniformOcean& setU(FloatType);
    UniformOcean& setV(FloatType);

private:
    FloatType sst0;
    FloatType sss0;
    FloatType mld0;
    FloatType u0;
    FloatType v0;
    FloatType qio0;
};

} /* namespace Nextsim */

#endif /* UNIFORMOCEAN_HPP */
