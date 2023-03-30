/*!
 * @file UniformOcean.hpp
 *
 * @date 30 Mar 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef UNIFORMOCEAN_HPP
#define UNIFORMOCEAN_HPP

#include "IOceanBoundary.hpp"

namespace Nextsim {

//* Ocean boundary values that are constant with space and time.
class UniformOcean : public IOceanBoundary {
public:
    UniformOcean();

    std::string getName() const override { return "UniformOcean"; }
    void setData(const ModelState::DataMap&) override;
    void updateBefore(const TimestepTime&) override { }
    void updateAfter(const TimestepTime&) override { }
    // TODO ^add the SlabOcean when it becomes available

    UniformOcean& setSST(double);
    UniformOcean& setSSS(double);
    UniformOcean& setMLD(double);
    UniformOcean& setU(double);
    UniformOcean& setV(double);

private:
    static double sst0;
    static double sss0;
    static double mld0;
    static double u0;
    static double v0;
};

} /* namespace Nextsim */

#endif /* UNIFORMOCEAN_HPP */
