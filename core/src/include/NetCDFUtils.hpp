#ifndef NETCDFUTILS_HPP
#define NETCDFUTILS_HPP

#include <ncDouble.h>
#include <ncFloat.h>
#include <ncInt64.h>

namespace Nextsim {

template <typename T> struct ToNetCDFType;

template <> struct ToNetCDFType<double> {
    static netCDF::NcDouble& get() { return netCDF::ncDouble; }
};

template <> struct ToNetCDFType<float> {
    static netCDF::NcFloat& get() { return netCDF::ncFloat; }
};

}

#endif
