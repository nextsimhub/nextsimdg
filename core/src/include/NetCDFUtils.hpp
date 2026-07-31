#ifndef NETCDFUTILS_HPP
#define NETCDFUTILS_HPP

#include <ncDouble.h>
#include <ncFloat.h>
#include <ncInt64.h>
#include <ncVar.h>

#include <numeric>
#include <stdexcept>
#include <vector>

namespace Nextsim {

template <typename T> struct ToNetCDFType;

template <> struct ToNetCDFType<double> {
    static netCDF::NcDouble& get() { return netCDF::ncDouble; }
};

template <> struct ToNetCDFType<float> {
    static netCDF::NcFloat& get() { return netCDF::ncFloat; }
};

namespace Details {
    template <typename ReadT, typename DestT>
    void readConvertNetCDFVar(const netCDF::NcVar& var, const std::vector<size_t>& start,
        const std::vector<size_t>& size, DestT* dest)
    {
        const size_t numElem
            = std::reduce(size.begin(), size.end(), static_cast<size_t>(1), std::multiplies<>());
        std::vector<ReadT> buf(numElem, 0.f);
        var.getVar(start, size, buf.data());
        DestT offset;
        DestT scale;
        try {
            const netCDF::NcVarAtt scaleAtt = var.getAtt("scale_factor");
            const netCDF::NcVarAtt offsetAtt = var.getAtt("add_offset");
            scaleAtt.getValues(&scale);
            offsetAtt.getValues(&offset);
        } catch (const netCDF::exceptions::NcException&) {
            // Ignore missing attributes
            offset = 0;
            scale = 1;
        }
        for (size_t i = 0; i < numElem; ++i) {
            dest[i] = static_cast<DestT>(buf[i] * scale + offset);
        }
    }

}

template <typename T>
void readNetCDFVar(const netCDF::NcVar& var, const std::vector<size_t>& start,
    const std::vector<size_t>& size, T* dest)
{
    static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>,
        "Currently only conversion to floating point types is supported during load.");

    const netCDF::NcType ncType = var.getType();

    if (ncType == ToNetCDFType<T>::get()) {
        var.getVar(start, size, dest);
    } else {
        if (ncType == netCDF::ncFloat) {
            Details::readConvertNetCDFVar<float>(var, start, size, dest);
        } else if (ncType == netCDF::ncDouble) {
            Details::readConvertNetCDFVar<double>(var, start, size, dest);
        } else if (ncType == netCDF::ncShort) {
            Details::readConvertNetCDFVar<short>(var, start, size, dest);
        } else {
            throw std::domain_error("Unsupported type of input field " + var.getName() + ".\n");
        }
    }
}

}

#endif
