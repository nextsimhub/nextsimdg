/*!
 * @file Fusi.hpp
 *
 * @date Aug 19, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef FUSI_HPP
#define FUSI_HPP

#include <cstdint>
#include <string>

/*!
 * A struct to act as a holder for floating point, unsigned int, string or integer data.
 * Fusi is an acronym of the types and means poured from in Latin as well as a
 * few melting related terms in other languages.
 */
struct Fusi {
    //    union Data {
    double f;
    std::int64_t i;
    std::uint64_t u;
    std::string s;
    //    } data;
    enum class Type { DOUBLE, INT, UINT, STRING, NONE } type;
    Fusi()
        : type(Type::NONE)
        , f(0.)
        , i(0)
        , u(0)
        , s()
    {
    }
    ~Fusi() = default;

    Fusi(int ii)
        : Fusi(static_cast<int64_t>(ii))
    {
    }

    Fusi(unsigned int uu)
        : Fusi(static_cast<uint64_t>(uu))
    {
    }

    Fusi(int64_t ii)
        : Fusi()
    {
        i = ii;
        type = Type::INT;
    }

    Fusi(uint64_t uu)
        : Fusi()
    {
        u = uu;
        type = Type::UINT;
    }

    Fusi(double ff)
        : Fusi()
    {
        f = ff;
        type = Type::DOUBLE;
    }

    Fusi(std::string ss)
        : Fusi()
    {
        s = ss;
        type = Type::STRING;
    }
};

#endif /* FUSI_HPP */
