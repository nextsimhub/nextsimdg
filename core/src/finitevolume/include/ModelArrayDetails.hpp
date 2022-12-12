/*!
 * @file ModelArrayDetails.hpp
 *
 * @date Oct 19, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODELARRAYDETAILS_HPP
#define MODELARRAYDETAILS_HPP

enum class Dimension { X, Y, Z, COUNT };

enum class Type {
    H,
    U,
    V,
    Z,
};

static ModelArray HField() { return ModelArray(Type::H); }
static ModelArray UField() { return ModelArray(Type::U); }
static ModelArray VField() { return ModelArray(Type::V); }
static ModelArray ZField() { return ModelArray(Type::Z); }

#endif /* MODELARRAYDETAILS_HPP */
