/*!
 * @file ModelArrayDetails.hpp
 *
 * @date Oct 19, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODELARRAYDETAILS_HPP
#define MODELARRAYDETAILS_HPP

enum class Dimension { X, Y, Z, XCG, YCG, DG, DGSTRESS, COUNT };

enum class Type {
    H,
    U,
    V,
    Z,
    DG,
    DGSTRESS,
    CG,
};

static ModelArray HField() { return ModelArray(Type::H); }
static ModelArray UField() { return ModelArray(Type::U); }
static ModelArray VField() { return ModelArray(Type::V); }
static ModelArray ZField() { return ModelArray(Type::Z); }
static ModelArray DGField() { return ModelArray(Type::DG); }
static ModelArray DGSField() { return ModelArray(Type::DGSTRESS); }
static ModelArray CGField() { return ModelArray(Type::CG); }

#endif /* MODELARRAYDETAILS_HPP */
