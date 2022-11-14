/*!
 * @file ModelArrayDetails.hpp
 *
 * @date Oct 19, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODELARRAYDETAILS_HPP
#define MODELARRAYDETAILS_HPP

// An inclusion file to detail the ModelArray dimensions and types for
// Discontinuous Galerkin models.

// Should be grouped with a consistent ModelArrayTypedefs.hpp and
// ModelArrayDetails.cpp

enum class Dimension { X, Y, Z, XVERTEX, YVERTEX, XCG, YCG, DG, DGSTRESS, NCOORDS, COUNT };

enum class Type {
    H,
    VERTEX,
    U,
    V,
    Z,
    DG,
    DGSTRESS,
    CG,
};

static ModelArray HField() { return ModelArray(Type::H); }
static ModelArray VertexField() { return ModelArray(Type::VERTEX); }
static ModelArray UField() { return ModelArray(Type::U); }
static ModelArray VField() { return ModelArray(Type::V); }
static ModelArray ZField() { return ModelArray(Type::Z); }
static ModelArray DGField() { return ModelArray(Type::DG); }
static ModelArray DGSField() { return ModelArray(Type::DGSTRESS); }
static ModelArray CGField() { return ModelArray(Type::CG); }

const static size_t nCoords;

#endif /* MODELARRAYDETAILS_HPP */
