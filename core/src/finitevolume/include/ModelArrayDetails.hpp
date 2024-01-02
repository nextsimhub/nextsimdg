/*!
 * @file ModelArrayDetails.hpp
 *
 * @date Oct 19, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODELARRAYDETAILS_HPP
#define MODELARRAYDETAILS_HPP

// An inclusion file to detail the ModelArray dimensions and types for
// finite volume models.

// Should be grouped with a consistent ModelArrayTypedefs.hpp and
// ModelArrayDetails.cpp

enum class Dimension { X, Y, Z, XVERTEX, YVERTEX, COUNT };

enum class Type {
    H,
    U,
    V,
    Z,
    VERTEX,
};

static ModelArray HField() { return ModelArray(Type::H); }
static ModelArray UField() { return ModelArray(Type::U); }
static ModelArray VField() { return ModelArray(Type::V); }
static ModelArray ZField() { return ModelArray(Type::Z); }
static ModelArray VertexField() { return ModelArray(Type::VERTEX); }

const static size_t nCoords;

#endif /* MODELARRAYDETAILS_HPP */
