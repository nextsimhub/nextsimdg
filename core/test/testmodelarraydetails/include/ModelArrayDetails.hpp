/*!
 * @file ModelArrayDetails.hpp
 *
 * @date Oct 19, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODELARRAYDETAILS_HPP
#define MODELARRAYDETAILS_HPP

enum class Dimension { X, Y, Z, U, COUNT };

enum class Type {
    ONED,
    TWOD,
    ZUFIELD,
    THREED,
    FOURD,
};

static ModelArray OneDField() { return ModelArray(Type::ONED); }
static ModelArray TwoDField() { return ModelArray(Type::TWOD); }
static ModelArray ZUField() { return ModelArray(Type::ZUFIELD); }
static ModelArray ThreeDField() { return ModelArray(Type::THREED); }
static ModelArray FourDField() { return ModelArray(Type::FOURD); }

#endif /* MODELARRAYDETAILS_HPP */
