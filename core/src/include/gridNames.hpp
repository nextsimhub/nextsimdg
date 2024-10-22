/*!
 * @file gridNames.hpp
 *
 * @date 07 Oct 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef GRIDNAMES_HPP
#define GRIDNAMES_HPP

#include <string>

namespace Nextsim {

// A set of canonical names for fields used in file reading and writing.
static const std::string hiceName = "hice";
static const std::string ciceName = "cice";
static const std::string hsnowName = "hsnow";
static const std::string ticeName = "tice";
static const std::string maskName = "mask";
static const std::string uName = "u";
static const std::string vName = "v";
static const std::string sstName = "sst";
static const std::string sssName = "sss";
static const std::string damageName = "damage";
static const std::string shearName = "shear";

static const std::string uWindName = "uwind";
static const std::string vWindName = "vwind";
static const std::string uOceanName = "uocean";
static const std::string vOceanName = "vocean";
static const std::string sshName = "ssh";
// Mixed layer depth
static const std::string mldName = "mld";

static const std::string uIOStressName = "uiostress";
static const std::string vIOStressName = "viostress";

static const std::string coordsName = "coords";
static const std::string latitudeName = "latitude";
static const std::string longitudeName = "longitude";
static const std::string xName = "x";
static const std::string yName = "y";
static const std::string gridAzimuthName = "grid_azimuth";

static const std::string mdiName = "missing_value";

static const std::string timeName = "time";

}

#endif /* GRIDNAMES_HPP */
