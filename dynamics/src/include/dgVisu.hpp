/*!
 * @file dgVisu.hpp
 * @date 1 Mar 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __DGVISU_HPP
#define __DGVISU_HPP

#include "ParametricMesh.hpp"
#include "cgVector.hpp"
#include "dgVector.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace Nextsim {

class VTK {
public:
    //! Puts together the file name by adding the index n with fixed width
    static std::string compose_vtkname(const std::string& fname, int n)
    {
        std::ostringstream ss;
        ss << fname << "." << std::setw(5) << std::setfill('0') << n << ".vtk";
        return ss.str();
    }

    //! Puts together the file name by adding the DGdegree and index n with fixed width
    static std::string compose_vtkname(const std::string& fname, int DG, int n)
    {
        std::ostringstream ss;
        ss << fname << DG << "." << std::setw(5) << std::setfill('0') << n << ".vtk";
        return ss.str();
    }

    static inline double cuttol(double v, double tol)
    {
        if ((v < tol) && (v > -tol))
            return 0;
        return v;
    }

    ////////////////////////////////////////////////// Output of CG - Vector (ParametricMesh)
    template <int CG>
    static void write_cgvector(
        const std::string& fname, const CGVector<CG>& v, const ParametricMesh& smesh)
    {
        // extract variable name
        std::string variableName
            = fname.substr(fname.find("/") + 1, fname.find_first_of(".") - fname.find("/") - 1);

        std::ofstream OUT(fname.c_str());
        if (!OUT.is_open()) {
            std::cerr << "Failed to open '" << fname << "'." << std::endl;
            assert(0);
        }

        // Structure Points
        OUT << "# vtk DataFile Version 2.0" << std::endl
            << "output generated by Nextsim (ParametricMesh)" << std::endl
            << "ASCII" << std::endl
            << "DATASET UNSTRUCTURED_GRID" << std::endl
            << "FIELD FieldData 1" << std::endl
            << "TIME 1 1 double" << std::endl
            << "0" << std::endl;
        OUT << "POINTS " << v.rows() << " DOUBLE" << std::endl;
        if (CG == 1) {
            assert(smesh.nnodes == v.rows());
            assert((smesh.nx + 1) * (smesh.ny + 1) == v.rows());
            for (size_t i = 0; i < smesh.nnodes; ++i)
                write_coords(
                    OUT, smesh.vertices(i, 0), smesh.vertices(i, 1), smesh.CoordinateSystem);

            OUT << "CELLS " << smesh.nelements << " " << smesh.nelements * 5 << std::endl;
            for (size_t iy = 0; iy < smesh.ny; ++iy)
                for (size_t ix = 0; ix < smesh.nx; ++ix)
                    OUT << "4"
                        << " " << iy * (smesh.nx + 1) + ix << " " << iy * (smesh.nx + 1) + ix + 1
                        << " " << (iy + 1) * (smesh.nx + 1) + ix + 1 << " "
                        << (iy + 1) * (smesh.nx + 1) + ix << std::endl;

            OUT << "CELL_TYPES " << smesh.nelements << std::endl;
            for (size_t i = 0; i < smesh.nelements; ++i)
                OUT << "9 ";
            OUT << std::endl;

            OUT << "POINT_DATA " << v.rows() << std::endl;
        } else if (CG == 2) {
            assert(static_cast<int>((2 * smesh.nx + 1) * (2 * smesh.ny + 1)) == v.rows());
            for (size_t iy = 0; iy < 2 * smesh.ny + 1; ++iy)
                for (size_t ix = 0; ix < 2 * smesh.nx + 1; ++ix) {
                    const auto v = smesh.coordinate<2>(ix, iy);
                    OUT << v[0] << " " << v[1] << " 0" << std::endl;
                }

            OUT << std::endl
                << "CELLS " << 4 * smesh.nelements << " " << 4 * smesh.nelements * 5 << std::endl;
            for (size_t iy = 0; iy < 2 * smesh.ny; ++iy)
                for (size_t ix = 0; ix < 2 * smesh.nx; ++ix) {
                    const size_t n0 = (2 * smesh.nx + 1) * iy + ix;
                    OUT << "4"
                        << " " << n0 << " " << n0 + 1 << " " << n0 + 2 * smesh.nx + 1 + 1 << " "
                        << n0 + 2 * smesh.nx + 1 << std::endl;
                }

            OUT << std::endl << "CELL_TYPES " << 4 * smesh.nelements << std::endl;
            for (size_t i = 0; i < 4 * smesh.nelements; ++i)
                OUT << "9 ";
            OUT << std::endl;

            OUT << "POINT_DATA " << v.rows() << std::endl;
        } else {
            std::cerr << "write_cgvector only for CG1 or CG2" << std::endl;
            abort();
        }

        OUT << "SCALARS v DOUBLE" << std::endl << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < v.rows(); ++i)
            OUT << cuttol(v(i), 1.e-20) << std::endl;
        OUT.close();
    }

    template <int CG>
    static void write_cg_velocity(const std::string& fname, const CGVector<CG>& vx,
        const CGVector<CG>& vy, const ParametricMesh& smesh)
    {
        assert(vx.rows() == vy.rows());
        // extract variable name
        std::string variableName
            = fname.substr(fname.find("/") + 1, fname.find_first_of(".") - fname.find("/") - 1);

        std::ofstream OUT(fname.c_str());
        if (!OUT.is_open()) {
            std::cerr << "Failed to open '" << fname << "'." << std::endl;
            assert(0);
        }

        // Structure Points
        OUT << "# vtk DataFile Version 2.0" << std::endl
            << "output generated by Nextsim (ParametricMesh)" << std::endl
            << "ASCII" << std::endl
            << "DATASET UNSTRUCTURED_GRID" << std::endl
            << "FIELD FieldData 1" << std::endl
            << "TIME 1 1 double" << std::endl
            << "0" << std::endl;
        OUT << "POINTS " << vx.rows() << " DOUBLE" << std::endl;
        if (CG == 1) {
            assert(static_cast<int>(smesh.nnodes) == vx.rows());
            assert(static_cast<int>((smesh.nx + 1) * (smesh.ny + 1)) == vx.rows());
            for (size_t i = 0; i < smesh.nnodes; ++i)
                write_coords(
                    OUT, smesh.vertices(i, 0), smesh.vertices(i, 1), smesh.CoordinateSystem);

            OUT << "CELLS " << smesh.nelements << " " << smesh.nelements * 5 << std::endl;
            for (size_t iy = 0; iy < smesh.ny; ++iy)
                for (size_t ix = 0; ix < smesh.nx; ++ix)
                    OUT << "4"
                        << " " << iy * (smesh.nx + 1) + ix << " " << iy * (smesh.nx + 1) + ix + 1
                        << " " << (iy + 1) * (smesh.nx + 1) + ix + 1 << " "
                        << (iy + 1) * (smesh.nx + 1) + ix << std::endl;

            OUT << "CELL_TYPES " << smesh.nelements << std::endl;
            for (size_t i = 0; i < smesh.nelements; ++i)
                OUT << "9 ";
            OUT << std::endl;

            OUT << "POINT_DATA " << vx.rows() << std::endl;
        } else if (CG == 2) {
            assert(static_cast<int>((2 * smesh.nx + 1) * (2 * smesh.ny + 1)) == vx.rows());
            for (size_t iy = 0; iy < 2 * smesh.ny + 1; ++iy)
                for (size_t ix = 0; ix < 2 * smesh.nx + 1; ++ix) {
                    const std::array<double, 2> v = smesh.coordinate<2>(ix, iy);
                    write_coords(OUT, v[0], v[1], smesh.CoordinateSystem);
                }

            OUT << std::endl
                << "CELLS " << 4 * smesh.nelements << " " << 4 * smesh.nelements * 5 << std::endl;
            for (size_t iy = 0; iy < 2 * smesh.ny; ++iy)
                for (size_t ix = 0; ix < 2 * smesh.nx; ++ix) {
                    const size_t n0 = (2 * smesh.nx + 1) * iy + ix;
                    OUT << "4"
                        << " " << n0 << " " << n0 + 1 << " " << n0 + 2 * smesh.nx + 1 + 1 << " "
                        << n0 + 2 * smesh.nx + 1 << std::endl;
                }

            OUT << std::endl << "CELL_TYPES " << 4 * smesh.nelements << std::endl;
            for (size_t i = 0; i < 4 * smesh.nelements; ++i)
                OUT << "9 ";
            OUT << std::endl;

            OUT << "POINT_DATA " << vx.rows() << std::endl;
        } else {
            std::cerr << "write_cgvector only for CG1 or CG2" << std::endl;
            abort();
        }

        OUT << "SCALARS vx DOUBLE" << std::endl << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < vx.rows(); ++i)
            OUT << cuttol(vx(i), 1.e-20) << std::endl;
        OUT << "SCALARS vy DOUBLE" << std::endl << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < vy.rows(); ++i)
            OUT << cuttol(vy(i), 1.e-20) << std::endl;

        OUT << "VECTORS v DOUBLE" << std::endl;
        for (int i = 0; i < vy.rows(); ++i)
            OUT << cuttol(vx(i), 1.e-20) << " " << cuttol(vy(i), 1.e-20) << " 0" << std::endl;

        OUT.close();
    }

    static void write_coords(std::ostream& OUT, double x, double y, COORDINATES CoordinateSystem)
    {
        if (CoordinateSystem == SPHERICAL) {
            OUT << Nextsim::EarthRadius * cos(y) * cos(x) << "\t"
                << Nextsim::EarthRadius * cos(y) * sin(x) << "\t" << Nextsim::EarthRadius * sin(y)
                << std::endl;
        } else
            OUT << x << "\t" << y << "\t0" << std::endl;
    }

    template <int DG>
    static void write_dg(
        const std::string& fname, const DGVector<DG>& v, const ParametricMesh& smesh)
    {
        // extract variable name
        std::string variableName
            = fname.substr(fname.find("/") + 1, fname.find_first_of(".") - fname.find("/") - 1);

        std::ofstream OUT(fname.c_str());
        assert(OUT.is_open());

        // Structure Points
        OUT << "# vtk DataFile Version 2.0" << std::endl
            << "dg(" << DG << ") output generated by Nextsim" << std::endl
            << "ASCII" << std::endl
            << "DATASET UNSTRUCTURED_GRID" << std::endl
            << "FIELD FieldData 1" << std::endl
            << "TIME 1 1 double" << std::endl
            << 0.0 << std::endl;
        OUT << "POINTS " << 4 * smesh.nx * smesh.ny << " DOUBLE" << std::endl;

        // local shift for indices within the element
        const size_t LEL[4] = { 0, 1, smesh.nx + 2, smesh.nx + 1 };

        size_t nid = 0; // id of first node
        for (size_t iy = 0; iy < smesh.ny; ++iy, ++nid) // also increase nid to add +1
            for (size_t ix = 0; ix < smesh.nx; ++ix, ++nid) {
                for (size_t j = 0; j < 4; ++j) // go along the 4 nodes of the element
                {
                    write_coords(OUT, smesh.vertices(nid + LEL[j], 0),
                        smesh.vertices(nid + LEL[j], 1), smesh.CoordinateSystem);
                }
            }
        OUT << "CELLS " << smesh.nx * smesh.ny << " " << 5 * smesh.nx * smesh.ny << std::endl;
        size_t ii = 0;
        for (size_t iy = 0; iy < smesh.ny; ++iy)
            for (size_t ix = 0; ix < smesh.nx; ++ix, ++ii)
                OUT << "4 " << 4 * ii << " " << 4 * ii + 1 << " " << 4 * ii + 2 << " " << 4 * ii + 3
                    << std::endl;

        OUT << "CELL_TYPES " << smesh.nx * smesh.ny << std::endl;
        for (ii = 0; ii < smesh.nx * smesh.ny; ++ii)
            OUT << "9 ";
        OUT << std::endl;

        OUT << "POINT_DATA " << 4 * smesh.nx * smesh.ny << std::endl
            << "SCALARS " << variableName << " DOUBLE " << std::endl
            << "LOOKUP_TABLE default" << std::endl;

        ii = 0;
        for (size_t iy = 0; iy < smesh.ny; ++iy)
            for (size_t ix = 0; ix < smesh.nx; ++ix, ++ii) {
                std::array<double, 4> interpolate = { v(ii, 0), v(ii, 0), v(ii, 0), v(ii, 0) };

                if (DG >= 3) {
                    interpolate[0] += -0.5 * v(ii, 1);
                    interpolate[1] += 0.5 * v(ii, 1);
                    interpolate[2] += 0.5 * v(ii, 1);
                    interpolate[3] += -0.5 * v(ii, 1);

                    interpolate[0] += -0.5 * v(ii, 2);
                    interpolate[1] += -0.5 * v(ii, 2);
                    interpolate[2] += 0.5 * v(ii, 2);
                    interpolate[3] += 0.5 * v(ii, 2);
                }
                if (DG >= 6) {
                    interpolate[0] += 1. / 6. * v(ii, 3);
                    interpolate[1] += 1. / 6. * v(ii, 3);
                    interpolate[2] += 1. / 6. * v(ii, 3);
                    interpolate[3] += 1. / 6. * v(ii, 3);

                    interpolate[0] += 1. / 6. * v(ii, 4);
                    interpolate[1] += 1. / 6. * v(ii, 4);
                    interpolate[2] += 1. / 6. * v(ii, 4);
                    interpolate[3] += 1. / 6. * v(ii, 4);

                    interpolate[0] += 1. / 4. * v(ii, 5);
                    interpolate[1] += -1. / 4. * v(ii, 5);
                    interpolate[2] += 1. / 4. * v(ii, 5);
                    interpolate[3] += -1. / 4. * v(ii, 5);
                }
                for (auto it : interpolate)
                    OUT << cuttol(it, 1.e-20) << std::endl;
            }
        OUT.close();
    }

    static void write_dg(
        const std::string& fname, const DGVector<6>& v, const ParametricMesh& smesh)
    {
        // extract variable name
        std::string variableName
            = fname.substr(fname.find("/") + 1, fname.find_first_of(".") - fname.find("/") - 1);

        std::ofstream OUT(fname.c_str());
        assert(OUT.is_open());

        // Structure Points
        OUT << "# vtk DataFile Version 2.0" << std::endl
            << "dg(" << 6 << ") output generated by Nextsim" << std::endl
            << "ASCII" << std::endl
            << "DATASET UNSTRUCTURED_GRID" << std::endl
            << "FIELD FieldData 1" << std::endl
            << "TIME 1 1 double" << std::endl
            << 0.0 << std::endl;
        OUT << "POINTS " << 3 * 3 * smesh.nx * smesh.ny << " DOUBLE"
            << std::endl; // add substructre
        size_t nid = 0; // id of first node

        const size_t sy = smesh.nx + 1; // shift one line up

        for (size_t iy = 0; iy < smesh.ny; ++iy, ++nid) // also increase nid to add +1
            for (size_t ix = 0; ix < smesh.nx; ++ix, ++nid) {
                // Print 9 points of the element incl. subsamples
                const Eigen::Matrix<Nextsim::FloatType, 4, 2> coords
                    = smesh.coordinatesOfElement(smesh.nx * iy + ix);

                write_coords(OUT, coords(0, 0), coords(0, 1), smesh.CoordinateSystem);
                write_coords(OUT, 0.5 * (coords(0, 0) + coords(1, 0)),
                    0.5 * (coords(0, 1) + coords(1, 1)), smesh.CoordinateSystem);
                write_coords(OUT, coords(1, 0), coords(1, 1), smesh.CoordinateSystem);

                write_coords(OUT, 0.5 * (coords(0, 0) + coords(2, 0)),
                    0.5 * (coords(0, 1) + coords(2, 1)), smesh.CoordinateSystem);
                write_coords(OUT,
                    0.25 * (coords(0, 0) + coords(1, 0) + coords(2, 0) + coords(3, 0)),
                    0.25 * (coords(0, 1) + coords(1, 1) + coords(2, 1) + coords(3, 1)),
                    smesh.CoordinateSystem);
                write_coords(OUT, 0.5 * (coords(1, 0) + coords(3, 0)),
                    0.5 * (coords(1, 1) + coords(3, 1)), smesh.CoordinateSystem);
                write_coords(OUT, coords(2, 0), coords(2, 1), smesh.CoordinateSystem);
                write_coords(OUT, 0.5 * (coords(2, 0) + coords(3, 0)),
                    0.5 * (coords(2, 1) + coords(3, 1)), smesh.CoordinateSystem);
                write_coords(OUT, coords(3, 0), coords(3, 1), smesh.CoordinateSystem);
            }
        OUT << "CELLS " << smesh.nx * smesh.ny << " " << 10 * smesh.nx * smesh.ny << std::endl;
        size_t ii = 0;
        for (size_t iy = 0; iy < smesh.ny; ++iy)
            for (size_t ix = 0; ix < smesh.nx; ++ix, ++ii)

                OUT << "9 " << 9 * ii << " " << 9 * ii + 2 << " " << 9 * ii + 8 << " " << 9 * ii + 6
                    << " " << 9 * ii + 1 << " " << 9 * ii + 5 << " " << 9 * ii + 7 << " "
                    << 9 * ii + 3 << " " << 9 * ii + 4 << std::endl;

        OUT << "CELL_TYPES " << smesh.nx * smesh.ny << std::endl;
        for (ii = 0; ii < smesh.nx * smesh.ny; ++ii)
            OUT << "28 ";
        OUT << std::endl;

        OUT << "POINT_DATA " << 9 * smesh.nx * smesh.ny << std::endl
            << "SCALARS " << variableName << " DOUBLE " << std::endl
            << "LOOKUP_TABLE default" << std::endl;

        ii = 0;
        for (size_t iy = 0; iy < smesh.ny; ++iy)
            for (size_t ix = 0; ix < smesh.nx; ++ix, ++ii) {

                const Eigen::Matrix<double, 1, 9> vt = v.row(ii) * PSILagrange<6, 3>;
                for (int i = 0; i < 9; ++i)
                    OUT << cuttol(vt(i), 1.e-20) << std::endl;
            }
        OUT.close();
    }

    template <int DG>
    static void write_dg(
        const std::string& fname, int n, const DGVector<DG>& v, const ParametricMesh& smesh)
    {
        write_dg(compose_vtkname(fname, DG, n), v, smesh);
    }
    template <int CGdegree>
    static void write_cg(
        const std::string& fname, int n, const CGVector<CGdegree>& v, const ParametricMesh& smesh)
    {
        write_cgvector(compose_vtkname(fname, CGdegree, n), v, smesh);
    }
    template <int CGdegree>
    static void write_cg_velocity(const std::string& fname, int n, const CGVector<CGdegree>& vx,
        const CGVector<CGdegree>& vy, const ParametricMesh& smesh)
    {
        write_cg_velocity(compose_vtkname(fname, CGdegree, n), vx, vy, smesh);
    }
};

} /* namespace Nextsim */

#endif /* __DGVISU_HPP */
