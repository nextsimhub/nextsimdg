/*!
 * @file DynamicsKernel.hpp
 *
 * @date 17 Feb 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#ifndef DYNAMICSKERNEL_HPP
#define DYNAMICSKERNEL_HPP

#include "Interpolations.hpp"
#include "ParametricMesh.hpp"
#include "ParametricTools.hpp"
#include "DGTransport.hpp"
#include "cgParametricMomentum.hpp"

#include "cgVector.hpp"
#include "dgVector.hpp"
#include "dgLimit.hpp"
#include "dgVisu.hpp"
#include "Tools.hpp"


#include "CGModelArray.hpp"
#include "DGModelArray.hpp"
#include "include/gridNames.hpp"
#include "include/Time.hpp"


#include <string>
#include <unordered_map>

namespace Nextsim {

template <int CGdegree, int DGadvection> class DynamicsKernel {
public:
    
    void initialisation() {

        //! Define the spatial mesh
        smesh = new Nextsim::ParametricMesh(Nextsim::CARTESIAN);
        // FIXME integrate the creation of the smesh based on restart file
        //smesh->readmesh("init_topaz128x128.smesh"); // file temporary committed
        smesh->readmesh("25km_NH.smesh");
        
        // output land mask
        Nextsim::DGVector<1> landmask(*smesh);
        for (size_t i=0;i<smesh->nelements;++i)
            landmask(i,0) = smesh->landmask[i];
        Nextsim::VTK::write_dg<1>("landmask", 0, landmask, *smesh);


        // output boundary info
        Nextsim::DGVector<1> boundary(*smesh);
        for (size_t j=0;j<4;++j)
            for (size_t i=0;i<smesh->dirichlet[j].size();++i)
            boundary(smesh->dirichlet[j][i],0) = 1+j;
        Nextsim::VTK::write_dg<1>("boundary", 0, boundary, *smesh);


        //! Initialize transport
        dgtransport = new Nextsim::DGTransport<DGadvection>(*smesh);
        dgtransport->settimesteppingscheme("rk2");

        //! Initialize momentum
        momentum = new Nextsim::CGParametricMomentum<CGdegree>(*smesh);


        //! initialize Forcing 
        Nextsim::Interpolations::Function2CG(*smesh, momentum->GetOceanx(), OceanX());
        Nextsim::Interpolations::Function2CG(*smesh, momentum->GetOceany(), OceanY());

        AtmForcingX = new AtmX ;
        AtmForcingY = new AtmY ;

        AtmForcingX->settime(0.0);
        AtmForcingY->settime(0.0);
        Nextsim::Interpolations::Function2CG(*smesh, momentum->GetAtmx(), *AtmForcingX);
        Nextsim::Interpolations::Function2CG(*smesh, momentum->GetAtmy(), *AtmForcingY);


        //resize CG and DG vectors
        hice.resize_by_mesh(*smesh);
        cice.resize_by_mesh(*smesh);

        u.resize_by_mesh(*smesh);
        v.resize_by_mesh(*smesh);
    }


    /*
        Temrorary Ocen and Wind Field
    */

    class OceanX : public Nextsim::Interpolations::Function {
    public:
        double operator()(double x, double y) const
        {
            return  0.01 * (2.0 * y / ( 121. * 25000.) - 1.0);
        }
    };
    class OceanY : public Nextsim::Interpolations::Function {
    public:
        double operator()(double x, double y) const
        {
            return  0.01 * (1.0 - 2.0 * x / (154 * 25000.) );
        }
    };

    struct AtmX : public Nextsim::Interpolations::Function {
        double time;

    public:
        void settime(double t) { time = t; }
        double operator()(double x, double y) const
        {
            constexpr double oneday = 24.0 * 60.0 * 60.0;
            double vmax_atm = 30.0 / exp(1.0);
            //! Center of cyclone (in m)
            double cMx = 154./2. * 25000. + 154./10.* 25000. * time / oneday;
            double cMy = 121./2. * 25000. + 121./10.* 25000. * time / oneday;

            //! scaling factor to reduce wind away from center
            double scale = exp(1.0) / 100.0 * exp(-0.01e-3 * sqrt(SQR(x - cMx) + SQR(y - cMy))) * 1.e-3;

            double alpha = 72.0 / 180.0 * M_PI;

            return -scale * vmax_atm * 50 * (cos(alpha) * (x - cMx) + sin(alpha) * (y - cMy));
        }
    };
    struct AtmY : public Nextsim::Interpolations::Function {
        double time;

    public:
        void settime(double t) { time = t; }
        double operator()(double x, double y) const
        {
            constexpr double oneday = 24.0 * 60.0 * 60.0;
            double vmax_atm = 30.0 / exp(1.0);
            //! Center of cyclone (in m)
            double cMx = 154./2. * 25000. + 154./10.* 25000. * time / oneday;
            double cMy = 121./2. * 25000. + 121./10.* 25000. * time / oneday;
            //! scaling factor to reduce wind away from center
            double scale = exp(1.0) / 100.0 * exp(-0.01e-3 * sqrt(SQR(x - cMx) + SQR(y - cMy))) * 1.e-3;

            double alpha = 72.0 / 180.0 * M_PI;

            return -scale * vmax_atm * 50 * (-sin(alpha) * (x - cMx) + cos(alpha) * (y - cMy));
        }
    };

    
        
    /*!
     * @brief Sets the data from a provided ModelArray.
     *
     * @details Given a name and a ModelArray, sets the data associated with that
     * name. In some special cases (hice, cice…) this is a special array used in
     * the dynamics calculations. In all other cases, these are added to the
     * container of name data fields to be advected. The provided ModelArray can be
     * of DG or DGSTRESS type, in which case all components of the DGVector are
     * filled, or any other type which only fills the DG0 finite volume element of
     * the dgVector. The behaviour is exactly that of the ma2dg() function defined
     * in the DGModelArray class.
     *
     * @param name The name of the data field to set.
     * @param data The ModelArray containing the data to be set.
     *
     */
    void setData(const std::string& name, const ModelArray& data)
    {

        // Special cases: hice, cice, (damage, stress) <- not yet implemented
        if (name == hiceName) {
            DGModelArray::ma2dg(data, hice);
        } else if (name == ciceName) {
            DGModelArray::ma2dg(data, cice);
        } else if (name == "u") {
            // FIXME take into account possibility to restart form CG
            //CGModelArray::ma2cg(data, u);
            DGVector<DGadvection> utmp(*smesh);
            DGModelArray::ma2dg(data, utmp);
            Nextsim::Interpolations::DG2CG(*smesh, u, utmp);
        } else if (name == "v") {
            //CGModelArray::ma2cg(data, v);
            DGVector<DGadvection> vtmp(*smesh);
            DGModelArray::ma2dg(data, vtmp);
            Nextsim::Interpolations::DG2CG(*smesh, v, vtmp);
        } else {
            // All other fields get shoved in a (labelled) bucket
            DGModelArray::ma2dg(data, advectedFields[name]);
        }

    }



    /*!
     * @brief Returns an HField ModelArray containing the DG0 finite volume
     * component of the named dynamics field.
     *
     * @param name the name of the requested field.
     *
     */
    ModelArray getDG0Data(const std::string& name)
    {
        HField data(ModelArray::Type::H);
        if (name == hiceName) {
            return DGModelArray::dg2ma(hice, data);
        } else if (name == ciceName) {
            return DGModelArray::dg2ma(cice, data);
        } else if (name == uName) {
            DGVector<DGadvection> utmp(*smesh);
            Nextsim::Interpolations::CG2DG(*smesh, utmp, u);
            return DGModelArray::dg2ma(utmp, data);
        } else if (name == vName) {
            DGVector<DGadvection> vtmp(*smesh);
            Nextsim::Interpolations::CG2DG(*smesh, vtmp, v);
            return DGModelArray::dg2ma(vtmp, data);
        } else {
            // Any other named field must exist
            return DGModelArray::dg2ma(advectedFields.at(name), data);
        }
    }

    /*!
     * @brief Returns a DG or DGSTRESS ModelArray containing the full DG data for
     * the named dynamics field.
     *
     * @param name the name of the requested field.
     */
    ModelArray getDGData(const std::string& name)
    {

        if (name == hiceName) {
            DGField data(ModelArray::Type::DG);
            data.resize();
            return DGModelArray::dg2ma(hice, data);
        } else if (name == ciceName) {
            DGField data(ModelArray::Type::DG);
            data.resize();
            return DGModelArray::dg2ma(cice, data);
        } else {
            ModelArray::Type type = fieldType.at(name);
            ModelArray data(type);
            data.resize();
            return DGModelArray::dg2ma(advectedFields.at(name), data);
        }
    }
    
    void update(const TimestepTime& tst) {

        static int step_number=0;
        
        //! Tmp Set forcing
        AtmForcingX->settime(step_number*tst.step.seconds());
        AtmForcingY->settime(step_number*tst.step.seconds());
        Nextsim::Interpolations::Function2CG(*smesh, momentum->GetAtmx(), *AtmForcingX);
        Nextsim::Interpolations::Function2CG(*smesh, momentum->GetAtmy(), *AtmForcingY);
        

        //! interpolates CG velocity to DG and reinits normal velocity
        int vtk_out = 60; // 1h = 30 (timestep = 120s)
        std::string resultsdir = "vtk";
        if (step_number % vtk_out == 0) {
            std::cout << "Save vtk at " << step_number / 30  << " h " << std::endl;
            
            Nextsim::VTK::write_dg<6>(resultsdir + "/h", step_number / vtk_out  , hice, *smesh);
            Nextsim::VTK::write_dg<6>(resultsdir + "/c", step_number / vtk_out  , cice, *smesh);
    
    
            Nextsim::VTK::write_cg_velocity(resultsdir + "/vel", step_number / vtk_out, 
                momentum->GetVx(), momentum->GetVy(), *smesh);
            Nextsim::VTK::write_cg_velocity(resultsdir + "/atm", step_number / vtk_out, 
                momentum->GetAtmx(), momentum->GetAtmy(), *smesh);
    
            Nextsim::VTK::write_dg(resultsdir + "/Shear", step_number / vtk_out, 
            Nextsim::Tools::Shear(*smesh, momentum->GetE11(), momentum->GetE12(), momentum->GetE22()), *smesh);

        }

        //! Perform transport step
        dgtransport->prepareAdvection(momentum->GetVx(), momentum->GetVy());

        dgtransport->step(tst.step.seconds(), cice);	    
        dgtransport->step(tst.step.seconds(), hice);

        //! Gauss-point limiting
        Nextsim::LimitMax(cice, 1.0);
        Nextsim::LimitMin(cice, 0.0);
        Nextsim::LimitMin(hice, 0.0);

        
        momentum->prepareIteration(hice, cice);
        //! Momentum
        for (size_t mevpstep = 0; mevpstep < NT_evp; ++mevpstep) {
	        momentum->mEVPStep(VP, NT_evp, alpha, beta, tst.step.seconds(), hice, cice);
	    }

        step_number++;
        
    };

private:
    DGVector<DGadvection> hice;
    DGVector<DGadvection> cice;
    CGVector<CGdegree> u;
    CGVector<CGdegree> v;

    Nextsim::DGTransport<DGadvection>* dgtransport;
    Nextsim::CGParametricMomentum<CGdegree>* momentum;

    AtmX* AtmForcingX;
    AtmY* AtmForcingY;

    Nextsim::ParametricMesh* smesh;
    
    //! Rheology-Parameters
    Nextsim::VPParameters VP;
    //! MEVP parameters
    double alpha = 1500.0;
    double beta = 1500.0;
    size_t NT_evp = 100;

    std::unordered_map<std::string, DGVector<DGadvection>> advectedFields;

    // A map from field name to the type of
    const static std::unordered_map<std::string, ModelArray::Type> fieldType;
};

} /* namespace Nextsim */

#endif /* DYNAMICSKERNEL_HPP */
