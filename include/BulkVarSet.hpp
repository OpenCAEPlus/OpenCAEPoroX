/*! \file    BulkVarSet.hpp
 *  \brief   BulkVarSet class declaration
 *  \author  Shizhe Li
 *  \date    Aug/19/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __BULKVARSET_HEADER__
#define __BULKVARSET_HEADER__


// OpenCAEPoroX header files
#include "OCPConst.hpp"
#include <vector>

using namespace std;


/// Type of content within a bulk
enum class BulkContent : USI
{
    /// all rock
    r,
    /// rock and fluid
    rf
};


/// Basic variable set of bulk, which contains all needed variables for a basic simulation
class BulkVarSet
{

    /////////////////////////////////////////////////////////////////////
    // General Information
    /////////////////////////////////////////////////////////////////////

public:
    /// Number of bulks
    OCP_USI nb;  
    /// Number of interior bulks
    OCP_USI nbI; 
    /// Number of phases
    USI     np;
    /// Number of components
    USI     nc;          

	/////////////////////////////////////////////////////////////////////
	// Basic Grid and Basic Rock Information (numBulk)
	/////////////////////////////////////////////////////////////////////

public:
    /// Size of cell along x-direction
    vector<OCP_DBL>    dx;
    /// Size of cell along y-direction
    vector<OCP_DBL>    dy;
    /// Size of cell along z-direction
    vector<OCP_DBL>    dz;
    /// Volume of bulk
    vector<OCP_DBL>    v;
    /// Depth of center of bulk
    vector<OCP_DBL>    depth;

    /// net to gross of bulk
    vector<OCP_DBL>    ntg;
    /// initial rock porosity(*ntg)
    vector<OCP_DBL>    poroInit;
    /// rock porosity(*ntg)
    vector<OCP_DBL>    poro;
    /// pore volume = Vgrid * poro
    vector<OCP_DBL>    rockVp;
    /// rock permeability along x-direction
    vector<OCP_DBL>    rockKx;
    /// rock permeability along x-direction
    vector<OCP_DBL>    rockKy;
    /// rock permeability along x-direction
    vector<OCP_DBL>    rockKz;
    /// sigma factor used in dual porosity matrix-fracture coupling term
    vector<OCP_DBL>    sigma;
    /// vertical dimension of a block of matrix material
    vector<OCP_DBL>    dzMtrx;
    /// Volume of rock
    vector<OCP_DBL>    vr;
    /// Enthalpy of rock
    vector<OCP_DBL>    Hr;

    /// poro at the previous time step
    vector<OCP_DBL>    lporo;   
    /// rockVp at the previous time step     
    vector<OCP_DBL>    lrockVp;
    /// vr at the previous time step         
    vector<OCP_DBL>    lvr;
    /// Hr at the previous time step        
    vector<OCP_DBL>    lHr;    

    /// d poro / d P
    vector<OCP_DBL>    poroP;
    /// d poro / d T
    vector<OCP_DBL>    poroT;
    /// d vr / d p
    vector<OCP_DBL>    vrP;
    /// dvr / dT
    vector<OCP_DBL>    vrT;
    /// dHr / dT
    vector<OCP_DBL>    HrT;

    /// poroP at the previous time step
    vector<OCP_DBL>    lporoP;
    /// poroT at the previous time step
    vector<OCP_DBL>    lporoT;
    /// vrp at the previous time step
    vector<OCP_DBL>    lvrP;
    /// vrT at the previous time step
    vector<OCP_DBL>    lvrT;
    /// HrT at the previous time step
    vector<OCP_DBL>    lHrT;

    /////////////////////////////////////////////////////////////////////
    // Basic Fluid Information (numBulk)
    /////////////////////////////////////////////////////////////////////

public:

    /// Index of oil, gas, water(neigative = inexisting)
    INT                 o, g, w;
    /// Index of wetting phase, r is one of o,g,w -- for test now
    INT                 r;
    /// Temperature
    vector<OCP_DBL>     T;
    /// Pressure
    vector<OCP_DBL>     P;
    /// Bubble point pressure
    vector<OCP_DBL>     Pb;
    /// Total fluid volume
    vector<OCP_DBL>     vf;
    /// Total moles of components
    vector<OCP_DBL>     Nt;
    /// Moles of component
    vector<OCP_DBL>     Ni;
    /// Existence of phase
    vector<OCP_BOOL>    phaseExist;
    /// Saturation of phase
    vector<OCP_DBL>     S;
    /// Volume of phase
    vector<OCP_DBL>     vj;
    /// molar fraction of component i in phase j
    vector<OCP_DBL>     xij;
    /// Moles density of phase
    vector<OCP_DBL>     xi;
    /// Mass density of phase
    vector<OCP_DBL>     rho;
    /// Viscosity of phase
    vector<OCP_DBL>     mu;
    /// Relative permeability of phase
    vector<OCP_DBL>     kr;
    /// Capillary pressure of phase
    vector<OCP_DBL>     Pc;
    /// Pressure of phase
    vector<OCP_DBL>     Pj;
    /// Internal energy of fluid
    vector<OCP_DBL>     Uf;
    /// Enthalpy of phase
    vector<OCP_DBL>     H;
         
    /// T at the previous time step
    vector<OCP_DBL>     lT;
    /// P at the previous time step     
    vector<OCP_DBL>     lP;
    /// vf at the previous time step          
    vector<OCP_DBL>     lvf;
    /// Nt at the previous time step        
    vector<OCP_DBL>     lNt;
    /// Ni at the previous time step        
    vector<OCP_DBL>     lNi; 
    /// phaseExist at the previous time step
    vector<OCP_BOOL>    lphaseExist;
    /// S at the previous time step
    vector<OCP_DBL>     lS;  
    /// vj at the previous time step        
    vector<OCP_DBL>     lvj;
    /// xij at the previous time step       
    vector<OCP_DBL>     lxij; 
    /// xi at the previous time step        
    vector<OCP_DBL>     lxi; 
    /// rho at the previous time step       
    vector<OCP_DBL>     lrho;        
    /// mu at the previous time step
    vector<OCP_DBL>     lmu; 
    /// kr at the previous time step        
    vector<OCP_DBL>     lkr;
    /// Pc at the previous time step        
    vector<OCP_DBL>     lPc;
    /// Pj at the previous time step        
    vector<OCP_DBL>     lPj; 
    /// Uf at the previous time step         
    vector<OCP_DBL>     lUf;
    /// H at the previous time step         
    vector<OCP_DBL>     lH;       

    /// d vf / d P
    vector<OCP_DBL>     vfP;
    /// d vf / d T      
    vector<OCP_DBL>     vfT;
    /// d vf / d Ni     
    vector<OCP_DBL>     vfi;
    /// d xi / d P      
    vector<OCP_DBL>     xiP;
    /// d xi / d T      
    vector<OCP_DBL>     xiT;
    /// d Xi / d xij
    vector<OCP_DBL>     xix;
    /// d rho  / d P 
    vector<OCP_DBL>     rhoP;
    /// d rho / d T  
    vector<OCP_DBL>     rhoT;
    /// d rho / d xij
    vector<OCP_DBL>     rhox;    
    /// d mu / d P      
    vector<OCP_DBL>     muP;
    /// d mu / d T      
    vector<OCP_DBL>     muT;
    /// d mu / d xij
    vector<OCP_DBL>     mux;
    /// d Pc / d S      
    vector<OCP_DBL>     dPcdS;
    /// d Kr / d S      
    vector<OCP_DBL>     dKrdS;  
    /// d Uf / d P      
    vector<OCP_DBL>     UfP;
    /// d Uf / d T      
    vector<OCP_DBL>     UfT;
    /// d Uf / d Ni     
    vector<OCP_DBL>     Ufi;
    /// d H / d T       
    vector<OCP_DBL>     HT;
    /// d H / d xij     
    vector<OCP_DBL>     Hx;    
                        
    /// vfP at the previous time step         
    vector<OCP_DBL>     lvfP;
    /// vfT at the previous time step        
    vector<OCP_DBL>     lvfT;
    /// vfi at the previous time step        
    vector<OCP_DBL>     lvfi;
    /// xiP at the previous time step        
    vector<OCP_DBL>     lxiP;
    /// xiT at the previous time step        
    vector<OCP_DBL>     lxiT;
    /// xix at the previous time step        
    vector<OCP_DBL>     lxix;
    /// rhoP at the previous time step       
    vector<OCP_DBL>     lrhoP;
    /// rhoT at the previous time step       
    vector<OCP_DBL>     lrhoT; 
    /// rhox at the previous time step       
    vector<OCP_DBL>     lrhox;    
    /// muP at the previous time step        
    vector<OCP_DBL>     lmuP; 
    /// muT at the previous time step        
    vector<OCP_DBL>     lmuT; 
    /// mux at the previous time step        
    vector<OCP_DBL>     lmux; 
    /// dPcdS at the previous time step     
    vector<OCP_DBL>     ldPcdS;
    /// dKrdS at the previous time step       
    vector<OCP_DBL>     ldKrdS; 
    /// UfP at the previous time step        
    vector<OCP_DBL>     lUfP;
    /// UfT at the previous time step        
    vector<OCP_DBL>     lUfT;
    /// Ufi at the previous time step        
    vector<OCP_DBL>     lUfi;
    /// HT at the previous time step      
    vector<OCP_DBL>     lHT; 
    /// Hx at the previous time step         
    vector<OCP_DBL>     lHx; 

    /// length of dSec_dPri.
    USI                 lendSdP;
    /// d Secondary variable / d Primary variable.
    vector<OCP_DBL>     dSec_dPri;  
    /// dSec_dPri at the previous time step 
    vector<OCP_DBL>     ldSec_dPri;    
   
    /// content type of bulk, all rock, rock and fluid
    vector<BulkContent> cType;
    /// initialization type
    string  initType;
    /// Initial temperature
    vector<OCP_DBL>     initT;
};



#endif /* end if __BulkVarSet_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/19/2023      Create file                          */
/*----------------------------------------------------------------------------*/
