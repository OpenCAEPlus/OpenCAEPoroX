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


enum class BulkContent : USI
{
    /// all rock
    r,
    /// rock and fluid
    rf
};


class BulkVarSet
{

    /////////////////////////////////////////////////////////////////////
    // General Information
    /////////////////////////////////////////////////////////////////////

public:
    /// Number of bulks
    OCP_USI nb;  
    /// Number of bulks inside
    OCP_USI nbI; 
    /// Number of phase
    USI     np;
    /// Number of component
    USI     nc;          

	/////////////////////////////////////////////////////////////////////
	// Basic Grid and Basic Rock Information (numBulk)
	/////////////////////////////////////////////////////////////////////

public:
    /// Size of cell in x-direction
    vector<OCP_DBL>    dx;
    /// Size of cell in y-direction
    vector<OCP_DBL>    dy;
    /// Size of cell in z-direction
    vector<OCP_DBL>    dz;
    /// Volume of bulk
    vector<OCP_DBL>    v;
    /// Depth of center of bulk
    vector<OCP_DBL>    depth;

    /// net to gross of bulk
    vector<OCP_DBL>    ntg;
    /// initial rock porosity(*ntg)
    vector<OCP_DBL>    poroInit;
    /// rock porosity(* ntg)
    vector<OCP_DBL>    poro;
    /// pore volume = Vgrid * poro
    vector<OCP_DBL>    rockVp;
    /// rock permeability along the x direction
    vector<OCP_DBL>    rockKx;
    /// rock permeability along the y direction
    vector<OCP_DBL>    rockKy;
    /// rock permeability along the z direction
    vector<OCP_DBL>    rockKz;
    /// sigma factor used in dual porosity matrix-fracture coupling term
    vector<OCP_DBL>    sigma;
    /// vertical dimension of a block of matrix material
    vector<OCP_DBL>    dzMtrx;
    /// Volume of rock
    vector<OCP_DBL>    vr;
    /// Enthalpy of rock
    vector<OCP_DBL>    Hr;

    /// last poro
    vector<OCP_DBL>    lporo;   
    /// last rockVp    
    vector<OCP_DBL>    lrockVp;
    /// last vr        
    vector<OCP_DBL>    lvr;
    /// last Hr        
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

    /// last poroP
    vector<OCP_DBL>    lporoP;
    /// last poroT
    vector<OCP_DBL>    lporoT;
    /// last vrp
    vector<OCP_DBL>    lvrP;
    /// last vrT
    vector<OCP_DBL>    lvrT;
    /// last HrT
    vector<OCP_DBL>    lHrT;

    /////////////////////////////////////////////////////////////////////
    // Basic Fluid Information (numBulk)
    /////////////////////////////////////////////////////////////////////

public:

    /// Index of oil, gas, water(neigative = inexisting)
    INT                 o, g, w;
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
         
    /// last T
    vector<OCP_DBL>     lT;
    /// last P          
    vector<OCP_DBL>     lP;
    /// last vf         
    vector<OCP_DBL>     lvf;
    /// last Nt         
    vector<OCP_DBL>     lNt;
    /// last Ni         
    vector<OCP_DBL>     lNi; 
    /// last phaseExist
    vector<OCP_BOOL>    lphaseExist;
    /// last S
    vector<OCP_DBL>     lS;  
    /// last vj         
    vector<OCP_DBL>     lvj;
    /// last xij        
    vector<OCP_DBL>     lxij; 
    /// last xi         
    vector<OCP_DBL>     lxi; 
    /// last rho        
    vector<OCP_DBL>     lrho;        
    /// last mu
    vector<OCP_DBL>     lmu; 
    /// last kr         
    vector<OCP_DBL>     lkr;
    /// last Pc         
    vector<OCP_DBL>     lPc;
    /// last Pj         
    vector<OCP_DBL>     lPj; 
    /// last Uf          
    vector<OCP_DBL>     lUf;
    /// last H          
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
                        
    /// last vfP        
    vector<OCP_DBL>     lvfP;
    /// last vfT        
    vector<OCP_DBL>     lvfT;
    /// last vfi        
    vector<OCP_DBL>     lvfi;
    /// last xiP        
    vector<OCP_DBL>     lxiP;
    /// last xiT        
    vector<OCP_DBL>     lxiT;
    /// last xix        
    vector<OCP_DBL>     lxix;
    /// last rhoP       
    vector<OCP_DBL>     lrhoP;
    /// last rhoT       
    vector<OCP_DBL>     lrhoT; 
    /// last rhox       
    vector<OCP_DBL>     lrhox;    
    /// last muP        
    vector<OCP_DBL>     lmuP; 
    /// last muT        
    vector<OCP_DBL>     lmuT; 
    /// last mux        
    vector<OCP_DBL>     lmux; 
    /// last dPcdS      
    vector<OCP_DBL>     ldPcdS;
    /// last dKrdS      
    vector<OCP_DBL>     ldKrdS; 
    /// last UfP        
    vector<OCP_DBL>     lUfP;
    /// last UfT        
    vector<OCP_DBL>     lUfT;
    /// last Ufi        
    vector<OCP_DBL>     lUfi;
    /// last HT      
    vector<OCP_DBL>     lHT; 
    /// last Hx         
    vector<OCP_DBL>     lHx; 

    /// length of dSec_dPri.
    USI                 lendSdP;
    /// d Secondary variable / d Primary variable.
    vector<OCP_DBL>     dSec_dPri;  
    /// last dSec_dPri
    vector<OCP_DBL>     ldSec_dPri;    
   
    /// content type of bulk, all rock, rock and fluid
    vector<BulkContent> cType;
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
