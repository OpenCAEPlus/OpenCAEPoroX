/*! \file    PreParamGridWell.hpp
 *  \brief   PreParamGridWell class declaration
 *  \author  Shizhe Li
 *  \date    Feb/15/2023
 *
 *  \note    The params used in OpenCAEPoroX is mostly compatible with Eclipse by SLB,
 *           but it has some own rules for easy to use. It is extensible and friendly.
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __PARAMGRID_HEADER__
#define __PARAMGRID_HEADER__

// Standard header files
#include <fstream>
#include <iostream>
#include <string>

// OpenCAEPoroX header files
#include "OCPConst.hpp"
#include "UtilInput.hpp"
#include "ParamWell.hpp"
#include "CornerGrid.hpp"
#include "GmshGrid.hpp"
#include "Output4Vtk.hpp"

using namespace std;


/// Active cell indicator and its index among active cells.
//  Note: GB_Pair contains two variables, which indicates if a grid cell is active or
//  not and its index among active cells if it is active.
class GB_Pair
{
public:
    /// Default constructor.
    GB_Pair() = default;
    /// Constructor with given information.
    GB_Pair(const OCP_BOOL& act, const OCP_ULL& i)
        : activity(act)
        , index(i) {};

    /// Return whether this cell is active or not.
    auto IsAct() const { return activity; }
    /// Return the index of this cell among active cells.
    auto GetId() const { return index; }

private:
    OCP_BOOL activity; ///< Activeness of a grid cell.
    OCP_ULL  index;    ///< Active index of grid if active
};


/// Effective area of intersection surfaces with neighboring cells.
class ConnPair
{
public:
    ConnPair() = default;
    ConnPair(const OCP_ULL& Id,
        const OCP_INT&      Wgt,
        const ConnDirect&   direct,
        const OCP_DBL&      AreaB,
        const OCP_DBL&      AreaE)
        : id(Id)
        , wgt(Wgt)
        , direction(direct)
        , areaB(AreaB)
        , areaE(AreaE) {};
    const auto& ID() const { return id; }
    const auto& WGT() const { return wgt; }
    const auto& Direct() const { return direction; }
    const auto& AreaB() const { return areaB; }
    const auto& AreaE() const { return areaE; }
    void SetTransMult(const OCP_DBL& var) { transMult *= var; }
    const auto& TransMult() const { return transMult; }

protected:
    /// index of a neighboring cell
    OCP_ULL    id;
    /// weight of edge
    OCP_INT    wgt;   
    /// direction of connection
    ConnDirect direction;
    /// Effective intersection area between this cell and the neighbor, self
    OCP_DBL    areaB; 
    /// Effective intersection area between this cell and the neighbor, neighbor
    OCP_DBL    areaE; 
    /// Transmissibility multipliers,
    OCP_DBL    transMult{ 1.0 };
};


/// Initial reservoir contions for each grid, for example, saturation, pressure.
class InitialReservoir
{
    friend class PreParamGridWell;
    friend class Reservoir;

public:
    /// check and fill in the missing data if necessary
    void CheckData(const OCP_USI& numGrid);

protected:
    /// Initial water saturation
    vector<OCP_DBL> swat;
    /// Initial water saturation in each grid and use saturation end point scaling
    vector<OCP_DBL> swatInit;
    /// if Pcow should be scaled.
    OCP_BOOL        scalePcow{ OCP_FALSE };
};


/// Input grid information and well geometry information from input file
/// Generate the connections between active grids and wells
/// Output connections to file
class PreParamGridWell
{

    friend class Partition;
    friend class Domain;
    friend class Reservoir;

    /////////////////////////////////////////////////////////////////////
    // Input params from input files
    /////////////////////////////////////////////////////////////////////

public:
    /// Input params files
    void InputFile(const string& myFile, const string& myWorkdir);

protected:
    void Input(const string& myFilename);
    /// Check Input
    void CheckInput();
    /// Post-process input param 
    void PostProcessInput();
    /// Input the problem model: isothermal or thermal
    void InputMODEL(ifstream& ifs);
    /// Input DUALPORO
    void InputDUALPORO() { DUALPORO = OCP_TRUE; }
    /// Input DPGRID
    void InputDPGRID() { DPGRID = OCP_TRUE; }
    /// Input dimensions
    void InputDIMENS(ifstream& ifs);
    /// Input EQUALS
    void InputEQUALS(ifstream& ifs);
    /// Input COPY
    void InputCOPY(ifstream& ifs);
    /// Input MULTIPLY
    void InputMULTIPLY(ifstream& ifs);
    /// Input grid property
    void InputGrid(ifstream& ifs, string& keyword);
    /// Input WELSPECS
    void InputWELSPECS(ifstream& ifs);
    /// Input COMPDAT
    void InputCOMPDAT(ifstream& ifs);
    /// Input INCLUDE
    void InputINCLUDE(ifstream& ifs);

#ifdef WITH_GMSH
    /// Input GMSH
    void InputGMSH(ifstream& ifs);
    /// Input Physical Property
    void InputGMSHPRO(ifstream& ifs);
#endif
    // Input tools
    /// Find pointer to the specified variable.
    vector<OCP_DBL>* FindPtr(const string& varName, const OCP_DBL&);
    vector<USI>* FindPtr(const string& varName, const USI&);
    /// It's used in InputEQUALS, assigning values in batches.
    template <typename T>
    void setVal(vector<T>& obj, const T& val, const vector<USI>& index);
    /// It's used in InputCOPY, copying the value of one variable to another.
    template <typename T>
    void CopyVal(vector<T>& obj, const vector<T>& src, const vector<USI>& index);
    /// It's used in InputMULTIPLY, multipling the value of a certain range of a
    /// variable by a coefficient.
    void MultiplyVal(vector<OCP_DBL>& obj, const OCP_DBL& val, const vector<USI>& index);


protected:

    /// current work directory
    string          workdir;
    /// Orthogonal or Corner grid              
    GridType        gridType{ GridType::orthogonal };
    /// Num of grids(numGridM + numGridF)
    OCP_ULL         numGrid{ 0 };
    /// Num of matrix grid
    OCP_ULL         numGridM{ 0 };
    /// Num of fracture grid
    OCP_ULL         numGridF{ 0 };
    /// thermal model or isothermal model
    OCPModel        model{ OCPModel::none };

    // Orthogonal grid
    /// Num of grids along x-direction.
    USI             nx;
    /// Num of grids along y-direction.
    USI             ny;
    /// Num of grids along z-direction. 
    USI             nz; 
    /// Size along the x - direction for each grid.
    vector<OCP_DBL> dx;   
    /// Size along the y - direction for each grid.
    vector<OCP_DBL> dy;
    /// Size along the z - direction for each grid.
    vector<OCP_DBL> dz; 
    /// Depth of the top surface of the uppermost grids.
    vector<OCP_DBL> tops;  

    // Corner point grid
    /// Lines of a corner-point grid.   
    vector<OCP_DBL> coord; 
    /// Z-values of a corner-point grid.
    vector<OCP_DBL> zcorn; 

    // General grid
#ifdef WITH_GMSH
    /// GMSH Grid
    GMSHGrid        gmshGrid;
#endif

    // Rock param
    /// Net to gross
    vector<OCP_DBL> ntg; 
    /// Porosity
    vector<OCP_DBL> poro; 
    /// permeability for x-direction
    vector<OCP_DBL> kx;
    /// permeability for y-direction
    vector<OCP_DBL> ky;
    /// permeability for z-direction
    vector<OCP_DBL> kz;
    /// transmissibility multipliers in Z-direction
    vector<OCP_DBL> multZ;

    /// Activity of grid from input file: numGridLocal: 0 = inactive, 1 = active.
    vector<USI>     ACTNUM; 
    /// Records the index of SAT region for each grid.  
    vector<USI>     SATNUM;  
    /// Records the index of PVT region for each grid.
    vector<USI>     PVTNUM;
    /// Records the index of ROCK region for each grid.
    vector<USI>     ROCKNUM;   
    /// boundary indicator
    vector<USI>     boundIndex;
    /// boundary area
    vector<OCP_DBL> boundArea;

    // Initial Condition
    InitialReservoir initR;

    // Well
    vector<WellParam> well;


    /////////////////////////////////////////////////////////////////////
    // Generate grid connections
    /////////////////////////////////////////////////////////////////////

public:
    void Setup();

protected:
    /// Setup grids
    void SetupGrid();
    /// Setup orthogonal grid.
    void SetupOrthogonalGrid();
    /// Calculate the depth and volume for orthogonal grid.
    void CalDepthVOrthogonalGrid();
    /// Setup the neighboring info for an orthogonal grid.
    void SetupActiveConnOrthogonalGrid();
    /// Setup connections for single media in orthogonal grid.
    void SetupActiveConnOrthogonalGridSM();
    /// Setup connections for dual pororsity in orthogonal grid.
    void SetupActiveConnOrthogonalGridDP();
    /// Output grid points for orthogonal grid.
    void OutputPointsOrthogonalGrid();
    

    /// Setup corner-point grid.
    void SetupCornerGrid();
    /// Setup dx,dy,dz,depth, v for a corner-point grid.
    void SetupBasicCornerGrid(const OCP_COORD& CoTmp);
    /// Setup the neighboring info for a corner-point grid.
    void SetupActiveConnCornerGrid(const OCP_COORD& CoTmp);
    /// Setup connections for single media in corner-point grid.
    void SetupActiveConnCornerGridSM(const OCP_COORD& CoTmp);
    /// Setup connections for dual pororsity in corner-point grid.
    void SetupActiveConnCornerGridDP(const OCP_COORD& CoTmp);
    /// Output grid points for corner-point grid.
    void OutputPointsCornerGrid(const OCP_COORD& mycord);

#ifdef WITH_GMSH
    /// Setup gmsh grid.
    void SetupGmshGrid();
    /// Setup depth, v for a gmsh grid
    void SetupBasicGmshGrid();
    /// Setup the neighboring info for a gmsh grid
    void SetupActiveConnGmshGrid();
    /// Output grid points for a gmsh grid
    void OutputPointsGmshGrid();
#endif

    // For Structral Grid
    /// Set location for grid: top, bottom, side or interior
    void SetLocationStructral();
    

    /// For General Grid
    /// Calculate the activity of grid cells
    void CalActiveGrid(const OCP_DBL& e1, const OCP_DBL& e2);
    /// Calculate the activity of grid cells for isothermal model
    void CalActiveGridIsoT(const OCP_DBL& e1, const OCP_DBL& e2);
    /// Calculate the activity of grid cells for Thermal model
    void CalActiveGridT(const OCP_DBL& e1, const OCP_DBL& e2);
    /// Setup Transmissibility multipliers
    void SetupTransMult();


protected:
    /// Volume of cells
    vector<OCP_DBL>          v; 
    /// Depth of center of grid cells
    vector<OCP_DBL>          depth;

    // Connections
    /// Neighboring information of active grid.
    vector<vector<ConnPair>> gNeighbor; 
    /// Num of neighbor
    vector<USI>              numNeighbor; 

    /// Num of active grid.
    OCP_ULL                  activeGridNum;
    /// Index mapping from active grid to all grid
    vector<OCP_ULL>          map_Act2All;
    /// Mapping from grid to active all grid
    vector<GB_Pair>          map_All2Act; 
    /// Num of fluid grids.
    OCP_ULL                  fluidGridNum;
    /// Mapping from all grid to fluid grid
    vector<GB_Pair>          map_All2Flu;


    /////////////////////////////////////////////////////////////////////
    // Dual Porosity Option
    /////////////////////////////////////////////////////////////////////

protected:
    /// If use dual porosity option. (matrix is in the front and followed by fracture grid)
    OCP_BOOL        DUALPORO{ OCP_FALSE };
    /// if true, then the default property of fracture will be copied from matrix(think of it)
    OCP_BOOL        DPGRID{ OCP_FALSE };
    /// sigma factor used in dual porosity matrix-fracture coupling term
    vector<OCP_DBL> sigma;
    /// vertical dimension of a block of matrix material
    vector<OCP_DBL> dzMtrx;

    /////////////////////////////////////////////////////////////////////
    // Generate connections between active grids and wells
    /////////////////////////////////////////////////////////////////////

protected:
    /// Setup connections between wells and active grids
    void SetupConnWellGrid();
    /// Return the index of bulk perforated by well
    OCP_ULL GetPerfLocation(const WellParam& well, const USI& p);

protected:
    /// Num of wells
    USI                     numWell;  
    /// Connections between wells and active grids
    vector<vector<OCP_ULL>> connWellGrid;


    /////////////////////////////////////////////////////////////////////
    // Weight of connections
    /////////////////////////////////////////////////////////////////////

protected:
	const OCP_INT WEIGHT_GG = 1;         ///< grid-grid
	const OCP_INT WEIGHT_GW = 1000000;   ///< grid-well

    /////////////////////////////////////////////////////////////////////
    // Output basic grid information and active grid connections
    /////////////////////////////////////////////////////////////////////

public:
    void OutputBaiscInfo() const;   ///< Calculate and return basic informations for grid



protected:
    OCP_BOOL ifUseVtk{ OCP_FALSE };

public:
    void FreeMemory();
};

#endif



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Feb/15/2023      Create file                          */
/*----------------------------------------------------------------------------*/