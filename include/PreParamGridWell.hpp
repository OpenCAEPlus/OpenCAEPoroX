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
#include "CornerGrid.hpp"

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
    GB_Pair(const OCP_BOOL& act, const OCP_USI& i)
        : activity(act)
        , index(i) {};

    /// Return whether this cell is active or not.
    OCP_BOOL IsAct() const { return activity; }
    /// Return the index of this cell among active cells.
    OCP_USI GetId() const { return index; }

private:
    OCP_BOOL activity; ///< Activeness of a grid cell.
    OCP_USI  index;    ///< Active index of grid if active
};


/// Effective area of intersection surfaces with neighboring cells.
class GPair
{
public:
    GPair() = default;
    GPair(const OCP_USI& Id,
        const OCP_INT& Wgt,
        const USI& direct,
        const OCP_DBL& AreaB,
        const OCP_DBL& AreaE)
        : id(Id)
        , wgt(Wgt)
        , areaB(AreaB)
        , areaE(AreaE) 
    {
        switch (direct)
        {
        case 0:
            direction = ConnDirect::x;
            break;
        case 1:
            direction = ConnDirect::y;
            break;
        case 2:
            direction = ConnDirect::z;
            break;
        default:
            OCP_ABORT("WRONG CONNECTION DIRECTION!");
            break;
        }
    };
    const auto& ID() const { return id; }
    const auto& WGT() const { return wgt; }
    const auto& Direct() const { return direction; }
    const auto& AreaB() const { return areaB; }
    const auto& AreaE() const { return areaE; }

protected:
    /// index of a neighboring cell
    OCP_USI    id;
    /// weight of edge
    OCP_INT    wgt;   
    /// direction of connection
    ConnDirect direction;
    /// Effective intersection area between this cell and the neighbor, self
    OCP_DBL    areaB; 
    /// Effective intersection area between this cell and the neighbor, neighbor
    OCP_DBL    areaE; 
};


class PreParamWell
{
   
public:
    PreParamWell(vector<string>& info) {
        name = info[0];
        if (info[1] != "DEFAULT") group = info[1];
        I = stoi(info[2]);
        J = stoi(info[3]);
        if (info[4] != "DEFAULT") depth = stod(info[4]);
    }
    // static infomation
    // WELSPECS
    string  name;             ///< Name of Well
    string  group{ "FEILD" }; ///< Group the well belongs to.
    USI     I;                ///< I index of well.
    USI     J;                ///< J index of well.
    OCP_DBL depth{ -1.0 };    ///< Depth of well.

    // COMPDAT ---- for all perforation.
    vector<USI>     I_perf;   ///< I-index of perforation in grid.
    vector<USI>     J_perf;   ///< J-index of perforation in grid.
    vector<USI>     K_perf;   ///< K-index of perforation in grid.
    vector<OCP_DBL> WI;       ///< Transmissibility connection factor.
    vector<OCP_DBL> diameter; ///< Diameter of perforations.
    vector<OCP_DBL> kh;
    vector<OCP_DBL> skinFactor; ///< Skin factor.
    vector<string>  direction;  ///< Direction of perforations.
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
    void InputFile(const string& myFile, const string& myWorkdir);

protected:
    void Input(const string& myFilename);
    void CheckInput();
    /// Input the problem model: isothermal or thermal
    void InputMODEL(ifstream& ifs);
    /// Input DUALPORO
    void InputDUALPORO() { DUALPORO = OCP_TRUE; }
    /// Input DPGRID
    void InputDPGRID() { DPGRID = OCP_TRUE; }
    void InputDIMENS(ifstream& ifs);
    void InputEQUALS(ifstream& ifs);
    void InputCOPY(ifstream& ifs);
    void InputMULTIPLY(ifstream& ifs);
    /// Input grid property
    void InputGrid(ifstream& ifs, string& keyword);
    void InputWELSPECS(ifstream& ifs);
    void InputCOMPDAT(ifstream& ifs);
    void InputINCLUDE(ifstream& ifs);
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
    USI             gridType{ ORTHOGONAL_GRID };
    /// Num of grids.
    OCP_USI         numGrid;
    /// thermal model or isothermal model
    USI             model{ 0 };
    /// If use dual porosity option. (matrix is in the front and followed by fracture grid)
    OCP_BOOL        DUALPORO{ OCP_FALSE };
    /// if true, then the default property of fracture will be copied from matrix
    OCP_BOOL        DPGRID{ OCP_FALSE };

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

    // Region
    /// Activity of grid from input file: numGridLocal: 0 = inactive, 1 = active.
    vector<USI>     ACTNUM; 
    /// Records the index of SAT region for each grid.  
    vector<USI>     SATNUM;  
    /// Records the index of PVT region for each grid.
    vector<USI>     PVTNUM;
    /// Records the index of ROCK region for each grid.
    vector<USI>     ROCKNUM;   

    // Initial Condition
    /// Initial water saturation in each grid.
    vector<OCP_DBL> Swat; 
    /// if Pcow should be scaled.
    OCP_BOOL        scalePcow{ OCP_FALSE }; 


    /// Grid Location
    vector<USI>     location;


    // Well
    vector<PreParamWell> well;


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
    

    /// Setup corner-point grid.
    void SetupCornerGrid();
    /// Setup dx,dy,dz,depth, v for a corner-point grid.
    void SetupBasicCornerGrid(const OCP_COORD& CoTmp);
    /// Setup the neighboring info for a corner-point grid.
    void SetupActiveConnCornerGrid(const OCP_COORD& CoTmp);


    // For Structral Grid
    /// Set location for grid: top, bottom, side or interior
    void SetLocationStructral();
    

    /// Calculate the activity of grid cells
    void CalActiveGrid(const OCP_DBL& e1, const OCP_DBL& e2);
    /// Calculate the activity of grid cells for isothermal model
    void CalActiveGridIsoT(const OCP_DBL& e1, const OCP_DBL& e2);
    /// Calculate the activity of grid cells for Thermal model
    void CalActiveGridT(const OCP_DBL& e1, const OCP_DBL& e2);


protected:

    vector<OCP_DBL> v;     ///< Volume of cells: numGridLocal.
    vector<OCP_DBL> depth; ///< Depth of center of grid cells: numGridLocal.

    // Connections
    vector<vector<GPair>> gNeighbor;    ///< Neighboring information of active grid.
    vector<USI>           numNeighbor;  ///< Num of neighbor

    // Active grid cells
    OCP_USI activeGridNum; ///< Num of active grid.
    vector<OCP_USI>
        map_Act2All; ///< Mapping from active grid to all grid: activeGridNum.
    vector<GB_Pair> map_All2Act; ///< Mapping from grid to active all grid: numGridLocal.
    // Fluid grid cells
    OCP_USI         fluidGridNum; ///< Num of fluid grids.
    vector<GB_Pair> map_All2Flu;  ///< Mapping from all grid to fluid grid: numGridLocal.


    /////////////////////////////////////////////////////////////////////
    // Generate connections between active grids and wells
    /////////////////////////////////////////////////////////////////////

protected:
    /// Setup connections between wells and active grids
    void SetupConnWellGrid();

protected:

    USI                     numWell;      ///< Num of wells
    vector<vector<OCP_USI>> connWellGrid; ///< Connections between wells and active grids


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
    /// Output grid points for orthogonal grid.
    void OutputPointsOrthogonalGrid();
    /// Output grid points for corner-point grid.
    void OutputPointsCornerGrid(const OCP_COORD& mycord);

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