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
        const USI& Direct,
        const OCP_DBL& AreaB,
        const OCP_DBL& AreaE)
        : id(Id)
        , wgt(Wgt)
        , direction(Direct)
        , areaB(AreaB)
        , areaE(AreaE) {};
    static OCP_BOOL lessG(const GPair& G1, const GPair& G2) { return G1.id < G2.id; }

    OCP_USI id;        ///< Id of a neighboring cell
    OCP_INT wgt;       ///< weight of edge
    USI     direction; ///< direction: 1-x, 2-y, 3-z
    OCP_DBL
        areaB; ///< Effective intersection area between this cell and the neighbor, self
    OCP_DBL areaE; ///< Effective intersection area between this cell and the neighbor,
                   ///< neighbor
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
    void InputDIMENS(ifstream& ifs);
    void InputEQUALS(ifstream& ifs);
    void InputCOPY(ifstream& ifs);
    void InputMULTIPLY(ifstream& ifs);
    void InputBasic(ifstream& ifs, string& keyword);
    void InputRegion(ifstream& ifs, const string& keyword);
    void InputWELSPECS(ifstream& ifs);
    void InputCOMPDAT(ifstream& ifs);
    void InputINCLUDE(ifstream& ifs);
    // Input tools
    /// Find pointer to the specified variable.
    vector<OCP_DBL>* FindPtr(const string& varName);
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

    string  workdir;  ///< current work directory

    USI     gridType{ ORTHOGONAL_GRID }; ///< Orthogonal or Corner grid
    OCP_USI numGrid;  ///< Num of grids.

    // Orthogonal grid
    USI             nx;    ///< Num of grids along x-direction.
    USI             ny;    ///< Num of grids along y-direction.
    USI             nz;    ///< Num of grids along z-direction.  
    vector<OCP_DBL> dx;    ///< Size along the x - direction for each grid.
    vector<OCP_DBL> dy;    ///< Size along the y - direction for each grid.
    vector<OCP_DBL> dz;    ///< Size along the z - direction for each grid.
    vector<OCP_DBL> tops;  ///< Depth of the top surface of the uppermost grids.

    // Corner point grid
    vector<OCP_DBL> coord; ///< Lines of a corner-point grid.
    vector<OCP_DBL> zcorn; ///< ZValues of a corner-point grid.

    // Rock param
    vector<OCP_DBL> ntg;   ///< Net to gross for each grid.
    vector<OCP_DBL> poro;  ///< Porosity for each grid.
    USI             model{ 0 }; ///< if thermal model will be used.

    // Grid activity
    vector<USI> ACTNUM;    ///< Indicate activity of grid from input file: numGridLocal.
                           ///< 0 = inactive, 1 = active.

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
    /// Setup the grid information and calculate the properties.
    void SetupInitGrid();
    /// Setup orthogonal grid.
    void SetupOrthogonalGrid();
    /// Calculate the depth and volume for orthogonal grid.
    void CalDepthVOrthogonalGrid();
    /// Setup the neighboring info for an orthogonal grid.
    void SetupActiveConnOrthogonalGrid();
    /// Output grid points for orthogonal grid.
    void OutputPointsOrthogonalGrid();

    /// Setup corner-point grid.
    void SetupCornerGrid();
    /// Setup dx,dy,dz,depth, v for a corner-point grid.
    void SetupBasicCornerGrid(const OCP_COORD& CoTmp);
    /// Setup the neighboring info for a corner-point grid.
    void SetupActiveConnCornerGrid(const OCP_COORD& CoTmp);
    /// Output grid points for corner-point grid.
    void OutputPointsCornerGrid(const OCP_COORD& mycord);

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