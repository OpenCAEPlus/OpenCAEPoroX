/*! \file    OCPTable.hpp
 *  \brief   OCPTable class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPTable_HEADER__
#define __OCPTable_HEADER__

// Standard header files
#include <iostream>
#include <vector>

// OpenCAEPoroX header files
#include "OCPConst.hpp"

using namespace std;

class OCPTable_Interpolate
{
public:
    OCPTable_Interpolate() = default;
    virtual void Interpolate(const vector<vector<OCP_DBL>>& table) = 0;
};


/// OCPTable is a Table class, which used to deal with everything about table
/// in OpenCAEPoroX such as PVT table, saturation table.
class OCPTable
{
public:
    /// Default constructor.
    OCPTable() = default;

    /// Construct from existing data
    OCPTable(const vector<vector<OCP_DBL>>& src);

    /// Setup tables from existing data of table.
    void Setup(const vector<vector<OCP_DBL>>& src);

    /// judge if table is empty.
    OCP_BOOL IsEmpty() const { return data.empty(); }

    /// return the column num of table.
    USI GetColNum() const { return nCol; }

    /// return the jth column in table to modify or use.
    vector<OCP_DBL>& GetCol(const USI& j) { return data[j]; }
    const vector<OCP_DBL>& GetCol(const USI& j) const { return data[j]; }

    /// interpolate the specified monotonically increasing column in table to evaluate
    /// all columns and return slope
    USI Eval_All(const USI&       j,
                 const OCP_DBL&   val,
                 vector<OCP_DBL>& outdata,
                 vector<OCP_DBL>& slope) const;
    USI Eval_All(const USI& j,
                 const OCP_DBL& val,
                 vector<OCP_DBL>& outdata) const;

    /// interpolate the specified monotonically increasing column in table to evaluate
    /// all columns, j = 0 here and index of returning date begins from 1
    USI Eval_All0(const OCP_DBL& val, vector<OCP_DBL>& outdata) const;

    /// interpolate the specified monotonically increasing column in table to evaluate
    /// the target column.
    OCP_DBL Eval(const USI& j, const OCP_DBL& val, const USI& destj) const;

    /// interpolate the specified monotonically increasing column in table to evaluate
    /// the target column, and return corresponding slope.
    OCP_DBL Eval(const USI& j, const OCP_DBL& val, const USI& destj, OCP_DBL& myK) const;

    /// interpolate the specified monotonically decreasing column in table to evaluate
    /// the target column.

    OCP_DBL Eval_Inv(const USI& j, const OCP_DBL& val, const USI& destj) const;

    /// Return the closest row away from specific val with given column j
    void GetCloseRow(const USI& j, const OCP_DBL& val, vector<OCP_DBL>& outdata) const;

    /// Display the data of table on screen.
    void Display() const;

protected:
    USI                     nRow; ///< number of rows of the table
    USI                     nCol; ///< number of columns of the table
    mutable USI             bId;  ///< the starting point of rows when interpolating
    vector<vector<OCP_DBL>> data; ///< data of the table, data[i] is the ith column.
};

#endif /* end if __OCPTABLE_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/