/*! \file    BulkConnArea.hpp
 *  \brief   BulkConnArea class declaration
 *  \author  Shizhe Li
 *  \date    Aug/24/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __BULKCONNAREA_HEADER__
#define __BULKCONNAREA_HEADER__


// OpenCAEPoroX header files
#include "BulkConnVarSet.hpp"
#include "Bulk.hpp"


#include <vector>

using namespace std;


class BulkConnAreaMethod
{
public:
    BulkConnAreaMethod() = default;
    virtual void CalArea(BulkConnPair& bp, const Bulk& bk) = 0;
};


class BulkConnAreaMethod01 : public BulkConnAreaMethod
{
public:
    BulkConnAreaMethod01() = default;
    void CalArea(BulkConnPair& bp, const Bulk& bk) override;
};



class BulkConnArea
{
public:
    BulkConnArea() = default;
    void Setup();
    void CalArea(BulkConnPair& bp, const Bulk& bk) { bcaM->CalArea(bp, bk); }

protected:
    /// method for area calculation of bulk connection
    BulkConnAreaMethod* bcaM;
};


#endif /* end if __BulkConnArea_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Aug/24/2023      Create file                          */
/*----------------------------------------------------------------------------*/
