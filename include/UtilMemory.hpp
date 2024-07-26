/** \file    UtilMemoery.hpp
 *  \brief   record the memory usuage
 *  \author  Shizhe Li
 *  \date    Jul/26/2024
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2024--present by the FASP++ team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __UTILMEMORY_HEADER__ /*-- allow multiple inclusions --*/
#define __UTILMEMORY_HEADER__ /**< indicate timing.hxx has been included --*/


#if defined(_CONSOLE) || defined(_WIN32) || defined(_WIN64)

#include <windows.h>
#include <psapi.h>
#include <iostream>

#else

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#endif


double OCPGetCurrentRSS();



#endif /*-- end if for __UTILMEMORY_HEADER__ --*/

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/26/2024      Create file                          */
/*----------------------------------------------------------------------------*/