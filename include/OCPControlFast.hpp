/*! \file    OCPControlFast.hpp
 *  \brief   OCPControlFast class declaration
 *  \author  Shizhe Li
 *  \date    Dec/05/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __OCPCONTROLFAST_HEADER__
#define __OCPCONTROLFAST_HEADER__

 // Standard header files
#include <vector>

// OpenCAEPoroX header files
#include "OCPConst.hpp"

using namespace std;


/// shortcut instructions from the command line
class FastControl
{
public:
    FastControl(const USI& argc, const char* optset[])
    {
        ifUse = OCP_FALSE;
        timeInit = timeMax = timeMin = -1.0;

        std::stringstream buffer;
        string            tmp;
        string            key;
        string            value;
        for (USI n = 3; n < argc; n++) {
            buffer << optset[n];
            buffer >> tmp;

            string::size_type pos = tmp.find_last_of('=');
            if (pos == string::npos) OCP_ABORT("Unknown Usage! See -h");

            key = tmp.substr(0, pos);
            value = tmp.substr(pos + 1, tmp.size() - pos);

            switch (Map_Str2Int(&key[0], key.size())) {

            case Map_Str2Int("method", 6):
                if (value == "FIM") {
                    method = OCPNLMethod::FIM;
                }
                else if (value == "IMPEC") {
                    method = OCPNLMethod::IMPEC;
                }
                else if (value == "AIMc") {
                    method = OCPNLMethod::AIMc;
                }
                else {
                    OCP_ABORT("Wrong method param in command line!");
                }
                ifUse = OCP_TRUE;
                if (method == OCPNLMethod::FIM || method == OCPNLMethod::AIMc) {
                    if (timeInit <= 0) timeInit = 1;
                    if (timeMax <= 0) timeMax = 10.0;
                    if (timeMin <= 0) timeMin = 0.1;
                }
                else {
                    if (timeInit <= 0) timeInit = 0.1;
                    if (timeMax <= 0) timeMax = 1.0;
                    if (timeMin <= 0) timeMin = 0.1;
                }
                break;

            case Map_Str2Int("dtInit", 6):
                timeInit = stod(value);
                break;

            case Map_Str2Int("dtMin", 5):
                timeMin = stod(value);
                break;

            case Map_Str2Int("dtMax", 5):
                timeMax = stod(value);
                break;

            case Map_Str2Int("verbose", 7):
                printLevel = OCP_MIN(OCP_MAX(stoi(value), PRINT_NONE), PRINT_ALL);
                break;

            default:
                OCP_ABORT("Unknown Options: " + key + "   See -h");
                break;
            }

            buffer.clear();
        }
    }


public:
    /// If use fastcontrol
    OCP_BOOL    ifUse{ OCP_FALSE };
    /// IMPEC, FIM or AIM
    OCPNLMethod method;
    /// length of the first time step beginning the next TSTEP
    OCP_DBL     timeInit;
    /// Maximum time step during running
    OCP_DBL     timeMax;
    /// Minimum time step during running
    OCP_DBL     timeMin;
    /// Print level
    USI         printLevel{ 0 };
};



#endif /* end if __OCPControlFast_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/05/2023      Create file                          */
/*----------------------------------------------------------------------------*/