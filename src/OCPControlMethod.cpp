/*! \file    ControlMethod.cpp
 *  \brief   ControlMethod class definition
 *  \author  Shizhe Li
 *  \date    Dec/05/2023
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPControlMethod.hpp"

void ControlMethod::SetCtrlParam(const ParamControl& CtrlParam)
{
	model = CtrlParam.model;
    for (const auto& m : CtrlParam.method) {
        if (m == "IMPEC") {
            method.push_back(OCPNLMethod::IMPEC);
        }
        else if (m == "FIM") {
            method.push_back(OCPNLMethod::FIM);
        }
        else if (m == "AIMc") {
            method.push_back(OCPNLMethod::AIMc);
        }
        else if (m == "FIMddm") {
            method.push_back(OCPNLMethod::FIMddm);
        }
        else {
            OCP_ABORT("Wrong method specified!");
        }
    }
    workDir = CtrlParam.workDir;
    lsFile  = CtrlParam.lsFile;

    if (method.size() == 0)  OCP_ABORT("METHOD is not input correctly!");
}


void ControlMethod::SetFastControl(const FastControl& fCtrl)
{
    method.clear();
    lsFile.clear();

    method.push_back(fCtrl.method);
    switch (method[0]) {
    case OCPNLMethod::IMPEC:
        lsFile.push_back("./csr.fasp");
        break;
    case OCPNLMethod::AIMc:
    case OCPNLMethod::FIM:
        lsFile.push_back("./bsr.fasp");
        break;
    default:
        OCP_ABORT("Wrong method specified from command line!");
        break;
    }
}


OCPNLMethod ControlMethod::InitMethod() const
{
    if (method.size() > 1) {
        wIndex = 1;
        return method[wIndex];
    }
    else {
        return method[0];
    }
}


OCPNLMethod ControlMethod::SwitchMethod() const
{
    if (method.size() > 1) {
        wIndex = (++wIndex) % (method.size());
        return method[wIndex];
    }
    else {
        return method[0];
    }
}


void ControlMethod::OutputInfo(const USI& pl) const
{
    if (CURRENT_RANK == MASTER_PROCESS) {
        switch (model) {
        case OCPModel::isothermal:
            if (pl >= PRINT_MIN) {
                cout << endl << "Dynamic simulation for isothermal models" << endl;
            }
            break;
        case OCPModel::thermal:
            if (pl >= PRINT_MIN) {
                cout << endl << "Dynamic simulation for thermal models" << endl;
            }
            break;
        default:
            OCP_ABORT("Wrong model type specified!");
        }

        cout << "\nDynamic simulation with ";
        cout << "main method ";
        switch (method[0])
        {
        case OCPNLMethod::FIM:
            cout << "FIM ";
            break;
        case OCPNLMethod::IMPEC:
            cout << "IMPEC ";
            break;
        case OCPNLMethod::AIMc:
            cout << "AIMc ";
            break;
        case OCPNLMethod::FIMddm:
            cout << "FIMddm ";
            break;
        default:
            break;
        }
        if (method.size() == 1) {
            cout << "!" << endl;
        }
        else {
            for (USI i = 1; i < method.size(); i++) {
                cout << "preconditioned by ";
                switch (method[i])
                {
                case OCPNLMethod::FIM:
                    cout << "FIM ";
                    break;
                case OCPNLMethod::IMPEC:
                    cout << "IMPEC ";
                    break;
                case OCPNLMethod::AIMc:
                    cout << "AIMc ";
                    break;
                case OCPNLMethod::FIMddm:
                    cout << "FIMddm ";
                    break;
                default:
                    break;
                }
            }
            cout << "!" << endl;
        }
    }
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Dec/05/2023      Create file                          */
/*----------------------------------------------------------------------------*/