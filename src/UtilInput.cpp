/*! \file    UtilInput.cpp
 *  \brief   Utilities for reading input file
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */
#include "UtilInput.hpp"

OCP_BOOL ReadLine(ifstream& ifs, vector<string>& result, const OCP_BOOL& no_slash)
{
    result.resize(0);
    string buf;

    while (!ifs.eof()) {
        getline(ifs, buf);
        if (buf.empty()) continue;
        while (buf[0] == ' ' || buf[0] == '\t' || buf[0] == '\r') buf.erase(0, 1);
        if (buf.empty() || buf[0] == '#') continue;
        if (buf.size() > 1 && (buf[0] == '-' && buf[1] == '-')) continue;

        break;
    }

    // file ends
    if (buf.empty()) return OCP_FALSE;


    // remove the string behind the '/'
    if (no_slash) {
        auto pos = buf.find_first_of('/');
        if (pos != string::npos) {
            buf.erase(pos);
            buf.push_back('/');
        }
    }

    // get rid of  '  and  ,
    for (auto& s : buf) {
        if (s == '\'' || s == ',') s = ' ';
    }

    istringstream tmp(buf);
    while (tmp >> buf) result.push_back(buf);

    return OCP_TRUE;
}


void DealDefault(vector<string>& result)
{
    vector<string> tmp;
    for (auto& str : result) {
        auto pos = str.find('*');
        if (pos == string::npos) {
            tmp.push_back(str);
        } else {
            USI    num = atoi(str.substr(0, pos).c_str()); // non number -> 0
            USI    len = str.size();
            string val = "DEFAULT";
            if (num == 0) {
                tmp.push_back(str);
            } else {
                if (pos != len - 1) {
                    val = str.substr(pos + 1, len - (pos + 1));
                }
                for (USI i = 0; i < num; i++) tmp.push_back(val);
            }
        }
    }
    swap(result, tmp);
}

vector<vector<string>> ExpandWellOptions(const vector<string>& result)
{
    vector<vector<string>> tmp;
    vector<string> new_ = result;

    string tstep_str = result[1];
    size_t pos = tstep_str.find('*');
    int num;
    if (pos == string::npos)
    {
        tmp.push_back(new_);
        return tmp;
    }
    else
    {
        num = stoi(tstep_str.substr(0, pos));
        string tstep = tstep_str.substr(pos+1, tstep_str.size() - (pos+1));
        new_[1] = tstep;
    }

    for (int i=0; i<num; ++i)
        tmp.push_back(new_);

    return tmp;
}



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/