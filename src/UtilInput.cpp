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
#include <regex>
#include <ctime>

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


int WhichDateFormat(const std::string& str)
{
    std::regex dateRegex("\\d{4}-\\d{2}-\\d{2}");
    if (std::regex_match(str, dateRegex))
        return 0;
    else
    {
        try
        {
            double num = std::stod(str);
            return 1;
        }
        catch (...)
        {
            OCP_ABORT("Wrong TIME format!");
        }
    }
}

int NumDaysBetweenDates(const string& begin, const string& end)
{
    std::tm tm1 = {};
    std::stringstream ss1(begin);
    ss1 >> std::get_time(&tm1, "%Y-%m-%d");

    std::tm tm2 = {};
    std::stringstream ss2(end);
    ss2 >> std::get_time(&tm2, "%Y-%m-%d");

    std::time_t time1 = std::mktime(&tm1);
    std::time_t time2 = std::mktime(&tm2);

    double diffSeconds = std::difftime(time2, time1);
    return diffSeconds / (60 * 60 * 24);
}



/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/