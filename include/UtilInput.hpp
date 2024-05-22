/*! \file    UtilInput.hpp
 *  \brief   Supply basic tools used to input files.
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#ifndef __UTILINPUT_HEADER__
#define __UTILINPUT_HEADER__

// Standard header files
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>

// OpenCAEPoroX header files
#include "OCPConst.hpp"

using namespace std;

/// TODO: Replace it with error log
#define ParamCheck1(exp)                                                               \
    std::cout << exp << " in " << __func__ << "() in " << __LINE__ << " in "           \
              << __FILE__;

/// Map_str2int is used to map string to integer, which used to match the keyword
/// efficiently in input file in the switch structure.
constexpr inline long long Map_Str2Int(const char* mystr, const USI& len)
{
    long long res = 0;
    long long t   = 100;
    for (USI i = 0; i < len; i++) {
        res += (int)mystr[len - 1 - i] * t;
        t *= 100;
    }
    return res;
}


/// ReadLine is the core function while inputting the file. It will capture the next
/// line which is meanningful, for example, not blank line or comments, and then gets
/// rid of some useless characters such as space, commas. Finally, the segments of rest
/// string will be stored in result. And if return OCP_FALSE, it indicates we have reach
/// the end of file.
OCP_BOOL ReadLine(ifstream& ifs, vector<string>& result, const OCP_BOOL& no_slash = OCP_TRUE);

/// DealDefault is used to deal with the expression with asterisk, for example
/// m*n  -> <n,...,n> size m ,  m* -> <DEFAULT,..., DEFAULT> size m.
void DealDefault(vector<string>& result);

/// Expand below
///     TIME 3*36.525 WIR  12000  BHP  12000
/// to
///     TIME 36.525 WIR  12000  BHP  12000
///     TIME 36.525 WIR  12000  BHP  12000
///     TIME 36.525 WIR  12000  BHP  12000
vector<vector<string>> ExpandWellOptions(const vector<string>& result);

// 0: like 1990-01-01
// 1: a double
int WhichDateFormat(const std::string& str);

int NumDaysBetweenDates(const string& str1, const string& str2);


/// DealData change a series of product of integers into two arrays.
/// For example, 8*1  16*2  8*3  16*4  -> obj <8, 16, 8, 16> & val <1, 2, 3, 4>.
template <typename T>
void DealData(const vector<string>& vbuf, vector<OCP_USI>& obj, vector<T>& region)
{
    obj.resize(0);
    region.resize(0);
    for (auto& str : vbuf) {
        auto pos = str.find('*');
        if (pos != string::npos) {
            USI     len = str.size();
            OCP_USI num = stoi(str.substr(0, pos));
            USI     val = stoi(str.substr(pos + 1, len - (pos + 1)));
            obj.push_back(num);
            region.push_back(val);
        }
    }
}

/// Check if the stream starts with @a comment_char. If so skip it.
inline void skip_comment_lines(std::istream &is, const char comment_char)
{
    while (1)
    {
        is >> std::ws;
        if (is.peek() != comment_char)
        {
            break;
        }
        is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}

/// Erase white-space characters in a string
inline std::vector<std::string> strip_split(std::string& str)
{
    // 去除首部空白字符
    size_t start = str.find_first_not_of(" \t\n\r");
    if (start != std::string::npos)
        str.erase(0, start);
    else
        str.clear(); // 如果字符串全是空白字符，则清空字符串

    // 去除尾部空白字符
    size_t end = str.find_last_not_of(" \t\n\r");
    if (end != std::string::npos)
        str.erase(end + 1);

    // 使用 std::istringstream 按空格分割字符串
    std::istringstream iss(str);
    std::vector<std::string> words;
    std::string tmp;
    while (iss >> tmp)
        words.push_back(tmp);

    // Erase comments
    int pos;
    for (pos=0; pos<words.size(); ++pos)
    {
        if (words[pos].find('#') == 0 || words[pos].find("--") == 0)
            break;
    }
    words.erase(words.begin() + pos, words.end());

    return words;
}

/// Save the current file pointer and get a line.
inline std::istream& GetLineSkipComments(std::istream& is, std::string& buff, const char comment_char = '#')
{
    skip_comment_lines(is, '#');
    return std::getline(is, buff);
}

/// Judge the str is a double or string
inline bool isRegularString(const std::string& str)
{
    std::istringstream iss(str);
    double num;
    iss >> std::noskipws >> num;

    // 如果转换成功，说明是double类型
    if (!iss.fail() && iss.eof()) {
        return false;
    } else {
        return true;
    }
}

inline bool isInteger(const std::string& s)
{
    if(s.empty() || ( (!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+')) )
        return false;

    char* p;
    strtol(s.c_str(), &p, 10);

    return (*p == 0);
}


#endif /* end if __UTILINPUT_HEADER__ */

/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/