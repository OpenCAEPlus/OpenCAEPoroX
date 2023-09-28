/*! \file    OCPTable.cpp
 *  \brief   OCPTable class declaration
 *  \author  Shizhe Li
 *  \date    Oct/01/2021
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "OCPTable.hpp"


OCPTable::OCPTable(const vector<vector<OCP_DBL>>& src) { this->Setup(src); }

void OCPTable::Setup(const std::vector<std::vector<OCP_DBL>>& src)
{
    data = src;
    nCol = data.size();
    nRow = data[0].size();
    bId  = nRow / 2;
}


/// Careful: the memory outdata and slope have not be allocated before
USI OCPTable::Eval_All(const USI&       j,
                       const OCP_DBL&   val,
                       vector<OCP_DBL>& outdata,
                       vector<OCP_DBL>& slope) const
{
    if (val >= data[j][bId]) {
        for (USI i = bId + 1; i < nRow; i++) {
            if (val < data[j][i]) {
                bId = i - 1;
                for (USI k = 0; k < nCol; k++) {
                    slope[k] = (data[k][bId + 1] - data[k][bId]) /
                               (data[j][bId + 1] - data[j][bId]);
                    outdata[k] = data[k][bId] + slope[k] * (val - data[j][bId]);
                }
                return bId;
            }
        }
        bId = nRow - 1;
    } else {
        for (OCP_INT i = bId - 1; i >= 0; i--) {
            if (val >= data[j][i]) {
                bId = i;
                for (USI k = 0; k < nCol; k++) {
                    slope[k] = (data[k][bId + 1] - data[k][bId]) /
                               (data[j][bId + 1] - data[j][bId]);
                    outdata[k] = data[k][bId] + slope[k] * (val - data[j][bId]);
                }
                return bId;
            }
        }
        bId = 0;
    }
    for (USI k = 0; k < nCol; k++) {
        slope[k] = 0;
        outdata[k] = data[k][bId];
    }
    return bId;
}

USI OCPTable::Eval_All(const USI& j, const OCP_DBL& val, vector<OCP_DBL>& outdata) const
{
    OCP_DBL slope = 0;
    if (val >= data[j][bId]) {
        for (USI i = bId + 1; i < nRow; i++) {
            if (val < data[j][i]) {
                bId = i - 1;
                for (USI k = 0; k < nCol; k++) {
                    slope = (data[k][bId + 1] - data[k][bId]) /
                            (data[j][bId + 1] - data[j][bId]);
                    outdata[k] = data[k][bId] + slope * (val - data[j][bId]);
                }
                return bId;
            }
        }
        bId = nRow - 1;
    }
    else {
        for (OCP_INT i = bId - 1; i >= 0; i--) {
            if (val >= data[j][i]) {
                bId = i;
                for (USI k = 0; k < nCol; k++) {
                    slope = (data[k][bId + 1] - data[k][bId]) /
                        (data[j][bId + 1] - data[j][bId]);
                    outdata[k] = data[k][bId] + slope * (val - data[j][bId]);
                }
                return bId;
            }
        }
        bId = 0;
    }
    for (USI k = 0; k < nCol; k++) {
        outdata[k] = data[k][bId];
    }
    return bId;
}

USI OCPTable::Eval_All0(const OCP_DBL& val, vector<OCP_DBL>& outdata) const
{
    const USI j    = 0;
    OCP_DBL   tmpk = 0;
    if (val >= data[j][bId]) {
        for (USI i = bId + 1; i < nRow; i++) {
            if (val < data[j][i]) {
                bId = i - 1;
                for (USI k = 1; k < nCol; k++) {
                    tmpk = (data[k][bId + 1] - data[k][bId]) /
                           (data[j][bId + 1] - data[j][bId]);
                    outdata[k - 1] = data[k][bId] + tmpk * (val - data[j][bId]);
                }
                return bId;
            }
        }
        bId = nRow - 1;
    } else {
        for (OCP_INT i = bId - 1; i >= 0; i--) {
            if (val >= data[j][i]) {
                bId = i;
                for (USI k = 1; k < nCol; k++) {
                    tmpk = (data[k][bId + 1] - data[k][bId]) /
                           (data[j][bId + 1] - data[j][bId]);
                    outdata[k - 1] = data[k][bId] + tmpk * (val - data[j][bId]);
                }
                return bId;
            }
        }
        bId = 0;
    }
    for (USI k = 1; k < nCol; k++) {
        outdata[k - 1] = data[k][bId];
    }
    return bId;
}


USI OCPTable::Eval_All0(const OCP_DBL& val, vector<OCP_DBL>& outdata, vector<OCP_DBL>& slope) const
{
    const USI j = 0;
    if (val >= data[j][bId]) {
        for (USI i = bId + 1; i < nRow; i++) {
            if (val < data[j][i]) {
                bId = i - 1;
                for (USI k = 1; k < nCol; k++) {
                    slope[k - 1] = (data[k][bId + 1] - data[k][bId]) /
                        (data[j][bId + 1] - data[j][bId]);
                    outdata[k - 1] = data[k][bId] + slope[k - 1] * (val - data[j][bId]);
                }
                return bId;
            }
        }
        bId = nRow - 1;
    }
    else {
        for (OCP_INT i = bId - 1; i >= 0; i--) {
            if (val >= data[j][i]) {
                bId = i;
                for (USI k = 1; k < nCol; k++) {
                    slope[k - 1] = (data[k][bId + 1] - data[k][bId]) /
                        (data[j][bId + 1] - data[j][bId]);
                    outdata[k - 1] = data[k][bId] + slope[k - 1] * (val - data[j][bId]);
                }
                return bId;
            }
        }
        bId = 0;
    }
    for (USI k = 1; k < nCol; k++) {
        slope[k - 1]   = 0;
        outdata[k - 1] = data[k][bId];
    }
    return bId;
}


OCP_DBL OCPTable::Eval(const USI& j, const OCP_DBL& val, const USI& destj) const
{
    if (val >= data[j][bId]) {
        for (USI i = bId + 1; i < nRow; i++) {
            if (val < data[j][i]) {
                bId       = i - 1;
                OCP_DBL k = (data[destj][bId + 1] - data[destj][bId]) /
                            (data[j][bId + 1] - data[j][bId]);
                return (data[destj][bId] + k * (val - data[j][bId]));
            }
        }
        bId = nRow - 1;
        
    } else {
        for (OCP_INT i = bId - 1; i >= 0; i--) {
            if (val >= data[j][i]) {
                bId       = i;
                OCP_DBL k = (data[destj][bId + 1] - data[destj][bId]) /
                            (data[j][bId + 1] - data[j][bId]);
                return (data[destj][bId] + k * (val - data[j][bId]));
            }
        }
        bId = 0;
    }
    return data[destj][bId];
}

OCP_DBL OCPTable::Eval(const USI& j, const OCP_DBL& val, const USI& destj, OCP_DBL& myK) const
{
    if (val >= data[j][bId]) {
        for (USI i = bId + 1; i < nRow; i++) {
            if (val < data[j][i]) {
                bId = i - 1;
                myK = (data[destj][bId + 1] - data[destj][bId]) /
                      (data[j][bId + 1] - data[j][bId]);
                return (data[destj][bId] + myK * (val - data[j][bId]));
            }
        }
        bId = nRow - 1;
    } else {
        for (OCP_INT i = bId - 1; i >= 0; i--) {
            if (val >= data[j][i]) {
                bId = i;
                myK = (data[destj][bId + 1] - data[destj][bId]) /
                      (data[j][bId + 1] - data[j][bId]);
                return (data[destj][bId] + myK * (val - data[j][bId]));
            }
        }
        bId = 0;       
    }
    return data[destj][bId];
}

OCP_DBL OCPTable::Eval_Inv(const USI& j, const OCP_DBL& val, const USI& destj) const
{
    if (val > data[j][bId]) {
        for (OCP_INT i = bId - 1; i >= 0; i--) {
            if (val <= data[j][i]) {
                bId       = i;
                OCP_DBL k = (data[destj][bId + 1] - data[destj][bId]) /
                            (data[j][bId + 1] - data[j][bId]);
                return (data[destj][bId] + k * (val - data[j][bId]));
            }
        }
        bId = 0;
    } else {
        for (USI i = bId + 1; i < nRow; i++) {
            if (val >= data[j][i]) {
                bId       = i;
                OCP_DBL k = (data[destj][bId] - data[destj][bId - 1]) /
                            (data[j][bId] - data[j][bId - 1]);
                return (data[destj][bId - 1] + k * (val - data[j][bId - 1]));
            }
        }
        bId = nRow - 1;      
    }
    return data[destj][bId];
}


void OCPTable::GetCloseRow(const USI& j, const OCP_DBL& val, vector<OCP_DBL>& outdata) const
{
    USI outrow = 0;
    if (val >= data[j][bId]) {
        outrow = nRow - 1;
        for (USI i = bId + 1; i < nRow; i++) {
            if (val < data[j][i]) {
                bId = i - 1;
                // choose the cloest one
                if ((val - data[j][bId]) > (data[j][bId + 1] - val)) outrow = bId + 1;
                else                                                 outrow = bId;
                break;
            }
        }
        bId = nRow - 1;
    }
    else {
        outrow = 0;
        for (OCP_INT i = bId - 1; i >= 0; i--) {
            if (val >= data[j][i]) {
                bId = i;
                // choose the cloest one
                if ((val - data[j][bId]) > (data[j][bId + 1] - val)) outrow = bId + 1;
                else                                                 outrow = bId;
                break;
            }
        }
        bId = 0;
    }
    for (USI k = 0; k < nCol; k++) {
        outdata[k] = data[k][outrow];
    }
}


void OCPTable::Display() const
{
    cout << "\n---------------------" << endl
         << "Pressure Distribution"
         << "\n---------------------" << endl;
    cout << "   Depth         "
         << "   Poil         "
         << "   Pgas         "
         << "   Pwat         " << endl;

    cout << fixed << setprecision(3);
    for (USI i = 0; i < nRow; i++) {
        for (USI j = 0; j < nCol; j++) {
            cout << data[j][i] << "\t";
        }
        cout << "\n";
    }
}


void OCPTable2::Setup(const Table2& tab)
{
    numtable = tab.data.size();
    ref      = tab.refData;
    ref.resize(numtable, 0);

    for (USI i = 0; i < numtable; i++) {
        tables.push_back(OCPTable(tab.data[i]));
    }

    nCol = tab.colNum;
    data1.resize(nCol);
    data2.resize(nCol);
    cdata1.resize(nCol);
    cdata2.resize(nCol);
}


void OCPTable2::Eval_All(const OCP_DBL& val1, const OCP_DBL& val2, const USI& j2, vector<OCP_DBL>& out) const
{
    if (numtable == 1)             tables[0].Eval_All(j2, val2, out);
    else if (val1 <= ref.front())  tables.front().Eval_All(j2, val2, out);
    else if (val1 >= ref.back())   tables.back().Eval_All(j2, val2, out);
    else {
    for (USI i = 0; i < numtable - 1; i++) {
        if (val1 <= ref[i + 1]) {
            const OCP_DBL w = (val1 - ref[i]) / (ref[i + 1] - ref[i]);
            tables[i].Eval_All(j2, val2, data1, cdata1);
            tables[i + 1].Eval_All(j2, val2, data2, cdata2);

            for (USI j = 0; j < nCol; j++) {
                out[j] = w * data1[j] + (1 - w) * data2[j];
            }
            break;
        }
    }
    }
}


void OCPTable2::Eval_All(const OCP_DBL& val1, const OCP_DBL& val2, const USI& j2, vector<OCP_DBL>& out,
    vector<OCP_DBL>& slope1, vector<OCP_DBL>& slope2) const
{
    fill(slope1.begin(), slope1.end(), 0.0);
    if (numtable == 1)             tables[0].Eval_All(j2, val2, out, slope2);
    else if (val1 <= ref.front())  tables.front().Eval_All(j2, val2, out, slope2);
    else if (val1 >= ref.back())   tables.back().Eval_All(j2, val2, out, slope2);
    else {
        for (USI i = 0; i < numtable - 1; i++) {
            if (val1 <= ref[i + 1]) {
                const OCP_DBL w = (val1 - ref[i]) / (ref[i + 1] - ref[i]);
                tables[i].Eval_All(j2, val2, data1, cdata1);
                tables[i + 1].Eval_All(j2, val2, data2, cdata2);

                for (USI j = 0; j < nCol; j++) {
                    out[j]    = w * data1[j] + (1 - w) * data2[j];
                    slope2[j] = w * cdata1[j] + (1 - w) * cdata2[j];
                    slope1[j] = (data2[j] - data1[j]) / (ref[i + 1] - ref[i]);
                }
                break;
            }
        }
    }
}


void OCPTable2::Eval_All0(const OCP_DBL& val1, const OCP_DBL& val2, vector<OCP_DBL>& out) const
{
    if (numtable == 1)             tables[0].Eval_All0(val2, out);
    else if (val1 <= ref.front())  tables.front().Eval_All0(val2, out);
    else if (val1 >= ref.back())   tables.back().Eval_All0(val2, out);
    else {
        for (USI i = 0; i < numtable - 1; i++) {
            if (val1 <= ref[i + 1]) {
                const OCP_DBL w = (val1 - ref[i]) / (ref[i + 1] - ref[i]);
                tables[i].Eval_All0(val2, data1, cdata1);
                tables[i + 1].Eval_All0(val2, data2, cdata2);

                for (USI j = 0; j < nCol - 1; j++) {
                    out[j] = w * data1[j] + (1 - w) * data2[j];
                }
                break;
            }
        }
    }
}


void OCPTable2::Eval_All0(const OCP_DBL& val1, const OCP_DBL& val2, vector<OCP_DBL>& out,
    vector<OCP_DBL>& slope1, vector<OCP_DBL>& slope2) const
{
    fill(slope1.begin(), slope1.end(), 0.0);
    if (numtable == 1)             tables[0].Eval_All0(val2, out, slope2);
    else if (val1 <= ref.front())  tables.front().Eval_All0(val2, out, slope2);
    else if (val1 >= ref.back())   tables.back().Eval_All0(val2, out, slope2);
    else {
        for (USI i = 0; i < numtable - 1; i++) {
            if (val1 <= ref[i + 1]) {
                const OCP_DBL w = (val1 - ref[i]) / (ref[i + 1] - ref[i]);
                tables[i].Eval_All0(val2, data1, cdata1);
                tables[i + 1].Eval_All0(val2, data2, cdata2);

                for (USI j = 0; j < nCol - 1; j++) {
                    out[j]    = w * data1[j] + (1 - w) * data2[j];
                    slope2[j] = w * cdata1[j] + (1 - w) * cdata2[j];
                    slope1[j] = (data2[j] - data1[j]) / (ref[i + 1] - ref[i]);
                }
                break;
            }
        }
    }
}


OCP_DBL OCPTable2::Eval(const OCP_DBL& val1, const OCP_DBL& val2, const USI& j2, const USI& destj) const
{
    if (numtable == 1)             return tables[0].Eval(j2, val2, destj);
    else if (val1 <= ref.front())  return tables.front().Eval(j2, val2, destj);
    else if (val1 >= ref.back())   return tables.back().Eval(j2, val2, destj);
    else {
        for (USI i = 0; i < numtable - 1; i++) {
            if (val1 <= ref[i + 1]) {
                const OCP_DBL w  = (val1 - ref[i]) / (ref[i + 1] - ref[i]);
                const OCP_DBL d1 = tables[i].Eval(j2, val2, destj);
                const OCP_DBL d2 = tables[i + 1].Eval(j2, val2, destj);

                return w * d1 + (1 - w) * d2;
            }
        }
    }
    OCP_ABORT("WRONG EVAL!");
}


OCP_DBL OCPTable2::Eval(const OCP_DBL& val1, const OCP_DBL& val2, const USI& j2, const USI& destj, OCP_DBL& myK1, OCP_DBL& myK2) const
{
    myK1 = 0;
    if (numtable == 1)             return tables[0].Eval(j2, val2, destj, myK2);
    else if (val1 <= ref.front())  return tables.front().Eval(j2, val2, destj, myK2);
    else if (val1 >= ref.back())   return tables.back().Eval(j2, val2, destj, myK2);
    else {
        for (USI i = 0; i < numtable - 1; i++) {
            if (val1 <= ref[i + 1]) {
                const OCP_DBL w  = (val1 - ref[i]) / (ref[i + 1] - ref[i]);
                const OCP_DBL d1 = tables[i].Eval(j2, val2, destj, myK1);
                const OCP_DBL d2 = tables[i + 1].Eval(j2, val2, destj, myK2);

                myK2 = w * myK1 + (1 - w) * myK2;
                myK1 = (d2 - d1) / (ref[i + 1] - ref[i]);

                return (w * d1 + (1 - w) * d2);
            }
        }
    }
    OCP_ABORT("WRONG EVAL!");
}


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Oct/01/2021      Create file                          */
/*  Chensong Zhang      Oct/15/2021      Format file                          */
/*----------------------------------------------------------------------------*/