/*! \file    UtilMemory.cpp
 *  \brief   UtilMemory class declaration
 *  \author  Shizhe Li
 *  \date    Jul/26/2024
 *
 *-----------------------------------------------------------------------------------
 *  Copyright (C) 2021--present by the OpenCAEPoroX team. All rights reserved.
 *  Released under the terms of the GNU Lesser General Public License 3.0 or later.
 *-----------------------------------------------------------------------------------
 */

#include "UtilMemory.hpp"


#if defined(_CONSOLE) || defined(_WIN32) || defined(_WIN64)

double OCPGetCurrentRSS() {
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
        return pmc.WorkingSetSize / (1024.0 * 1024.0 * 1024.0);  // Convert to GB
    }
    return 0.0;  // If not found
}


#else

double OCPGetCurrentRSS() {
    std::ifstream file("/proc/self/status");
    std::string line;
    while (std::getline(file, line)) {
        if (line.find("VmRSS:") == 0) {
            std::istringstream iss(line);
            std::string key, value, unit;
            iss >> key >> value >> unit;
            return std::stoull(value) * 1024 / (1024.0 * 1024.0 * 1024.0);  // Convert to GB
        }
    }
    return 0.0;  // If not found
}

#endif


/*----------------------------------------------------------------------------*/
/*  Brief Change History of This File                                         */
/*----------------------------------------------------------------------------*/
/*  Author              Date             Actions                              */
/*----------------------------------------------------------------------------*/
/*  Shizhe Li           Jul/26/2024      Create file                          */
/*----------------------------------------------------------------------------*/