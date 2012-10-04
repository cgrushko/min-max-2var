/*
 * This file is part of minmax2var - A two-variable linar min-max solver.
 * 
 * Copyright (C) 2012  Carmi Grushko, carmi.grushko@gmail.com
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. 
 */

#include "performance-data.h"
#include <time.h>
#include <string>
#include <stdlib.h>

namespace hough_min_max {

PerformanceData performanceDataSingleton;

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
static const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://www.cplusplus.com/reference/clibrary/ctime/strftime/
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
}

PerformanceData::PerformanceData() 
  #ifdef PERFORMANCE_DATA_3D
  : data_file_3d((getUserDirectory() + "/" + PERFORMANCE_DATA_3D).c_str(), std::ios::app)
  #endif
  #ifdef PERFORMANCE_DATA_2D
  , data_file_2d((getUserDirectory() + "/" + PERFORMANCE_DATA_2D).c_str(), std::ios::app)
  #endif
{
  start_new_measurement();
}

void PerformanceData::add_scans_count_3d(size_t scansCount) {
  #ifdef PERFORMANCE_DATA_3D
  data_file_3d << scansCount << ' ';
  #endif
}

void PerformanceData::add_scans_count_2d(size_t scansCount) {
  #ifdef PERFORMANCE_DATA_2D
  data_file_2d << scansCount << ' ';
  #endif
}

void PerformanceData::start_new_measurement() {
  #ifdef PERFORMANCE_DATA_3D
  data_file_3d << '\n' << currentDateTime() << '\n';
  #endif
  #ifdef PERFORMANCE_DATA_2D
  data_file_2d << '\n' << currentDateTime() << '\n';
  #endif
}

} // namespace hough_min_max
