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
//  Created by Carmi Grushko on 7/25/12.

#ifndef minmax2var_performance_data_h
#define minmax2var_performance_data_h

#include <fstream>
#include "settings.h"
#include "misc.h"

namespace hough_min_max {

// Do NOT instatiate
class PerformanceData {
#ifdef PERFORMANCE_DATA_3D
  std::ofstream data_file_3d;
#endif
#ifdef PERFORMANCE_DATA_2D
  std::ofstream data_file_2d;
#endif
public:
  PerformanceData();
  void add_scans_count_3d(size_t scansCount);
  void add_scans_count_2d(size_t scansCount);
  void start_new_measurement();
};

extern PerformanceData performanceDataSingleton;

} // namespace hough_min_max

#endif
