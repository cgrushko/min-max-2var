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
//  Created by Carmi Grushko on 8/4/12.

#ifndef __minmax2var__openmp_reductor__
#define __minmax2var__openmp_reductor__

#include "settings.h"

#ifdef PARALLEL

#include <vector>

using std::vector;

#include <omp.h>

// ConstIterator should be RandomAccess for performance reasons
template <class ConstIterator, class Body>
class OpenMPReductor {
  vector<Body> data;
  size_t cache_line_size;
  Body* body;
public:
  OpenMPReductor(int cache_line_size = 1024)
  : data((omp_get_max_threads() - 1) * cache_line_size) // -1 because the user-supplied one is also used
  , cache_line_size(cache_line_size) {
  }

  void reduce(const ConstIterator& begin, const ConstIterator& end, Body* body) {
    size_t n = end - begin;
    
    // Initialize
    for (size_t i = 0; i < data.size(); i += cache_line_size) {
      data[i] = *body;
    }
    if (data.size() < (omp_get_max_threads() - 1) * cache_line_size) {
      data.resize((omp_get_max_threads() - 1) * cache_line_size, *body);
    }
    
    // Process
    size_t thread_count_global;
#pragma omp parallel
    {
      size_t thread_count = omp_get_num_threads();
      size_t thread_id = omp_get_thread_num();
      size_t work_per_thread = n / thread_count;
      Body* current_body = NULL;
      size_t begin_i, end_i;
      
      if (thread_id == 0) {
        begin_i = (thread_count - 1) * work_per_thread;
        end_i = n;
        current_body = body;
      } else {
        size_t temp = thread_id - 1;
        current_body = &data[temp * cache_line_size];
        begin_i = temp * work_per_thread;
        end_i = (temp+1) * work_per_thread;
      }
      
      for (size_t i = begin_i; i < end_i; ++i) {
        (*current_body)(begin + i);
      }
      
#pragma omp master
      {
        thread_count_global = thread_count;
      }
    }
    
    // Reduce
    for(size_t i=0; i < cache_line_size * (thread_count_global - 1); i += cache_line_size) {
      body->join(data[i]);
    }
  }

};

#endif

#endif /* defined(__minmax2var__openmp_reductor__) */
