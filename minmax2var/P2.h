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
//  Created by Carmi Grushko on 6/19/12.

#ifndef minmax2var_P2_h
#define minmax2var_P2_h

#include <stdio.h>
#include "my_utilities.h"

namespace hough_min_max {
  
  struct P2 {
    union {
      struct {
        double x, y;
      };
      double coord[2];
    };
    P2(double x = 0, double y = 0) : x(x), y(y) {
    }
  };
  
  typedef P2 V2;
  
  inline double dot(const P2& p, const P2& q) {
    return p.x * q.x + p.y * q.y;
  }
  
  inline P2 subtract(const P2& p, const P2& rhs) {
    return P2(p.x - rhs.x, p.y - rhs.y);
  }
  
  // p0 and p1 must define N; if (N * (p - p0)) is not big enough, use robust calculation
  inline bool PtIsBelowLine(const P2& p0, const P2& p1, const V2& N, const P2& p) {
#ifdef FILTER_EPA
    double result = dot(N, subtract(p, p0));
    if (fabs(result) > FILTER_EPA)
      return result < 0;
#endif
    double det = ::orient2d(p0.coord, p1.coord, p.coord);
    return (p0.x > p1.x) ? (det > 0) : (det < 0);
  }
  
  inline bool PtIsOnLine(const P2& p0, const P2& p1, const P2& p) {
    return ::orient2d(p0.coord, p1.coord, p.coord) == 0;
  }

  inline double cross(const P2& p, const P2& q) {
    return filtered_det2(p.coord, q.coord);
//    return p.x*q.y - p.y*q.x;
  }
  
  inline P2 rotate90(const P2& p) {
    return P2(-p.y, p.x);
  }
  
  inline void printPt(const P2& p) {
    printf("%.15lf, %.15lf", p.x, p.y);
  }
  
} // namespace hough_min_max

#endif
