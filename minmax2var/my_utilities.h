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
//  Created by Carmi Grushko on 7/1/12.

#ifndef minmax2var_my_utilities_h
#define minmax2var_my_utilities_h

#include <limits>
#include <math.h>
#include "settings.h"

extern "C" double orient2d(const double *pa, const double *pb, const double *pc);
extern "C" double orient2dfast(const double *pa, const double *pb, const double *pc);
extern "C" double orient3d(const double *pa, const double *pb, const double *pc, const double* pd);

namespace hough_min_max {

  inline double filtered_orient2d(const double a[], const double b[], const double c[]) {
#ifdef FILTER_EPA
    double result = ::orient2dfast(a, b, c);
    if (fabs(result) > FILTER_EPA)
      return result;
#endif
    return ::orient2d(a, b, c);
  }

  inline double filtered_det2(const double a[], const double b[]) {
    static const double zero2[2] = {0, 0};    
    return filtered_orient2d(a, b, zero2);
  }
  
  struct Result {
    double u, v, value;
    bool operator<(const Result& rhs) const {
      return value < rhs.value;
    }
	Result() : u(0), v(0), value(std::numeric_limits<double>::infinity()) {
    }
  };
  
  struct Rational {
    union {
      struct {
        double num, denom;
      };
      double elements[2];
    };

    Rational(double num = 0, double denom = 1) : num(num), denom(denom) {
      internal_verification(denom != 0);
    }
    bool operator<(const Rational& rhs) const {
      bool reverse_sign = (denom < 0) != (rhs.denom < 0);
      return (filtered_det2(elements, rhs.elements) > 0) != reverse_sign;
//      
//      if ((denom < 0) != (rhs.denom < 0))
//        return num * rhs.denom > denom * rhs.num + h * denom * rhs.denom;
//      
//      return num * rhs.denom < denom * rhs.num - h * denom * rhs.denom;
    }
    bool operator<=(const Rational& rhs) const {
//      return (fabs(num - rhs.num) < h && fabs(denom - rhs.denom) < h) || *this < rhs;
      return filtered_det2(elements, rhs.elements) == 0 || *this < rhs;
    }
  };
  
} // namespace hough_min_max

#endif
