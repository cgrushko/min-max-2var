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

#ifndef minmax2var_P3_h
#define minmax2var_P3_h

#include <stdio.h>
#include <inttypes.h>
#include <iostream>
#include <iomanip>
#include "P2.h"
#include "my_utilities.h"

namespace hough_min_max {
  
  const double h = 1e-9;
  
  struct P3 {
    enum {X = 0, Y = 1, Z = 2};
    
    union {
      struct {
        double x, y, z;
      };
      double coord[3];
    };
    
    P3(double x=0, double y=0, double z=0) : z(z), x(x), y(y)  {
    }
    
    static const P3 zero;
  };
  
  typedef P3 V3;
  
  inline P3 makePt(double x, double y, double z) {
    P3 p;
    p.x = x;
    p.y = y;
    p.z = z;
    return p;
  }
  
  inline double dot(const P3& p, const P3& q) {
    return p.x * q.x + p.y * q.y + p.z * q.z;
  }
  
  inline P3 cross(const P3& p, const P3& q) {
    return makePt(p.y*q.z - p.z*q.y, p.z*q.x - p.x*q.z, p.x*q.y - p.y*q.x);
  }
  
  inline P3 negate(const P3& p) {
    return makePt(-p.x, -p.y, -p.z);
  }
  
  inline P3 subtract(const P3& p, const P3& rhs) {
    return makePt(p.x - rhs.x, p.y - rhs.y, p.z - rhs.z);
  }
  
  // Discard a coordinate, make p into a P2 point. Keep .z coordinate
  template <int Coordinate>
  inline P2 project(const P3& p) {
    return P2(p.coord[Coordinate], p.z);
  }
  
  // Returns true iff cross(p1-p0, p2-p0) is facing down (negative z)
  inline bool normalFacesDown(const P3& p0, const P3& p1, const P3& p2) {
    return filtered_orient2d(p0.coord, p1.coord, p2.coord) < 0;
  }
  
  inline bool PtIsBelowPlane(const P3& p0, const P3& p1, const P3& p2, const V3& N, const P3& new_point) {
#ifdef FILTER_EPA
    double result = dot(N, subtract(new_point, p0));
    if (fabs(result) > FILTER_EPA)
      return result < 0;
#endif
    double det = ::orient3d(p0.coord, p1.coord, p2.coord, new_point.coord);
    return normalFacesDown(p0, p1, p2) ? (det < 0) : (det > 0);
  }
  
  // Points are P3, but are theoretically mapped into P2 points, with Coordinate and z.
  // p0 and p1 must define N; if (N * (p - p0)) is not big enough, use robust calculation
  template <int Coordinate>
  inline bool PtIsBelowLine(const P3& p0, const P3& p1, const V2& N, const P3& p) {
    return PtIsBelowLine(project<Coordinate>(p0),
                         project<Coordinate>(p1),
                         N,
                         project<Coordinate>(p));
  }


  inline void printPt(const P3& p, bool python_list = false) {
    if (python_list)
      printf("[%.16g, %.16g, %.16g]\n", p.x, p.y, p.z);
    else
      printf("%.16g %.16g %.16g\n", p.x, p.y, p.z);
  }

  inline void printPt_hex(const P3& p, std::ostream& stream, bool python_list = false) {
    uint64_t x = *reinterpret_cast<const uint64_t*>(&p.x);
    uint64_t y = *reinterpret_cast<const uint64_t*>(&p.y);
    uint64_t z = *reinterpret_cast<const uint64_t*>(&p.z);
    if (python_list)
      stream << "[";
    stream << std::hex << "0x" << x << (python_list ? ", " : " ");
    stream << std::hex << "0x" << y << (python_list ? ", " : " ");
    stream << std::hex << "0x" << z;
    if (python_list)
      stream << "]";
    stream << std::endl;
  }
  
  inline void print_problem_hex(const P3 points[], int n, std::ostream& stream) {
    for (int i=0; i<n; i++) {
      printPt_hex(points[i], stream);
    }
    stream << std::flush;
  }
  
  inline P3 readPt_hex(std::istream& stream) {
    uint64_t x, y, z;
    stream >> std::hex >> x >> std::hex >> y >> std::hex >> z;
    return P3(*reinterpret_cast<double*>(&x), *reinterpret_cast<double*>(&y), *reinterpret_cast<double*>(&z));
  }
  
  inline void printPt_bin(const P3& p, std::ostream& stream) {
    stream.write((char*)&p.x, sizeof(p.x));
    stream.write((char*)&p.y, sizeof(p.y));
    stream.write((char*)&p.z, sizeof(p.z));
  }
  
  inline P3 readPt_bin(std::istream& stream) {
    P3 result;
    stream.read((char*)&result.x, sizeof(result.x));
    stream.read((char*)&result.y, sizeof(result.y));
    stream.read((char*)&result.z, sizeof(result.z));
    return result;
  }
} // namespace hough_min_max

#endif
