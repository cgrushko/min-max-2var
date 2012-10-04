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
//  Created by Carmi Grushko on 8/5/12.

#ifndef minmax2var_external_solver_h
#define minmax2var_external_solver_h

#include "P3.h"
#include <string>
#include <math.h>
#include <ostream>
#include <stdlib.h>

class ExternalSolver {
public:
  double time; // Subclasses should update this to reflect total run time in ms.
  double u, v, value; // Subclasses should update this in solve()
  bool bounded; // Subclasses should update this in solve()
public:
  // n - number of constraints
  virtual void init(int n) = 0;

  ExternalSolver() : time(0) { 
  }
  
  virtual void solve(const hough_min_max::P3 points[]) = 0;
  virtual std::string name() const = 0;
  
  virtual void print_status(std::ostream& ostr) const {
    if (bounded) {
      ostr << name() << "\t: U, V, Value: " << u << " " << v << " " << value << std::endl;
    } else {
      ostr << name() << "\t: Unbounded" << std::endl;
    }
  }
  
  virtual bool assert_solution(bool bounded, double u, double v, double value) {
    static const double h = 1e-5;
    return (bounded == this->bounded) && (fabs(value-this->value) < h);
  }
  
  virtual ~ExternalSolver() {
  }
};

class NullExternalSolver : public ExternalSolver {
  virtual void init(int n){
  }
  virtual bool assert_solution(bool bounded, double u, double v, double value) {
    return true;
  }
  virtual void solve(const hough_min_max::P3 points[]) {
  }
  virtual void print_status(std::ostream& ostr) const {
  }
  virtual std::string name() const {
    return "";
  }
};

#endif
