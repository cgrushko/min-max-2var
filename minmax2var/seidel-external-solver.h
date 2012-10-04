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
//  Created by Carmi Grushko on 9/6/12.

#ifndef __minmax2var__seidel_external_solver__
#define __minmax2var__seidel_external_solver__

#include "external_solver.h"
#include "settings.h"

#ifdef EXTERNAL_SOLVER_SEIDEL
#include "lp.h"

class Seidel_Solver : public ExternalSolver {
private:
  int halves_count;
  static int numberOfExtraHalves();
  void preparePermutationArrays();
  void fillHalvesArray(const hough_min_max::P3 points[]);
  void addConstraint(FLOAT a, FLOAT b, FLOAT c, FLOAT d);
protected:
  int *next, *prev, *perm;
  FLOAT *work;
  
  // input to linprog
  FLOAT *halves;
  
  // = halves + <something>, where <something> is constraints such as
  // x>0, y>0, w>0 (the points are projective) and x+y < 1
  FLOAT *input_halves;
  
  FLOAT n_vec[4];
  FLOAT d_vec[4];
  FLOAT opt[4];
  int n;
  int d;
public:
  Seidel_Solver();
  virtual void init(int n);
  virtual void solve(const hough_min_max::P3 points[]);
  virtual std::string name() const;
  virtual ~Seidel_Solver();
};

#else
typedef NullExternalSolver Seidel_Solver;
#endif

#endif /* defined(__minmax2var__seidel_external_solver__) */
