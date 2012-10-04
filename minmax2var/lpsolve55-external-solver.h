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
//  Created by Carmi Grushko on 6/21/12.

#ifndef minmax2var_lp_solver_stuff_h
#define minmax2var_lp_solver_stuff_h

#include "external_solver.h"
#include "settings.h"

#ifdef EXTERNAL_SOLVER_LPSOLVE55
#include "lp_lib.h"

void setConstraint(lprec* lp, int row, double a, double b, double c);

lprec* setup_unconstrained_lp(int N);

class lpsolve55_Solver : public ExternalSolver {
protected:
  lprec* lp;
  int n;
public:
  lpsolve55_Solver();
  virtual void init(int n);
  virtual void solve(const hough_min_max::P3 points[]);
  virtual std::string name() const;
  virtual ~lpsolve55_Solver();
};

// A faster, but less stable version, with fall-back to the
// safe version.
// NOTE: This class assumes the problem must be bounded and feasible !

#ifndef ETA_PFI_LIBRARY
#warning "ETA_PFI_LIBRARY is undefined - lpsolve55-etaPFI will not be used. The define should contain the path of bfp_etaPFI.{dll/so}. The Makefile present in the root directory should take care of that, if building from command-line. Otherwise, set it in the project's preprocessor directorives section (Visual Studio, Xcode, Eclipse, etc.)"
#else
class lpsolve55_Solver_etaPFI : public lpsolve55_Solver {
  virtual void init(int n);
  virtual void solve(const hough_min_max::P3 points[]);
  virtual std::string name() const;
};
#endif
#else
typedef NullExternalSolver lpsolve55_Solver;
typedef NullExternalSolver lpsolve55_Solver_etaPFI;
#endif

#endif
