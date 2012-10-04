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
//  Created by Carmi Grushko on 6/29/12.

#include "solver2d.h"
#include "P2.h"
#include "P3.h"

#include <limits>

#include "performance-data.h"
#include "openmp-reductor.h"

namespace hough_min_max {
  namespace solver_2d {
    
    template <int Coordinate>
    V2 lineNormalFrom2Points(const P3& p0, const P3& p1) {
      V2 res(-p1.z + p0.z, p1.coord[Coordinate] - p0.coord[Coordinate]);
      if (res.y < 0) {
        res.x = -res.x;
        res.y = -res.y;
      }
      return res;
    }

    template <int Coordinate>
    bool assertDualSolution(const P3 ptArray[], int n, const P3& p0, const P3& p1) {
      for (int i=0; i<n; i++) {
        const P3& p = ptArray[i];
        double p0_coord[2] = {p0.coord[Coordinate], p0.z};
        double p1_coord[2] = {p1.coord[Coordinate], p1.z};
        double p_coord[2] = {p.coord[Coordinate], p.z};
        release_assert((p0_coord[0] > p1_coord[0]) !=
                       (filtered_orient2d(p0_coord, p1_coord, p_coord) >= 0));
      }
      return true;
    }
    
    typedef std::pair<bool, bool> PBB;
    
    void find_max_on_corners(P3 const *a, int n, Result *result_corner_x, Result *result_corner_y, Result *result_corner_0) {
      if (n == 0)
        return;
      result_corner_x->u = 1;  // corner (1,0)
      result_corner_y->v = 1;  // corner (0,1)
      
      result_corner_x->value = a[0].x - a[0].z;
      result_corner_y->value = a[0].y - a[0].z;
      result_corner_0->value = a[0].z;

      for (int i=0; i < n; ++i) {
        if (a[i].z < result_corner_0->value) {
          result_corner_0->value = a[i].z;
        }
        if (a[i].x - a[i].z > result_corner_x->value) {
          result_corner_x->value = a[i].x - a[i].z;
        }
        if (a[i].y - a[i].z > result_corner_y->value) {
          result_corner_y->value = a[i].y - a[i].z;
        }
      }

      result_corner_0->value = -result_corner_0->value;
    }

    template <int Coordinate> // 0 for .x, 1 for .y
    bool findInitialGuess(const P3* ptArray, int n, const P3** p0, const P3** p1) {
      // Easy: find one with positive x (or y), one with negative x (resp. y)
      int negative_i = -1;
      int positive_i = -1;
      for (int i=0; i<n && (negative_i == -1 || positive_i == -1); i++) {
        if (negative_i == -1 && ptArray[i].coord[Coordinate] < 0) {
          negative_i = i;
          continue;
        }
        if (positive_i == -1 && ptArray[i].coord[Coordinate] > 0)
          positive_i = i;
      }
      
      if (negative_i == -1 || positive_i == -1) {
        for (int i=0; i<n && (negative_i == -1 || positive_i == -1); i++) {
          if (negative_i == -1 && ptArray[i].coord[Coordinate] <= 0)
            negative_i = i;
          if (positive_i == -1 && ptArray[i].coord[Coordinate] >= 0)
            positive_i = i;
        }      
      }
      
      if (negative_i != -1 && positive_i != -1 && positive_i != negative_i) {
        *p0 = &ptArray[negative_i];
        *p1 = &ptArray[positive_i];
        return true;
      }
      return false;
    }
    
    template <int Coordinate> // 0 for .x, 1 for .y
    void dualSolutionToMinMaxSolution(const P3& p0, const P3& p1, double* u, double* val) {
      double Nz = p1.coord[Coordinate] - p0.coord[Coordinate];
      *u = (p1.z - p0.z) / Nz;
      *val = - (p1.coord[Coordinate]*p0.z - p1.z*p0.coord[Coordinate]) / Nz;
    }
    
    template <int Coordinate>
    void newPointIsBelowLine(const P3** p0, const P3** p1, const P3& p) {
      if ((p.coord[Coordinate] < 0) == (*p0)->coord[Coordinate] < 0) {
        *p0 = &p;
      } else {
        *p1 = &p;
      }
    }

    void letItFall_on_edges_worker(const P3& new_point, bool run_x, bool run_y, const P3** p0_x, const P3** p1_x, const P3** p0_y, const P3** p1_y, bool *changed_x, bool *changed_y, V2* N_x, V2* N_y) {
      if (run_x) {
        if (PtIsBelowLine<P3::X>(**p0_x, **p1_x, *N_x, new_point)) {
          newPointIsBelowLine<P3::X>(p0_x, p1_x, new_point);
          *changed_x = true;
          *N_x = lineNormalFrom2Points<P3::X>(**p0_x, **p1_x);
        }
      }
      
      if (run_y) {
        if (PtIsBelowLine<P3::Y>(**p0_y, **p1_y, *N_y, new_point)) {
          newPointIsBelowLine<P3::Y>(p0_y, p1_y, new_point);
          *changed_y = true;
          *N_y = lineNormalFrom2Points<P3::Y>(**p0_y, **p1_y);
        }
      }
    }
    
    PBB letItFall_on_edges(const P3* ptArray, int n, bool run_x, bool run_y, const P3** p0_x, const P3** p1_x, const P3** p0_y, const P3** p1_y) {
      V2 N_x = run_x ? lineNormalFrom2Points<P3::X>(**p0_x, **p1_x) : V2();
      V2 N_y = run_y ? lineNormalFrom2Points<P3::Y>(**p0_y, **p1_y) : V2();
      bool changed_x = false;
      bool changed_y = false;
      
      for (int i=0; i<n; i++) {
        letItFall_on_edges_worker(ptArray[i], run_x, run_y, p0_x, p1_x, p0_y, p1_y, &changed_x, &changed_y, &N_x, &N_y);
      }
      return PBB(changed_x, changed_y);
    }
#ifdef PARALLEL
    typedef std::pair<P3 const *, P3 const *> PPP;
    typedef std::pair<V2, V2> PVV;
    template <int Coordinate>
    void letItFall_worker(const P3& new_point, const P3** p0, const P3** p1, bool *changed, V2 *N);
    
    struct LetItFallBody_on_edges {
      PVV N;
      PBB changed;
      PBB run;
      PPP p_x;
      PPP p_y;
      
      void init(const PPP& p_x, const PPP& p_y, const PBB& run) {
        this->p_x = p_x;
        this->p_y = p_y;
        this->run = run;
        changed.first = changed.second = false;
        N.first = lineNormalFrom2Points<P3::X>(*(p_x.first), *(p_x.second));
        N.second = lineNormalFrom2Points<P3::Y>(*(p_y.first), *(p_y.second));
      }
      
      void operator()(const P3* iterator) {
        letItFall_on_edges_worker(*iterator, run.first, run.second, &p_x.first, &p_x.second, &p_y.first, &p_y.second, &changed.first, &changed.second, &N.first, &N.second);
      }
      
      void join( LetItFallBody_on_edges& rhs ) {
        if (run.first) {
          letItFall_worker<P3::X>(*rhs.p_x.first, &p_x.first, &p_x.second, &changed.first, &N.first);
          letItFall_worker<P3::X>(*rhs.p_x.second, &p_x.first, &p_x.second, &changed.first, &N.first);
        }
        if (run.second) {
          letItFall_worker<P3::Y>(*rhs.p_y.first, &p_y.first, &p_y.second, &changed.second, &N.second);
          letItFall_worker<P3::Y>(*rhs.p_y.second, &p_y.first, &p_y.second, &changed.second, &N.second);
        }
      }
    };
#endif

    void solve_on_edges(const P3* ptArray, int n, Result *result_x, Result* result_y) {
      PBB result;
      const P3* p0_x;
      const P3* p1_x;
      const P3* p0_y;
      const P3* p1_y;
      result.first = findInitialGuess<P3::X>(ptArray, n, &p0_x, &p1_x);
      result.second = findInitialGuess<P3::Y>(ptArray, n, &p0_y, &p1_y);
      
      if (!result.first && !result.second) {
        result_x->u = std::numeric_limits<double>::infinity();
        result_y->v = std::numeric_limits<double>::infinity();
        return;
      }

#ifdef PARALLEL
      static OpenMPReductor<const P3*, LetItFallBody_on_edges> reductor;
      LetItFallBody_on_edges body;
      body.init(PPP(p0_x, p1_x), PPP(p0_y, p1_y), result);
#endif
      
      PBB c1;
      PBB run = result;
      size_t counter = 0;
      do {
        counter++;
#ifndef PARALLEL
        c1 = letItFall_on_edges(ptArray, n, run.first, run.second, &p0_x, &p1_x, &p0_y, &p1_y);
#else
        body.changed = PBB(false, false);
        reductor.reduce(ptArray, ptArray + n, &body);
        c1 = body.changed;
#endif
        run.first = run.first && c1.first;
        run.second = run.second && c1.second;
      } while (run.first || run.second);
#ifdef PARALLEL
      p0_x = body.p_x.first;
      p1_x = body.p_x.second;
      p0_y = body.p_y.first;
      p1_y = body.p_y.second;
#endif
      
      performanceDataSingleton.add_scans_count_2d(counter);

      if (result.first) {
        dualSolutionToMinMaxSolution<P3::X>(*p0_x, *p1_x, &result_x->u, &result_x->value);
      } else {
        result_x->u = std::numeric_limits<double>::infinity();
      }
      if (result.second) {
        dualSolutionToMinMaxSolution<P3::Y>(*p0_y, *p1_y, &result_y->v, &result_y->value);
      } else {
        result_y->v = std::numeric_limits<double>::infinity();
      }
    }
    
    template <int Coordinate>
    void letItFall_worker(const P3& new_point, const P3** p0, const P3** p1, bool *changed, V2 *N) {
      if (PtIsBelowLine<Coordinate>(**p0, **p1, *N, new_point)) {
        newPointIsBelowLine<Coordinate>(p0, p1, new_point);
        *changed = true;
        *N = lineNormalFrom2Points<Coordinate>(**p0, **p1);
      }
    }
    
    template <int Coordinate>
    bool letItFall(const P3* ptArray, int n, const P3** p0, const P3** p1) {
      V2 N = lineNormalFrom2Points<Coordinate>(**p0, **p1);
      bool changed = false;
      
      for (int i=0; i<n; i++) {
        letItFall_worker<Coordinate>(ptArray[i], p0, p1, &changed, &N);
      }
      return changed;
    }    
        
#ifdef PARALLEL
    template <int Coordinate>
    struct LetItFallBody {
      V2 N;
      bool changed;
      P3 const *p0, *p1;
      
      void init(const P3* p0, const P3* p1) {
        this->p0 = p0; this->p1 = p1;
        changed = false;
        N = lineNormalFrom2Points<Coordinate>(*p0, *p1);
      }
      
      void operator()(const P3* iterator) {
        letItFall_worker<Coordinate>(*iterator, &p0, &p1, &changed, &N);
      }
      
      void join( LetItFallBody& rhs ) {
        letItFall_worker<Coordinate>(*rhs.p0, &p0, &p1, &changed, &N);
        letItFall_worker<Coordinate>(*rhs.p1, &p0, &p1, &changed, &N);
      }
    };
#endif
    
    void solve_on_x(const P3* ptArray, int n, Result *result_x) {
      P3 const *p0;
      P3 const *p1;
      if (!findInitialGuess<P3::X>(ptArray, n, &p0, &p1)) {
        result_x->u = std::numeric_limits<double>::infinity();
        return;
      }

#ifdef PARALLEL
      static OpenMPReductor<const P3*, LetItFallBody<P3::X> > reductor;
      LetItFallBody<P3::X> body;
      body.init(p0, p1);
#endif
      bool c1;
      size_t counter = 0;
      do {
        counter++;
#ifndef PARALLEL
        c1 = letItFall<P3::X>(ptArray, n, &p0, &p1);
#else
        body.changed = false;
        reductor.reduce(ptArray, ptArray + n, &body);
        c1 = body.changed;
#endif
      } while (c1);
#ifdef PARALLEL
      p0 = body.p0;
      p1 = body.p1;
#endif

      performanceDataSingleton.add_scans_count_2d(counter);
      internal_verification(assertDualSolution<P3::X>(ptArray, n, *p0, *p1));

      dualSolutionToMinMaxSolution<P3::X>(*p0, *p1, &result_x->u, &result_x->value);
    }
    
  } // namespace solver_2d
} // namespace hough_min_max

