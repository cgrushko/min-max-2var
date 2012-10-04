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
//  Created by Carmi Grushko on 7/8/12.

#ifndef minmax2var_settings_h
#define minmax2var_settings_h

// Preferably, set these in Makefile/Xcode project
//#define EXTERNAL_SOLVER_LPSOLVE55
//#define EXTERNAL_SOLVER_CGAL
//#define EXTERNAL_SOLVER_GLPK
//#define EXTERNAL_SOLVER_SEIDEL
//#define EXTERNAL_SOLVER_MOSEK

//#define DEBUG_PRINTF

// Verifies the solution is feasible, etc.
// Turn off for benchmarks, turn on for test runs.
//#define INTERNAL_VERIFICATIONS

// Record performance data: how many iterations each problem took.
// 1) Number of scans hough_min_max::solve3d() performed for each call. 
//#define PERFORMANCE_DATA_3D "hough-min-max-solver-perf-data-3d.txt"
// 2) Number of scans hough_min_max::solve2d() performed for each call.
//#define PERFORMANCE_DATA_2D "hough-min-max-solver-perf-data-2d.txt"

// Use OpenMP to parallelize the solver
//#define PARALLEL

// Define to compute fast and cheap geometric predicates (e.g., point above plane) and falling
// back to slow and safe extended-precision arithmetic if result is smaller in magnitude than
// value
#define FILTER_EPA 1e-5

#ifdef DEBUG_PRINTF
#define dbg_printf(...) fprintf(stderr, __VA_ARGS__)
#define printPts_p printPt(*p0); printPt(*p1); printPt(*p2);
#define printPts printPt(p0); printPt(p1); printPt(p2);
#define dbg_printPt(p) printPt(p)
#else
#define dbg_printf(...)
#define printPts_p
#define printPts
#define dbg_printPt(p)
#endif

#include <stdlib.h>
#define release_assert(e) \
(__builtin_expect(!(e), 0) ? ((void)printf ("%s:%u: failed assertion `%s'\n", __FILE__, __LINE__, #e), abort()) : (void)0)

#ifdef INTERNAL_VERIFICATIONS
#define internal_verification(e) release_assert(e)
#else
#define internal_verification(e)
#endif

#endif
