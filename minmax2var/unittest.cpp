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
//  Created by Carmi Grushko on 6/18/12.

#include <stdio.h>
#include <stdlib.h>
#include "P3.h"
#include "P2.h"

#include "cgal-external-solver.h"
#include "lpsolve55-external-solver.h"
#include "glpk-external-solver.h"
#include "seidel-external-solver.h"
#include "mosek-external-solver.h"

#include <algorithm>
#include <cmath>
#include <sys/time.h>
#include <vector>
#include <fstream>
#include <tr1/memory>
#include <sstream>
#include <map>
#include <set>
#ifdef __WIN32__
#include <Winsock2.h>
#else
#include <arpa/inet.h>
#endif

#include "solver-main.h"
#include "solver2d.h"
#include "my_utilities.h"
#include "settings.h"
#include "performance-data.h"

using std::vector;
using std::tr1::shared_ptr;

using hough_min_max::P3;
using hough_min_max::V3;
using hough_min_max::performanceDataSingleton;

int participating_external_solver = -1;

void assertPrimalSolution(P3 ptArray[], int n, double u, double v, double val) {
  double static const h = 1e-9;
  for (int i=0; i<n; i++) {
    double z = ptArray[i].x * u + ptArray[i].y * v - ptArray[i].z;
    release_assert(z - h <= val);
  }
}

template <int Coordinate>
void assertPrimalSolution2D(P3 ptArray[], int n, double u, double val) {
  double static const h = 1e-9;
  for (int i=0; i<n; i++) {
    double z = ptArray[i].coord[Coordinate] * u - ptArray[i].z;
    release_assert(z - h <= val);
  }
}

shared_ptr<ExternalSolver> create_external_solver(int idx) {
  switch (idx) {
    case 0: return shared_ptr<ExternalSolver>(new lpsolve55_Solver());
    case 1: return shared_ptr<ExternalSolver>(new lpsolve55_Solver_etaPFI());
    case 2: return shared_ptr<ExternalSolver>(new CGAL_Solver());
    case 3: return shared_ptr<ExternalSolver>(new GLPK_Solver_primal<false>());
    case 4: return shared_ptr<ExternalSolver>(new GLPK_Solver_dual<false>());
    case 5: return shared_ptr<ExternalSolver>(new GLPK_Solver_primal<true>());
    case 6: return shared_ptr<ExternalSolver>(new GLPK_Solver_dual<true>());
    case 7: return shared_ptr<ExternalSolver>(new Seidel_Solver());
    case 8: return shared_ptr<ExternalSolver>(new MOSEK_Solver(MOSEK_Solver_Type::IPM));
    case 9: return shared_ptr<ExternalSolver>(new MOSEK_Solver(MOSEK_Solver_Type::PrimalSimplex));
    case 10: return shared_ptr<ExternalSolver>(new MOSEK_Solver(MOSEK_Solver_Type::DualSimplex));
  }
  return shared_ptr<ExternalSolver>();
}

vector<shared_ptr<ExternalSolver> > create_external_solvers_vector() {
  if (participating_external_solver == -2)
    return vector<shared_ptr<ExternalSolver> >();
  if (participating_external_solver == -1) {
    vector<shared_ptr<ExternalSolver> > result;
    shared_ptr<ExternalSolver> external_solver;
    int i = 0;
    while (1) {
      external_solver = create_external_solver(i++);
      if (external_solver == NULL)
        break;
      result.push_back(external_solver);
    }
    return result;
  } else {
    return vector<shared_ptr<ExternalSolver> >(1, create_external_solver(participating_external_solver));
  }
}

bool load_data_from_file(const char* filename, std::vector<P3>* points, bool read_in_hex = false) {
  if (!read_in_hex) {
    FILE* f = fopen(filename, "r");
    if (f == NULL)
      return false;
    while(!feof(f)) {
      double x, y, z;
      if (fscanf(f, "%lf %lf %lf", &x, &y, &z) != 3)
        abort();
      points->push_back(P3(x, y, z));
    }
  } else {
    std::ifstream f(filename, std::ios::in);
    if (!f)
      return false;
    while (f) {
      points->push_back(hough_min_max::readPt_hex(f));
    }
  }
  return true;
}

typedef std::map<std::string, std::set<int> > ProblemSkipMap;

ProblemSkipMap load_skip_map(const char* filename) {
  ProblemSkipMap result;
  std::ifstream f(filename, std::ios::in);
  std::string line, word;
  std::getline(f, line);
  while (f) {
    std::stringstream ss(line);
    ss >> word;
    if (word != "#") {
      int problem_number;
      ss >> problem_number;
      while (ss) {
        result[word].insert(problem_number);
        ss >> problem_number;
      }
    }
    std::getline(f, line);
  }
  return result;
}

bool is_big_endian() {
  return htonl(43) == 43;
}

bool test3d_by_loading_from_file(const char* filename, bool read_in_hex) {
  std::vector<P3> points_vec;
  bool loaded = load_data_from_file(filename, &points_vec, read_in_hex);
  if (!loaded)
    return false;
  P3* points = points_vec.data();
  int k = (int)points_vec.size();
  
  vector<shared_ptr<ExternalSolver> > external_solvers = create_external_solvers_vector();
  for (int i=0; i<external_solvers.size(); i++)
    external_solvers[i]->init(k);
  
  for (int i=0; i<k; i++) {
    dbg_printPt(points[i]);
  }

  double hough_u, hough_v, hough_value;
  std::vector<P3> points_copy(points, points+k);
  hough_min_max::solve(points, k, &hough_u, &hough_v, &hough_value);
  printf("hough:   U, V, Value: %lf %lf %.7lf\n", hough_u, hough_v, hough_value);
  assertPrimalSolution(points_copy.data(), k, hough_u, hough_v, hough_value);
  
  for (int i=0; i<external_solvers.size(); i++) {
    external_solvers[i]->solve(points_copy.data());
    external_solvers[i]->print_status(std::cout);
    if (!external_solvers[i]->assert_solution(true, hough_u, hough_v, hough_value)) {
      hough_min_max::print_problem_hex(points_copy.data(), k, std::cout);
      abort();
    }
  }
  return true;
}

void test1(std::string directory) {
  int i = 1;
  while (1) {
    std::stringstream ss;
    ss << directory << i << "test-in.txt";
    if (!test3d_by_loading_from_file(ss.str().c_str(), true))
      break;
    i++;
  }
}


// Test whole system
void test4(int k, int start_from_problem_number = 0, int reporting_step = 100) {
  timeval t1, t2;
  double hough_time = 0;
  const static int progress_indiciators_in_line = 20;
  
  P3* points = NULL;
  if (posix_memalign((void**)&points, 16, sizeof(P3)*k) != 0) {
    printf("Error in posix_memalign\n");
    abort();
  }

  vector<shared_ptr<ExternalSolver> > external_solvers = create_external_solvers_vector();
  for (int i=0; i<external_solvers.size(); i++)
    external_solvers[i]->init(k);

  for (int j=0; j<start_from_problem_number; j++) {
    for (int i=0; i<k; i++) {
      points[i].x = (rand() / (double)RAND_MAX) - 0.5;
      points[i].y = (rand() / (double)RAND_MAX) - 0.5;
      points[i].z = (rand() / (double)RAND_MAX) - 0.5;
    }
  }
  int printed_number = 1;
  for (int j=0; j<10000; j++) {
    if (j % reporting_step == 0) {
      if (printed_number % progress_indiciators_in_line == 1)
        printf("Problem number: ");
      printf("%d ", j);
      if (printed_number % progress_indiciators_in_line == 0)
        printf("\n");
      printed_number++;
    }

    for (int i=0; i<k; i++) {
      points[i].x = (rand() / (double)RAND_MAX) - 0.5;
      points[i].y = (rand() / (double)RAND_MAX) - 0.5;
      points[i].z = ((rand() / (double)RAND_MAX) - 0.5)*20;
      dbg_printPt(points[i]);
    }

    double hough_u, hough_v, hough_value;
    std::vector<P3> points_copy(points, points+k);
    gettimeofday(&t1, NULL);
    hough_min_max::solve(points, k, &hough_u, &hough_v, &hough_value);
    gettimeofday(&t2, NULL);    
    hough_time += (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
    hough_time += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms    
//    printf("hough:   U, V, Value: %lf %lf %lf\n", hough_u, hough_v, hough_value);
    //    assertPrimalSolution(points, k, u, v, value);
    //      assertSolution(points, k, p0, p1, p2);
    
    for (int i=0; i<external_solvers.size(); i++) {
      external_solvers[i]->solve(points_copy.data());
//      external_solvers[i]->print_status(std::cout);
      if (!external_solvers[i]->assert_solution(true, hough_u, hough_v, hough_value)) {
        hough_min_max::print_problem_hex(points_copy.data(), k, std::cout);
        abort();
      }
    }

  }
  printf("Hough total time (ms): %lf\n", hough_time);

  for (int i=0; i<external_solvers.size(); i++) {
    double time = external_solvers[i]->time;
    std::string name = external_solvers[i]->name();
    printf("%s total time (ms): %lf\n", name.c_str(), time);
    printf("  Speed-up: %lf\n", time / hough_time);
  }
  
  free(points);
}

void test_lps_from_GMDS(const char* filename, std::set<int> problem_skip_list, int start_from_problem_number = 0, int reporting_step = 100) {
  release_assert(!is_big_endian());
  
  timeval t1, t2;
  double hough_time = 0;
  const static int progress_indiciators_in_line = 20;

  vector<shared_ptr<ExternalSolver> > external_solvers = create_external_solvers_vector();
  vector<P3> points_vector;
  std::ifstream file(filename);
  int problem_counter = 0;
  int printed_number = 1;
  while (1) {
    problem_counter++;
    if (problem_counter >= start_from_problem_number) {
      if (problem_counter % reporting_step == 0) {
        if (printed_number % progress_indiciators_in_line == 1)
          printf("Problem number: ");
        printf("%d ", problem_counter);
        if (printed_number % progress_indiciators_in_line == 0)
          printf("\n");
        printed_number++;
      }
    }

    int k;
    char dummy;
    file >> k;
    if (!file)
      break;
    file.get(dummy);
    release_assert(dummy = ' ');
    
    points_vector.resize(k);
    P3* points = points_vector.data();
    
    for (int i=0; i<k; i++)
      points[i] = hough_min_max::readPt_bin(file);
    file.get(dummy);
    release_assert(dummy == '\n');

    if (problem_counter < start_from_problem_number)
      continue;
    if (problem_skip_list.find(problem_counter) != problem_skip_list.end())
      continue;

    for (int i=0; i<external_solvers.size(); i++)
      external_solvers[i]->init(k);
    
    double hough_u, hough_v, hough_value;
    std::vector<P3> points_copy(points, points+k);
    gettimeofday(&t1, NULL);
    hough_min_max::solve(points, k, &hough_u, &hough_v, &hough_value);
    gettimeofday(&t2, NULL);
    hough_time += (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
    hough_time += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
//    printf("hough:   U, V, Value: %lf %lf %lf\n", hough_u, hough_v, hough_value);
    //    assertPrimalSolution(points, k, u, v, value);
    //      assertSolution(points, k, p0, p1, p2);
    
    for (int i=0; i<external_solvers.size(); i++) {
      external_solvers[i]->solve(points_copy.data());
//      external_solvers[i]->print_status(std::cout);
      if (!external_solvers[i]->assert_solution(true, hough_u, hough_v, hough_value)) {
        fprintf(stderr, "Solution differs with %s. on problem %d\n", external_solvers[i]->name().c_str(), problem_counter);
        external_solvers[i]->print_status(std::cout);
        printf("hough:   U, V, Value: %lf %lf %lf\n", hough_u, hough_v, hough_value);
        //hough_min_max::print_problem_hex(points_copy.data(), k, std::cout);
        abort();
      }
    }
  }
  if (problem_counter % progress_indiciators_in_line != 0)
    printf("\n");
  printf("Hough total time (ms): %lf\n", hough_time);
  
  for (int i=0; i<external_solvers.size(); i++) {
    if (external_solvers[i]->name().empty())
      continue;
    double time = external_solvers[i]->time;
    std::string name = external_solvers[i]->name();
    printf("%s total time (ms): %lf\n", name.c_str(), time);
    printf("  Speed-up: %lf\n", time / hough_time);
  }
}


// Test 2D solvers
void test5() {
  
  int k = 4452;
  P3* points = NULL;
  if (posix_memalign((void**)&points, 16, sizeof(P3)*k) != 0) {
    printf("Error in posix_memalign\n");
    abort();
  }
#ifdef LPSOLVE
  lprec* rec_x = setup_unconstrained_lp(k);
  lprec* rec_y = setup_unconstrained_lp(k);
#endif
  
  for (int j=0; j<0; j++) { 
    for (int i=0; i<k; i++) {
      points[i].x = (rand() / (double)RAND_MAX) - 0.5;
      points[i].y = (rand() / (double)RAND_MAX) - 0.5;
      points[i].z = (rand() / (double)RAND_MAX) - 0.5;
    }
  }
  for (int j=0; j<100000; j++) { 
    if (j % 1 == 0)
      printf("%d\n", j);
    for (int i=0; i<k; i++) {
      points[i].x = (rand() / (double)RAND_MAX) - 0.5;
      points[i].y = (rand() / (double)RAND_MAX) - 0.5;
      points[i].z = ((rand() / (double)RAND_MAX) - 0.5)*20;
      dbg_printPt(points[i]);
#ifdef LPSOLVE
      setConstraint(rec_x, i - 1, points[i].x, 0, points[i].z);
      setConstraint(rec_y, i - 1, 0, points[i].y, points[i].z);
#endif
    }
    
#ifdef LPSOLVE
//    print_lp(rec_x);
//    print_lp(rec_y);
    bool bounded_x = solve(rec_x) == 0;
    bool bounded_y = solve(rec_y) == 0;
    REAL var_x[3];
    REAL var_y[3];
    get_variables(rec_x, var_x);
    get_variables(rec_y, var_y);
    printf("lpsolve_x: U, V, Value: %lf %lf %lf\n", var_x[0], var_x[1], var_x[2]);
    printf("lpsolve_y: U, V, Value: %lf %lf %lf\n", var_y[0], var_y[1], var_y[2]);
#endif
    hough_min_max::Result result_x1, result_x2, result_y;
    hough_min_max::solver_2d::solve_on_edges(points, k, &result_x1, &result_y);
    hough_min_max::solver_2d::solve_on_x(points, k, &result_x2);
    printf("hough_x1:   U, Value: %lf %lf\n", result_x1.u, result_x1.value);
    printf("hough_x2:   U, Value: %lf %lf\n", result_x2.u, result_x2.value);
    printf("hough_y:    V, Value: %lf %lf\n", result_y.v, result_y.value);
    //    assertPrimalSolution(points, k, u, v, value);
    //      assertSolution(points, k, p0, p1, p2);
#ifdef LPSOLVE
    release_assert(!bounded_x == !std::isfinite(result_x1.u));
    release_assert(!bounded_x == !std::isfinite(result_x2.u));
    release_assert(!bounded_y == !std::isfinite(result_y.v));
    
    static const double h = 1e-5;
    if (bounded_x) { 
      release_assert(fabs(result_x1.u-var_x[0]) < h);
      release_assert(fabs(result_x1.value-var_x[2]) < h);
      release_assert(fabs(result_x2.u-var_x[0]) < h);
      release_assert(fabs(result_x2.value-var_x[2]) < h);
    }

    if (bounded_y) {
      release_assert(fabs(result_y.v-var_y[1]) < h);
      release_assert(fabs(result_y.value-var_y[2]) < h);
    }
#endif
  }
  
#ifdef LPSOLVE
  delete_lp(rec_x);
  delete_lp(rec_y);
#endif
  
  free(points);
}

void test2d_by_loading_from_file() {
  
  std::vector<P3> points_vec;
  load_data_from_file("/Users/n/temp/temp.txt", &points_vec);
  P3* points = points_vec.data();
  int k = (int)points_vec.size();
  
#ifdef LPSOLVE
  lprec* rec_x = setup_unconstrained_lp(k);
  lprec* rec_y = setup_unconstrained_lp(k);
#endif

  for (int i=0; i<k; i++) {
    dbg_printPt(points[i]);
#ifdef LPSOLVE
    setConstraint(rec_x, i - 1, points[i].x, 0, points[i].z);
    setConstraint(rec_y, i - 1, 0, points[i].y, points[i].z);
#endif
  }
  
#ifdef LPSOLVE
  //    print_lp(rec_x);
  //    print_lp(rec_y);
  bool bounded_x = solve(rec_x) == 0;
  bool bounded_y = solve(rec_y) == 0;
  REAL var_x[3];
  REAL var_y[3];
  get_variables(rec_x, var_x);
  get_variables(rec_y, var_y);
  printf("lpsolve_x: U, V, Value: %lf %lf %lf\n", var_x[0], var_x[1], var_x[2]);
  printf("lpsolve_y: U, V, Value: %lf %lf %lf\n", var_y[0], var_y[1], var_y[2]);
#endif
  hough_min_max::Result result_x1, result_x2, result_y;
  hough_min_max::solver_2d::solve_on_edges(points, k, &result_x1, &result_y);
  hough_min_max::solver_2d::solve_on_x(points, k, &result_x2);
  printf("hough_x1:   U, Value: %lf %lf\n", result_x1.u, result_x1.value);
  printf("hough_x2:   U, Value: %lf %lf\n", result_x2.u, result_x2.value);
  printf("hough_y:    V, Value: %lf %lf\n", result_y.v, result_y.value);
  assertPrimalSolution2D<P3::X>(points, k, result_x2.u, result_x2.value);
  assertPrimalSolution2D<P3::X>(points, k, result_x1.u, result_x1.value);
  assertPrimalSolution2D<P3::Y>(points, k, result_y.v, result_y.value);
  //      assertSolution(points, k, p0, p1, p2);
#ifdef LPSOLVE
  release_assert(!bounded_x == !std::isfinite(result_x1.u));
  release_assert(!bounded_x == !std::isfinite(result_x2.u));
  release_assert(!bounded_y == !std::isfinite(result_y.v));
  
  static const double h = 1e-5;
  if (bounded_x) {
    release_assert(fabs(result_x1.u-var_x[0]) < h);
    release_assert(fabs(result_x1.value-var_x[2]) < h);
    release_assert(fabs(result_x2.u-var_x[0]) < h);
    release_assert(fabs(result_x2.value-var_x[2]) < h);
  }
  
  if (bounded_y) {
    release_assert(fabs(result_y.v-var_y[1]) < h);
    release_assert(fabs(result_y.value-var_y[2]) < h);
  }
#endif

#ifdef LPSOLVE
  delete_lp(rec_x);
  delete_lp(rec_y);
#endif

}

void benchmark1() {
  timeval t1, t2;
  double hough_time = 0;
  
  int k = 4452;
  P3* points;// = new P3[k];
  if (posix_memalign((void**)&points, 16, sizeof(P3)*k) != 0) {
    printf("Error in posix_memalign\n");
    abort();
  }
  
  for (int j=0; j<10000; j++) {
    if (j % 1000 == 0)
      printf("%d\n", j);
    for (int i=0; i<k; i++) {
      points[i].x = (rand() / (double)RAND_MAX) - 0.5;
      points[i].y = (rand() / (double)RAND_MAX) - 0.5;
      points[i].z = ((rand() / (double)RAND_MAX) - 0.5)*20;
    }
    
    double u, v, value;
    std::vector<P3> points_copy(points, points+k);
    gettimeofday(&t1, NULL);
    hough_min_max::solve(points, k, &u, &v, &value);
    gettimeofday(&t2, NULL);
    hough_time += (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
    hough_time += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
    
  }
  printf("Hough total time (ms): %lf\n", hough_time);
  
  free(points);
}

/**
 * Syntax: ./a.out [number of solver] [skip-map]
 * number of external solver: integer, or -1 (=nothing) for all of them, or -2 for none
 */
int main(int argc, const char * argv[])
{
  srand(10);
  setvbuf(stdout, NULL, _IONBF, 0);
  ProblemSkipMap problem_skip_map;
  if (argc > 1)
    participating_external_solver = atoi(argv[1]);
  if (argc > 2)
    problem_skip_map = load_skip_map(argv[2]);

  hough_min_max::init();
  
  printf("Random 1000\n"); test4(1000);
  performanceDataSingleton.start_new_measurement();
  printf("Random 2000\n"); test4(2000);
  performanceDataSingleton.start_new_measurement();
  printf("Random 4000\n"); test4(4000);
  performanceDataSingleton.start_new_measurement();
  printf("Random 8000\n"); test4(8000);

  performanceDataSingleton.start_new_measurement();
  printf("GMDS 1110\n"); test_lps_from_GMDS("GMDS-data/1110.dat", problem_skip_map["GMDS-1110"]);
  performanceDataSingleton.start_new_measurement();
  printf("GMDS 2224\n"); test_lps_from_GMDS("GMDS-data/2224.dat", problem_skip_map["GMDS-2224"]);
  performanceDataSingleton.start_new_measurement();
  printf("GMDS 4450\n"); test_lps_from_GMDS("GMDS-data/4450.dat", problem_skip_map["GMDS-4450"]);
  performanceDataSingleton.start_new_measurement();
  printf("GMDS 8902\n"); test_lps_from_GMDS("GMDS-data/8902.dat", problem_skip_map["GMDS-8902"]);

  printf("Done\n");
  return 0;
}

