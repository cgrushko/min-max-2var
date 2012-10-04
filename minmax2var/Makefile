uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
uname_M := $(shell sh -c 'uname -m 2>/dev/null || echo not')
UNAME = $(uname_S)_$(uname_M)

# Define these, either here or on the command line to gmake,
# to turn-on external solvers (such as GLPK or MOSEK) in order
# to compare+benchmark them against this solver.
#
EXTERNAL_SOLVER_LPSOLVE55 = 1
EXTERNAL_SOLVER_CGAL = 1
EXTERNAL_SOLVER_SEIDEL = 1
EXTERNAL_SOLVER_GLPK = 1
EXTERNAL_SOLVER_MOSEK = 1

ifdef EXTERNAL_SOLVER_LPSOLVE55
ifeq ($(UNAME), Darwin_x86_64)
LPSOLVE55_PLATFORM = osx64
endif
ifeq ($(UNAME), Linux_x86_64)
LPSOLVE55_PLATFORM = ux64
endif
DEFINES += -DEXTERNAL_SOLVER_LPSOLVE55 -DETA_PFI_LIBRARY=\"./lp_solve_5.5/bfp/bfp_etaPFI/bin/$(LPSOLVE55_PLATFORM)/libbfp_etaPFI.so\"
OBJS += lpsolve55-external-solver.o
LIBS += ./lp_solve_5.5/lpsolve55/bin/$(LPSOLVE55_PLATFORM)/liblpsolve55.a -ldl 
INCLUDE_PATH += -I./lp_solve_5.5/
endif

ifdef EXTERNAL_SOLVER_CGAL
DEFINES += -DEXTERNAL_SOLVER_CGAL
OBJS += cgal-external-solver.o
LIBS += -lCGAL -lgmp -lboost_thread
endif

ifdef EXTERNAL_SOLVER_SEIDEL
DEFINES += -DEXTERNAL_SOLVER_SEIDEL
DEFINES += -DDOUBLE
OBJS += seidel-external-solver.o
LIBS += -lseidel_$(UNAME)
INCLUDE_PATH += -I./seidel-lp-solver/
LIB_PATH += -L./seidel-lp-solver/
endif

ifdef EXTERNAL_SOLVER_GLPK
DEFINES += -DEXTERNAL_SOLVER_GLPK
OBJS += glpk-external-solver.o
LIBS += -lglpk
endif

ifdef EXTERNAL_SOLVER_MOSEK
DEFINES += -DEXTERNAL_SOLVER_MOSEK
OBJS += mosek-external-solver.o
LIBS += -lmosek64
endif

LD=$(CXX)
CXXFLAGS = -O3 -DNDEBUG -g -flto -fuse-linker-plugin -march=native $(INCLUDE_PATH) $(DEFINES)
CFLAGS = $(CXXFLAGS)
LDFLAGS = -O3 -flto -fuse-linker-plugin -g $(LIB_PATH) 
OBJS += solver-main.o solver2d.o solver3d.o performance-data.o predicates.o P3.o   
UNIT_TESTS_OBJS = unittest.o 

unittest: $(OBJS) $(UNIT_TESTS_OBJS)
	$(LD) $(LDFLAGS) $(OBJS) $(UNIT_TESTS_OBJS) $(LIBS)

solver-main.o: solver-main.cpp solver-main.h P3.h my_utilities.h \
  solver2d.h solver3d.h settings.h
solver2d.o: solver2d.cpp solver2d.h P3.h my_utilities.h openmp-reductor.h settings.h
solver3d.o: solver3d.cpp P2.h P3.h solver3d.h my_utilities.h openmp-reductor.h settings.h

unittest.o: unittest.cpp P3.h P2.h solver-main.h \
 solver2d.h my_utilities.h settings.h

performance-data.o: performance-data.cpp settings.h

predicates.o: predicates.c

P3.o: P3.cpp settings.h

cgal-external-solver.o: cgal-external-solver.cpp settings.h cgal-external-solver.h

lpsolve55-external-solver.o: lpsolve55-external-solver.cpp settings.h lpsolve55-external-solver.h

glpk-external-solver.o: glpk-external-solver.cpp settings.h glpk-external-solver.h

seidel-external-solver.o: seidel-external-solver.cpp settings.h seidel-external-solver.h

clean:
	rm -f $(OBJS) $(UNIT_TESTS_OBJS) libhoughsolver.a

.PHONY: clean
