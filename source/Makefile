OBJS = keys.o readkey.o genutil.o mt19937.o inputparams.o dynnetworkutil.o stackutil.o dynamicprocesses.o browndyn.o

FC = gfortran
FLAGS = -O3
DEBUGFLAGS = -ffast-math -funroll-loops
#DEBUGFLAGS = -g -fbacktrace -ffpe-trap=zero,overflow,invalid -fcheck=all -fbounds-check -Wall -Wextra
LIBS =  -lblas -llapack

OSGLIBS =  -lblas -llapack -lgfortran -lm
OSGFLAGS = -static
OSGEXE = ../dynnetwork.exe

# executable name
EXE = ../dynnetwork.exe

.SUFFIXES: .o .f90

# general suffix rules
%.o: %.f90
	$(FC) $(FLAGS) $(DEBUGFLAGS) -c $<

%.o: tools/%.f90
	$(FC) $(FLAGS) $(DEBUGFLAGS) -c $<

default: $(OBJS) main.o
	$(FC) $(FLAGS) $(DEBUGFLAGS) $(OBJS) main.o -o $(EXE) $(LIBS)

osg: $(OBJS) main.o
	$(FC) $(FLAGS) $(OSGFLAGS) $(OBJS) main.o -o $(OSGEXE) $(OSGLIBS)

clean:
	rm *.o *.mod

veryclean: clean
	rm *~ $(EXE)


# ------------
# individual file dependencies
# ------------
chainutil.o: keys.o mt19937.o genutil.o
readkey.o: inputparams.o keys.o genutil.o mt19937.o
genutil.o: mt19937.o
#browndyn.o: mt19937.o chainutil.o
dynnetworkutil.o: stackutil.o genutil.o
dynamicprocesses.o: stackutil.o dynnetworkutil.o
browndyn.o: dynnetworkutil.o dynamicprocesses.o
