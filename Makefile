CXX = g++
CXXFLAGS = -O2 -Wall -fPIC -std=c++11 -fopenmp
#CXXFLAGS = -g -Wall -fPIC -std=c++11 -fopenmp
INCLUDES = -ILRModel -Ispline123 -Ilib -I/usr/include/eigen3
PYINCLUDES = $(shell python3 -m pybind11 --includes)

#ROOTINCLUDES = -I$(shell root-config --incdir)
#ROOTFLAGS = $(shell root-config --cflags)
#ROOTLIBS = $(shell root-config --libs) -lMinuit2
MINUITDIR = ../Minuit2
ROOTINCLUDES = -I$(MINUITDIR)/inc
ROOTLIBS = -L$(MINUITDIR)/build/lib -lMinuit2 -lMinuit2Math
#ROOTINCLUDES = -I../Minuit2/inc
#ROOTLIBS = ../Minuit2/build/lib/libMinuit2.a ../Minuit2/build/lib/libMinuit2Math.a

OBJS123 = bsfit123.o bspline123d.o json11.o profileHist.o wformula.o
OBJSLRM = compress.o lrfaxial.o lrformula1.o lrformulav.o lrformulaxy.o lrfaxial3d.o lrfxy.o lrfxyz.o lrfcomp.o lrf.o lrfio.o lrmodel.o transform.o

libmercury: $(OBJS123) $(OBJSLRM) reconstructor.o reconstructor_mp.o
	$(CXX) -shared -o libmercury.so $(OBJS123) $(OBJSLRM) reconstructor.o reconstructor_mp.o $(ROOTLIBS) -fopenmp

liblrm: $(OBJS123) $(OBJSLRM)
	$(CXX) -shared -o liblrm.so $(OBJS123) $(OBJSLRM) -fopenmp

pymercury: $(OBJS123) $(OBJSLRM) reconstructor.o reconstructor_mp.o wraprec.o
	$(CXX) -shared -o mercury.so $(OBJS123) $(OBJSLRM) reconstructor.o reconstructor_mp.o wraprec.o $(ROOTLIBS) -fopenmp

pylrm: $(OBJS123) $(OBJSLRM) wraplrm.o
	$(CXX) -shared -o lrmodel.so $(OBJS123) $(OBJSLRM) wraplrm.o -fopenmp

pyaxial: $(OBJS123) $(OBJSLRM) wrapaxial.o
	$(CXX) -shared -o lrfaxial.so $(OBJS123) $(OBJSLRM) wrapaxial.o -fopenmp

pyformula1: $(OBJS123) $(OBJSLRM) wrapformula1.o
	$(CXX) -shared -o lrformula1.so $(OBJS123) $(OBJSLRM) wrapformula1.o -fopenmp

pyformulav: $(OBJS123) $(OBJSLRM) wrapformulav.o
	$(CXX) -shared -o lrformulav.so $(OBJS123) $(OBJSLRM) wrapformulav.o -fopenmp

pyformulaxy: $(OBJS123) $(OBJSLRM) wrapformulaxy.o
	$(CXX) -shared -o lrformulaxy.so $(OBJS123) $(OBJSLRM) wrapformulaxy.o -fopenmp

test: $(OBJS123) $(OBJSLRM) test.o 
	$(CXX) -o test $(OBJS123) $(OBJSLRM) test.o -fopenmp

testxy: $(OBJS123) $(OBJSLRM) testxy.o 
	$(CXX) -o testxy $(OBJS123) $(OBJSLRM) testxy.o -fopenmp

json11.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) lib/json11.cpp

wformula.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) lib/wformula.cpp
	
profileHist.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) spline123/profileHist.cpp

bspline123d.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) spline123/bspline123d.cpp
	
bsfit123.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) spline123/bsfit123.cpp

compress.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) LRModel/compress.cpp

transform.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) LRModel/transform.cpp

lrf.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) LRModel/lrf.cpp

lrfaxial.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) LRModel/lrfaxial.cpp

lrformula1.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) LRModel/lrformula1.cpp

lrformulav.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) LRModel/lrformulav.cpp

lrformulaxy.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) LRModel/lrformulaxy.cpp

lrfaxial3d.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) LRModel/lrfaxial3d.cpp

lrfxy.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) LRModel/lrfxy.cpp

lrfxyz.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) LRModel/lrfxyz.cpp	

lrfcomp.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) LRModel/lrfcomp.cpp

lrmodel.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) LRModel/lrmodel.cpp		
	
lrfio.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) LRModel/lrfio.cpp

reconstructor.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(ROOTINCLUDES) reconstructor.cpp	

reconstructor_mp.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(ROOTINCLUDES) reconstructor_mp.cpp		

wraplrm.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(PYINCLUDES) LRModel/wraplrm.cpp

wrapaxial.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(PYINCLUDES) LRModel/wrapaxial.cpp

wrapformula1.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(PYINCLUDES) LRModel/wrapformula1.cpp

wrapformulav.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(PYINCLUDES) LRModel/wrapformulav.cpp

wrapformulaxy.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(PYINCLUDES) LRModel/wrapformulaxy.cpp

wraprec.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(PYINCLUDES) $(ROOTINCLUDES) wraprec.cpp

test.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) tests/test.cpp

testxy.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) tests/testxy.cpp
	
clean:
	rm -f *.o *.so


