CXX = g++
CXXFLAGS = -O2 -Wall -fPIC -std=c++11 -fopenmp
INCLUDES = -ILRModel -Ispline123 -Ilib -I/usr/include/eigen3
PYINCLUDES = $(shell python3 -m pybind11 --includes)
ROOTINCLUDES = -I$(shell root-config --incdir)
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs) -lMinuit2
OBJS123 = bsfit123.o bspline123d.o json11.o profileHist.o
OBJSLRM = compress.o lrfaxial.o lrfaxial3d.o lrfxy.o lrfxyz.o lrfcomp.o lrf.o lrfio.o lrmodel.o transform.o

pymercury: $(OBJS123) $(OBJSLRM) reconstructor.o reconstructor_mp.o wraprec.o
	$(CXX) -shared -o mercury.so $(OBJS123) $(OBJSLRM) reconstructor.o reconstructor_mp.o wraprec.o $(ROOTLIBS) -fopenmp

pylrm: $(OBJS123) $(OBJSLRM) wraplrm.o
	$(CXX) -shared -o lrmodel.so $(OBJS123) $(OBJSLRM) wraplrm.o -fopenmp

test: $(OBJS123) $(OBJSLRM) test.o 
	$(CXX) -o test $(OBJS123) $(OBJSLRM) test.o -fopenmp

json11.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) lib/json11.cpp
	
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

wraprec.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(PYINCLUDES) $(ROOTINCLUDES) wraprec.cpp

test.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) test.cpp
	
clean:
	rm -f *.o *.so

