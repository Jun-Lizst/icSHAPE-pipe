CXX        = g++

CXXFLAGS    = -O3 -std=c++0x -Wall -lPsBL -lhts

all: sam2tab calc_sliding_shape countRT

clean:
	rm *.o || true
	rm sam2tab || true
	rm calc_sliding_shape || true
	rm countRT || true

sam2tab: sam2tab.cpp
	$(CXX) sam2tab.cpp $(CXXFLAGS) -o sam2tab

calc_sliding_shape: calc_sliding_shape.cpp sliding_shape.o
	$(CXX) calc_sliding_shape.cpp sliding_shape.o $(CXXFLAGS) -o calc_sliding_shape

countRT: countRT.cpp sliding_shape.o
	$(CXX) countRT.cpp sliding_shape.o $(CXXFLAGS) -o countRT

sliding_shape.o: sliding_shape.cpp sliding_shape.h
	$(CXX) -c sliding_shape.cpp $(CXXFLAGS) -o sliding_shape.o
