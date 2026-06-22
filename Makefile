CXX = g++
CXXFLAGS = -fopenmp -Wall -Wextra -O3 -march=native

CUDA_CXX = nvcc
CUDA_CXXFLAGS = -Xcompiler="-fopenmp -Wall -Wextra -O3 -march=native" -DUSE_CUDA -x cu -arch=sm_70 -std=c++14

SRCS = main/main.cpp
TARGET = build/md_simulation
TARGET_CUDA = build/md_simulation_cuda

all: $(TARGET)

$(TARGET): $(SRCS) libs/system.h
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET)

cuda: $(SRCS) system.h
	$(CUDA_CXX) $(CUDA_CXXFLAGS) $(SRCS) -o $(TARGET_CUDA)

clean:
	rm -f $(TARGET) $(TARGET_CUDA) trajectory.csv

.PHONY: all cuda clean