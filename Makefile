CXX = g++
CXXFLAGS = -fopenmp -Wall -Wextra -I./libs

# Все .cpp файлы из папки main
SRCS = $(wildcard main/*.cpp)
TARGET = build/my_program

all: $(TARGET)

$(TARGET): $(SRCS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -rf build

.PHONY: all run clean
