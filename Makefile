# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++17 -Wall

# Detect OS (macOS or Linux)
UNAME_S := $(shell uname -s)

# Set OS-specific flags
ifeq ($(UNAME_S), Darwin)  # macOS
    INCLUDES = -I/opt/homebrew/include
    LIBS = -L/opt/homebrew/lib -lhdf5 -lhdf5_cpp
else  # Linux (Ubuntu/Debian)
    INCLUDES = -I/usr/include/hdf5/serial -I/usr/local/include
    LIBS = -L/usr/lib/x86_64-linux-gnu -lhdf5 -lhdf5_cpp
endif

# Source and output files
SRC = collapse_hdf5.cpp
OBJ = $(SRC:.cpp=.o)
TARGET = collapse_hdf5

# Default build target
all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $(TARGET) $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)