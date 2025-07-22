# Makefile for Morphogenesis Simulation

# Compiler settings
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++98 -O2
LDFLAGS = 

# Directories
SRC_DIR = .
OBJ_DIR = obj
BIN_DIR = bin

# Target executable
TARGET = $(BIN_DIR)/morphogenesis

# Source files
SOURCES = main.cpp Cell.cpp Energy.cpp IO.cpp Simulation.cpp
OBJECTS = $(SOURCES:%.cpp=$(OBJ_DIR)/%.o)

# Header files
HEADERS = Constants.h Cell.h Energy.h IO.h Simulation.h

# Default target
all: $(TARGET)

# Create target executable
$(TARGET): $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(OBJECTS) -o $@ $(LDFLAGS)

# Compile source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS) | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create directories
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Clean build files
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

# Clean and rebuild
rebuild: clean all

# Install target (optional)
install: $(TARGET)
	cp $(TARGET) /usr/local/bin/

# Run target (for testing)
run: $(TARGET)
	./$(TARGET) test_run 1.0 0.3 0.1 0.05 1.2 1.5 0.08 0.04

# Debug build
debug: CXXFLAGS += -g -DDEBUG
debug: $(TARGET)

# Help
help:
	@echo "Available targets:"
	@echo "  all      - Build the executable (default)"
	@echo "  clean    - Remove build files"
	@echo "  rebuild  - Clean and rebuild"
	@echo "  debug    - Build with debug symbols"
	@echo "  run      - Build and run with test parameters"
	@echo "  install  - Install to /usr/local/bin"
	@echo "  help     - Show this help message"

.PHONY: all clean rebuild install run debug help 