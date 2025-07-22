# Morphogenesis Simulation - Reorganized Code Structure

## Overview

This is a Monte Carlo simulation for modeling cell dynamics and morphogenesis in tissue self-assembly. The original single-file implementation has been reorganized into a modular, object-oriented structure for better maintainability and extensibility.

## Physics/Biology Background

The simulation implements a **Cellular Potts Model** for tissue morphogenesis with:

- **Hamiltonian**: `H = Σ J(σ(n),σ(m)) + λ Σ (a_i-A)²`
- **Cell Types**: 
  - Type A (Sender cells)
  - Type B (Receiver cells) 
  - Type C (Sticky cells, activated B cells)
  - Type D (Responsive cells, activated A cells)
- **Gene Switching**: Dynamic cell type transitions based on neighbor interactions
- **Energy Minimization**: Metropolis Monte Carlo algorithm

## Code Organization

### Header Files

- **`Constants.h`**: All simulation constants, enums, and global parameters
- **`Cell.h`**: Cell class definition with properties and behaviors
- **`Energy.h`**: Energy calculation and Metropolis algorithm
- **`IO.h`**: File input/output operations and parameter parsing
- **`Simulation.h`**: Main simulation class and Monte Carlo logic

### Implementation Files

- **`main.cpp`**: Clean main function that orchestrates the simulation
- **`Cell.cpp`**: Cell class implementation
- **`Energy.cpp`**: Energy calculation functions
- **`IO.cpp`**: File I/O operations (to be implemented)
- **`Simulation.cpp`**: Main simulation loop (to be implemented)

## Improvements Made

### 1. **Modular Structure**
- Separated concerns into logical modules
- Each class has a single responsibility
- Clear interfaces between components

### 2. **Better Organization**
- Consolidated all constants in one place
- Used enums for type safety
- Eliminated global variables where possible

### 3. **Object-Oriented Design**
- `Cell` class encapsulates cell properties
- `EnergyCalculator` handles energy computations
- `IOManager` manages file operations
- `MorphogenesisSimulation` orchestrates the simulation

### 4. **Improved Maintainability**
- Clear function and variable names
- Proper encapsulation
- Removed dead/commented code
- Added proper documentation

### 5. **Build System**
- Created Makefile for easy compilation
- Proper dependency management
- Debug and release builds

## Building the Code

```bash
# Build the simulation
make

# Clean build files
make clean

# Build with debug symbols
make debug

# Run with test parameters
make run

# Show help
make help
```

## Usage

```bash
./bin/morphogenesis output_name alpha_stiff fB kon koff alpha_stiffer alpha_stiffest kon2 koff2
```

### Parameters:
1. `output_name`: Name for output files
2. `alpha_stiff`: Cell adherence (stiff cells)
3. `fB`: Fraction of B=receiver cells
4. `kon`: Switch rate ON B→C
5. `koff`: Switch rate OFF C→B
6. `alpha_stiffer`: Cell adherence (stiffer cells)
7. `alpha_stiffest`: Cell adherence (stiffest cells)
8. `kon2`: Switch rate ON A→D
9. `koff2`: Switch rate OFF D→A

## Key Features

### Cell Types and Transitions
- **A (Sender)** → **D (Responsive)** when near C cells
- **B (Receiver)** → **C (Sticky)** when near A cells
- Dynamic switching based on neighbor interactions

### Energy Components
- **Interaction Energy**: J(σ(n),σ(m)) = α if different cells, 0 if same
- **Area Constraint**: λ(a_i - A)² penalizes deviation from target area
- **Metropolis Acceptance**: exp(-ΔE/kT) probability

### Output Files
- **Configuration**: Cell positions and types at each time step
- **Center of Mass**: Tracking cell movement
- **Statistics**: Cell counts and boundary interactions

## Future Enhancements

1. **Complete Implementation**: Finish IO.cpp and Simulation.cpp
2. **Template Support**: Make grid size and cell types configurable
3. **Performance Optimization**: Use more efficient data structures
4. **Visualization**: Add real-time visualization support
5. **Parameter Validation**: Add input validation and error handling
6. **Unit Tests**: Add comprehensive test suite

## Files Structure

```
Morphogenesis/
├── Constants.h          # Global constants and enums
├── Cell.h/.cpp         # Cell class definition and implementation
├── Energy.h/.cpp       # Energy calculation functions
├── IO.h/.cpp          # File input/output operations
├── Simulation.h/.cpp   # Main simulation class
├── main.cpp           # Clean main function
├── Makefile           # Build configuration
├── README_REORGANIZED.md # This documentation
└── gwitch21.cp        # Original monolithic implementation
```

This reorganized structure provides a solid foundation for future development and makes the codebase much more maintainable and extensible. 