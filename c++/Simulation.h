#ifndef SIMULATION_H
#define SIMULATION_H

#include "Cell.h"
#include "Energy.h"
#include "IO.h"
#include "Constants.h"
#include <vector>
#include <random>

class MorphogenesisSimulation {
public:
    MorphogenesisSimulation(const SimulationParameters& params);
    ~MorphogenesisSimulation();
    
    // Main simulation control
    void run();
    
    // Initialization
    void initializeGrid();
    void initializeCells();
    void initializePlaidPattern();
    
    // Monte Carlo steps
    void performMonteCarloSweep();
    bool attemptCellMove(int x, int y);
    
    // Gene switching
    bool attemptGeneSwitch(int cell_id, int neighbor_cell_id);
    
    // Boundary analysis
    void calculateBoundaryInteractions();
    void updateBoundaryStatistics();
    
    // Center of mass calculations
    void calculateCenterOfMass();
    void handlePeriodicBoundaryConditions();
    
    // Statistics and output
    void outputConfiguration(int time_step);
    void outputStatistics(int time_step);
    
private:
    // Simulation state
    std::vector<std::vector<int> > grid_;
    std::vector<Cell> cells_;
    int current_time_;
    
    // Parameters
    SimulationParameters params_;
    
    // Helper classes
    EnergyCalculator energy_calculator_;
    IOManager io_manager_;
    
    // Random number generation
    std::mt19937 rng_;
    std::uniform_real_distribution<double> uniform_dist_;
    std::uniform_int_distribution<int> grid_dist_;
    std::uniform_int_distribution<int> direction_dist_;
    
    // Utility functions
    void applyPeriodicBoundary(int& x, int& y) const;
    std::pair<int, int> getRandomNeighbor(int x, int y);
    CellType determineCellType(int cell_id) const;
    void updateCellAlpha(int cell_id, CellType new_type);
    
    // Statistics tracking
    std::vector<double> boundary_ratios_;
    void resetBoundaryStatistics();
    void aggregateBoundaryStatistics();
};

#endif // SIMULATION_H 