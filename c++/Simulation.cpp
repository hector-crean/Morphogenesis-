#include "Simulation.h"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>

MorphogenesisSimulation::MorphogenesisSimulation(const SimulationParameters& params) 
    : grid_(GRID_SIZE, std::vector<int>(GRID_SIZE, 0)),
      cells_(TOTAL_CELLS),
      current_time_(0),
      params_(params),
      energy_calculator_(),
      io_manager_(params),
      rng_(time(NULL)),
      uniform_dist_(0.0, 1.0),
      grid_dist_(0, GRID_SIZE * GRID_SIZE - 1),
      direction_dist_(0, 7),
      boundary_ratios_(16, 0.0) {
}

MorphogenesisSimulation::~MorphogenesisSimulation() {
}

void MorphogenesisSimulation::run() {
    io_manager_.printSimulationHeader();
    
    // Initialize simulation
    initializeGrid();
    initializeCells();
    
    // Output initial configuration
    outputConfiguration(0);
    
    // Main simulation loop
    for (current_time_ = 1; current_time_ < MAX_MONTE_CARLO_SWEEPS; ++current_time_) {
        performMonteCarloSweep();
        
        // Output every 10 time steps
        if (current_time_ % 10 == 0) {
            outputConfiguration(current_time_);
            outputStatistics(current_time_);
        }
    }
}

void MorphogenesisSimulation::initializeGrid() {
    // Initialize grid with plaid pattern
    initializePlaidPattern();
}

void MorphogenesisSimulation::initializeCells() {
    // Initialize cells with random types based on fraction_B_cells
    for (int i = 0; i < TOTAL_CELLS; ++i) {
        CellType type = (uniform_dist_(rng_) < params_.fraction_B_cells) ? TYPE_B_RECEIVER : TYPE_A_SENDER;
        double alpha = (type == TYPE_B_RECEIVER) ? ALPHA_SOFT : params_.alpha_stiff;
        
        cells_[i] = Cell(i, type, alpha);
        cells_[i].setLambda(LAMBDA_AREA_CONSTRAINT);
    }
    
    // Calculate initial areas
    for (int x = 0; x < GRID_SIZE; ++x) {
        for (int y = 0; y < GRID_SIZE; ++y) {
            cells_[grid_[x][y]].incrementArea();
        }
    }
    
    // Calculate initial center of mass
    calculateCenterOfMass();
}

void MorphogenesisSimulation::initializePlaidPattern() {
    const int cell_size = 4;  // 4x4 cells in plaid pattern
    int cell_id = 0;
    
    for (int block_x = 0; block_x < GRID_SIZE / cell_size; ++block_x) {
        for (int block_y = 0; block_y < GRID_SIZE / cell_size; ++block_y) {
            for (int i = 0; i < cell_size; ++i) {
                for (int j = 0; j < cell_size; ++j) {
                    int x = block_x * cell_size + i;
                    int y = block_y * cell_size + j;
                    grid_[x][y] = cell_id;
                }
            }
            cell_id++;
        }
    }
}

void MorphogenesisSimulation::performMonteCarloSweep() {
    int attempts = GRID_SIZE * GRID_SIZE;
    
    for (int attempt = 0; attempt < attempts; ++attempt) {
        // Pick random site
        int site_index = grid_dist_(rng_);
        int x = site_index % GRID_SIZE;
        int y = site_index / GRID_SIZE;
        
        // Attempt cell move
        attemptCellMove(x, y);
    }
    
    // Update boundary statistics
    updateBoundaryStatistics();
}

bool MorphogenesisSimulation::attemptCellMove(int x, int y) {
    // Get random neighbor
    std::pair<int, int> neighbor = getRandomNeighbor(x, y);
    int nx = neighbor.first;
    int ny = neighbor.second;
    
    // Only attempt move if different cells
    if (grid_[x][y] == grid_[nx][ny]) {
        return false;
    }
    
    // Attempt gene switching first
    attemptGeneSwitch(grid_[x][y], grid_[nx][ny]);
    
    // Calculate energy change
    double energy_old = energy_calculator_.computeLocalEnergy(x, y, grid_, cells_);
    
    // Record old cell and update areas
    int old_cell = grid_[x][y];
    int new_cell = grid_[nx][ny];
    
    cells_[old_cell].decrementArea();
    cells_[new_cell].incrementArea();
    
    // Change cell
    grid_[x][y] = new_cell;
    
    // Calculate new energy
    double energy_new = energy_calculator_.computeLocalEnergy(x, y, grid_, cells_);
    
    // Metropolis test
    if (energy_calculator_.metropolisAccept(energy_old, energy_new)) {
        // Accept move
        return true;
    } else {
        // Reject move - revert changes
        grid_[x][y] = old_cell;
        cells_[old_cell].incrementArea();
        cells_[new_cell].decrementArea();
        return false;
    }
}

bool MorphogenesisSimulation::attemptGeneSwitch(int cell_id, int neighbor_cell_id) {
    Cell& cell = cells_[cell_id];
    Cell& neighbor = cells_[neighbor_cell_id];
    
    // B->C transition when B cell is near A cell
    if (cell.getType() == TYPE_B_RECEIVER && neighbor.getType() == TYPE_A_SENDER) {
        if (cell.attemptGeneSwitch(params_.switch_rate_B_to_C)) {
            updateCellAlpha(cell_id, TYPE_C_STICKY);
            return true;
        }
    }
    
    // C->B transition
    if (cell.getType() == TYPE_C_STICKY) {
        if (cell.attemptGeneSwitch(params_.switch_rate_C_to_B)) {
            updateCellAlpha(cell_id, TYPE_B_RECEIVER);
            return true;
        }
    }
    
    // A->D transition when A cell is near C cell
    if (cell.getType() == TYPE_A_SENDER && neighbor.getType() == TYPE_C_STICKY) {
        if (cell.attemptGeneSwitch(params_.switch_rate_A_to_D)) {
            updateCellAlpha(cell_id, TYPE_D_RESPONSIVE);
            return true;
        }
    }
    
    // D->A transition
    if (cell.getType() == TYPE_D_RESPONSIVE) {
        if (cell.attemptGeneSwitch(params_.switch_rate_D_to_A)) {
            updateCellAlpha(cell_id, TYPE_A_SENDER);
            return true;
        }
    }
    
    return false;
}

void MorphogenesisSimulation::calculateBoundaryInteractions() {
    // Reset boundary counts
    for (int i = 0; i < TOTAL_CELLS; ++i) {
        cells_[i].resetBoundaries();
    }
    
    // Calculate boundary interactions
    for (int x = 0; x < GRID_SIZE; ++x) {
        for (int y = 0; y < GRID_SIZE; ++y) {
            int cell_id = grid_[x][y];
            
            // Check all 8 neighbors
            for (int dir = 0; dir < 8; ++dir) {
                int nx = x + NEIGHBOR_DIRECTIONS[dir][0];
                int ny = y + NEIGHBOR_DIRECTIONS[dir][1];
                
                applyPeriodicBoundary(nx, ny);
                
                int neighbor_id = grid_[nx][ny];
                
                if (cell_id != neighbor_id) {
                    // Determine boundary type based on cell types
                    CellType cell_type = cells_[cell_id].getType();
                    CellType neighbor_type = cells_[neighbor_id].getType();
                    
                    BoundaryType boundary_type = static_cast<BoundaryType>(
                        static_cast<int>(cell_type) * 4 + static_cast<int>(neighbor_type));
                    
                    cells_[cell_id].incrementBoundary(boundary_type);
                }
            }
        }
    }
}

void MorphogenesisSimulation::updateBoundaryStatistics() {
    calculateBoundaryInteractions();
    // Additional boundary statistics processing can be added here
}

void MorphogenesisSimulation::calculateCenterOfMass() {
    // Calculate center of mass for each cell
    for (int i = 0; i < TOTAL_CELLS; ++i) {
        double com_x = 0.0, com_y = 0.0;
        int count = 0;
        
        for (int x = 0; x < GRID_SIZE; ++x) {
            for (int y = 0; y < GRID_SIZE; ++y) {
                if (grid_[x][y] == i) {
                    com_x += x;
                    com_y += y;
                    count++;
                }
            }
        }
        
        if (count > 0) {
            com_x /= count;
            com_y /= count;
            cells_[i].setCOM(com_x, com_y);
            cells_[i].setUnwrappedCOM(com_x, com_y);
        }
    }
}

void MorphogenesisSimulation::handlePeriodicBoundaryConditions() {
    // Handle periodic boundary conditions for center of mass
    // This is a simplified version - the original code had more complex PBC handling
    calculateCenterOfMass();
}

void MorphogenesisSimulation::outputConfiguration(int time_step) {
    clock_t current_time = clock();
    double cpu_time = static_cast<double>(current_time) / CLOCKS_PER_SEC;
    
    io_manager_.printTimeStep(time_step, cpu_time);
    io_manager_.writeConfiguration(time_step, grid_, cells_);
    
    calculateCenterOfMass();
    io_manager_.writeCenterOfMass(time_step, cells_);
}

void MorphogenesisSimulation::outputStatistics(int time_step) {
    io_manager_.printCellCounts(cells_);
    io_manager_.writeStatistics(time_step, cells_);
}

void MorphogenesisSimulation::applyPeriodicBoundary(int& x, int& y) const {
    if (x >= GRID_SIZE) x = 0;
    if (y >= GRID_SIZE) y = 0;
    if (x < 0) x = GRID_SIZE - 1;
    if (y < 0) y = GRID_SIZE - 1;
}

std::pair<int, int> MorphogenesisSimulation::getRandomNeighbor(int x, int y) {
    int direction = direction_dist_(rng_);
    int nx = x + NEIGHBOR_DIRECTIONS[direction][0];
    int ny = y + NEIGHBOR_DIRECTIONS[direction][1];
    
    applyPeriodicBoundary(nx, ny);
    
    return std::make_pair(nx, ny);
}

CellType MorphogenesisSimulation::determineCellType(int cell_id) const {
    return cells_[cell_id].getType();
}

void MorphogenesisSimulation::updateCellAlpha(int cell_id, CellType new_type) {
    cells_[cell_id].setType(new_type);
    
    double alpha;
    switch (new_type) {
        case TYPE_A_SENDER:
            alpha = params_.alpha_stiff;
            break;
        case TYPE_B_RECEIVER:
            alpha = ALPHA_SOFT;
            break;
        case TYPE_C_STICKY:
            alpha = params_.alpha_stiffest;
            break;
        case TYPE_D_RESPONSIVE:
            alpha = params_.alpha_stiffer;
            break;
        default:
            alpha = params_.alpha_stiff;
    }
    
    cells_[cell_id].setAlpha(alpha);
}

void MorphogenesisSimulation::resetBoundaryStatistics() {
    std::fill(boundary_ratios_.begin(), boundary_ratios_.end(), 0.0);
}

void MorphogenesisSimulation::aggregateBoundaryStatistics() {
    // Aggregate boundary statistics from all cells
    // This is a placeholder for more sophisticated boundary analysis
    resetBoundaryStatistics();
} 