#include "Energy.h"
#include <cmath>
#include <cstdlib>

EnergyCalculator::EnergyCalculator() {
}

double EnergyCalculator::computeInteractionEnergy(int cell_id_a, int cell_id_b, double alpha) const {
    return interactionStrength(cell_id_a, cell_id_b, alpha);
}

double EnergyCalculator::computeLocalEnergy(int x, int y, 
                                          const std::vector<std::vector<int> >& grid,
                                          const std::vector<Cell>& cells) const {
    double energy = 0.0;
    
    // Interaction energy with neighbors
    std::vector<std::pair<int, int> > neighbors = getNeighbors(x, y);
    int center_cell_id = grid[x][y];
    
    for (size_t i = 0; i < neighbors.size(); ++i) {
        int nx = neighbors[i].first;
        int ny = neighbors[i].second;
        int neighbor_cell_id = grid[nx][ny];
        
        energy += computeInteractionEnergy(center_cell_id, neighbor_cell_id, 
                                         cells[center_cell_id].getAlpha());
    }
    
    // Area constraint energy
    energy += computeAreaConstraintEnergy(cells);
    
    return energy;
}

double EnergyCalculator::computeAreaConstraintEnergy(const std::vector<Cell>& cells) const {
    double energy = 0.0;
    
    for (size_t i = 0; i < cells.size(); ++i) {
        int area_diff = cells[i].getArea() - TARGET_CELL_AREA;
        energy += cells[i].getLambda() * area_diff * area_diff;
    }
    
    return energy;
}

double EnergyCalculator::computeTotalEnergy(const std::vector<std::vector<int> >& grid,
                                          const std::vector<Cell>& cells) const {
    double energy = 0.0;
    
    // Interaction energy
    for (int x = 0; x < GRID_SIZE; ++x) {
        for (int y = 0; y < GRID_SIZE; ++y) {
            std::vector<std::pair<int, int> > neighbors = getNeighbors(x, y);
            int center_cell_id = grid[x][y];
            
            for (size_t i = 0; i < neighbors.size(); ++i) {
                int nx = neighbors[i].first;
                int ny = neighbors[i].second;
                int neighbor_cell_id = grid[nx][ny];
                
                energy += 0.5 * computeInteractionEnergy(center_cell_id, neighbor_cell_id,
                                                       cells[center_cell_id].getAlpha());
            }
        }
    }
    
    // Area constraint energy
    energy += computeAreaConstraintEnergy(cells);
    
    return energy;
}

bool EnergyCalculator::metropolisAccept(double energy_old, double energy_new, double temperature) const {
    double delta_energy = energy_new - energy_old;
    
    if (delta_energy <= 0.0) {
        return true;
    }
    
    double probability = exp(-delta_energy / temperature);
    double random_value = static_cast<double>(rand()) / RAND_MAX;
    
    return random_value < probability;
}

double EnergyCalculator::interactionStrength(int cell_id_a, int cell_id_b, double alpha) {
    if (cell_id_a == cell_id_b) {
        return 0.0;
    }
    return alpha;
}

void EnergyCalculator::applyPeriodicBoundary(int& x, int& y) const {
    if (x >= GRID_SIZE) x = 0;
    if (y >= GRID_SIZE) y = 0;
    if (x < 0) x = GRID_SIZE - 1;
    if (y < 0) y = GRID_SIZE - 1;
}

std::vector<std::pair<int, int> > EnergyCalculator::getNeighbors(int x, int y) const {
    std::vector<std::pair<int, int> > neighbors;
    
    for (int i = 0; i < 8; ++i) {
        int nx = x + NEIGHBOR_DIRECTIONS[i][0];
        int ny = y + NEIGHBOR_DIRECTIONS[i][1];
        
        // Apply periodic boundary conditions
        if (nx >= GRID_SIZE) nx = 0;
        if (ny >= GRID_SIZE) ny = 0;
        if (nx < 0) nx = GRID_SIZE - 1;
        if (ny < 0) ny = GRID_SIZE - 1;
        
        neighbors.push_back(std::make_pair(nx, ny));
    }
    
    return neighbors;
} 