#ifndef ENERGY_H
#define ENERGY_H

#include "Cell.h"
#include "Constants.h"
#include <vector>

class EnergyCalculator {
public:
    EnergyCalculator();
    
    // Core energy functions
    double computeInteractionEnergy(int cell_id_a, int cell_id_b, double alpha) const;
    double computeLocalEnergy(int x, int y, const std::vector<std::vector<int> >& grid, 
                              const std::vector<Cell>& cells) const;
    double computeAreaConstraintEnergy(const std::vector<Cell>& cells) const;
    double computeTotalEnergy(const std::vector<std::vector<int> >& grid, 
                              const std::vector<Cell>& cells) const;
    
    // Metropolis algorithm
    bool metropolisAccept(double energy_old, double energy_new, double temperature = 1.0) const;
    
    // Utility functions
    static double interactionStrength(int cell_id_a, int cell_id_b, double alpha);
    
private:
    // Apply periodic boundary conditions
    void applyPeriodicBoundary(int& x, int& y) const;
    
    // Get neighbors with periodic boundary conditions
    std::vector<std::pair<int, int> > getNeighbors(int x, int y) const;
};

#endif // ENERGY_H 