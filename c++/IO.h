#ifndef IO_H
#define IO_H

#include "Cell.h"
#include "Constants.h"
#include <string>
#include <vector>
#include <fstream>

struct SimulationParameters {
    std::string output_name;
    double alpha_stiff;
    double alpha_stiffer;
    double alpha_stiffest;
    double fraction_B_cells;
    double switch_rate_B_to_C;
    double switch_rate_C_to_B;
    double switch_rate_A_to_D;
    double switch_rate_D_to_A;
};

class IOManager {
public:
    IOManager(const SimulationParameters& params);
    ~IOManager();
    
    // Configuration output
    void writeConfiguration(int time_step, const std::vector<std::vector<int> >& grid, 
                           const std::vector<Cell>& cells);
    
    // Center of mass output
    void writeCenterOfMass(int time_step, const std::vector<Cell>& cells);
    
    // Statistics output
    void writeStatistics(int time_step, const std::vector<Cell>& cells);
    
    // Boundary analysis output
    void writeBoundaryAnalysis(int time_step, const std::vector<Cell>& cells);
    
    // Parse command line arguments
    static SimulationParameters parseCommandLine(int argc, char* argv[]);
    
    // Print simulation info
    void printSimulationHeader() const;
    void printTimeStep(int time_step, double cpu_time) const;
    void printCellCounts(const std::vector<Cell>& cells) const;
    
private:
    SimulationParameters params_;
    std::ofstream config_file_;
    std::ofstream com_file_;
    std::ofstream stats_file_;
    
    // Helper functions
    std::string generateFileName(const std::string& prefix, int time_step = -1) const;
    void createOutputDirectory() const;
};

#endif // IO_H 