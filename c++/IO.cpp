#include "IO.h"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cstdlib>
#include <sys/stat.h>
#include <sys/types.h>

IOManager::IOManager(const SimulationParameters& params) : params_(params) {
    createOutputDirectory();
    
    // Open initial files
    std::string config_filename = generateFileName(params_.output_name + "_t", 0);
    config_file_.open(config_filename.c_str());
    if (!config_file_) {
        throw std::runtime_error("Cannot open configuration file: " + config_filename);
    }
    
    std::string com_filename = generateFileName("com_" + params_.output_name);
    com_file_.open(com_filename.c_str());
    if (!com_file_) {
        throw std::runtime_error("Cannot open COM file: " + com_filename);
    }
}

IOManager::~IOManager() {
    if (config_file_.is_open()) {
        config_file_.close();
    }
    if (com_file_.is_open()) {
        com_file_.close();
    }
    if (stats_file_.is_open()) {
        stats_file_.close();
    }
}

void IOManager::writeConfiguration(int time_step, const std::vector<std::vector<int> >& grid, 
                                  const std::vector<Cell>& cells) {
    if (time_step > 0) {
        config_file_.close();
        std::string filename = generateFileName(params_.output_name + "_t", time_step);
        config_file_.open(filename.c_str());
        if (!config_file_) {
            throw std::runtime_error("Cannot open configuration file: " + filename);
        }
    }
    
    for (int x = 0; x < GRID_SIZE; ++x) {
        for (int y = 0; y < GRID_SIZE; ++y) {
            int cell_id = grid[x][y];
            config_file_ << (x + 0.5) << " " << (y + 0.5) << " " << cell_id << " "
                        << cells[cell_id].getAlpha() << " " << cells[cell_id].getType() << std::endl;
        }
    }
}

void IOManager::writeCenterOfMass(int time_step, const std::vector<Cell>& cells) {
    com_file_ << "#t" << time_step << std::endl;
    
    for (size_t i = 0; i < cells.size(); ++i) {
        com_file_ << i << " " << cells[i].getUnwrappedCOMX() + 0.5 << " " 
                  << cells[i].getUnwrappedCOMY() + 0.5 << std::endl;
    }
    
    com_file_ << std::endl;
    com_file_ << std::endl;
}

void IOManager::writeStatistics(int time_step, const std::vector<Cell>& cells) {
    // This can be extended to write various statistics
    // For now, just a placeholder
    (void)time_step;
    (void)cells;
}

void IOManager::writeBoundaryAnalysis(int time_step, const std::vector<Cell>& cells) {
    // This can be extended to write boundary analysis
    // For now, just a placeholder
    (void)time_step;
    (void)cells;
}

SimulationParameters IOManager::parseCommandLine(int argc, char* argv[]) {
    if (argc != 10) {
        throw std::runtime_error("Usage: ./morphogenesis output_name alpha_stiff fB kon koff alpha_stiffer alpha_stiffest kon2 koff2");
    }
    
    SimulationParameters params;
    
    try {
        params.output_name = std::string(argv[1]);
        params.alpha_stiff = atof(argv[2]);
        params.fraction_B_cells = atof(argv[3]);
        params.switch_rate_B_to_C = atof(argv[4]);
        params.switch_rate_C_to_B = atof(argv[5]);
        params.alpha_stiffer = atof(argv[6]);
        params.alpha_stiffest = atof(argv[7]);
        params.switch_rate_A_to_D = atof(argv[8]);
        params.switch_rate_D_to_A = atof(argv[9]);
    } catch (const std::exception& e) {
        throw std::runtime_error("Error parsing command line arguments: " + std::string(e.what()));
    }
    
    return params;
}

void IOManager::printSimulationHeader() const {
    std::cout << "## ATTENTION! HAVE YOU CREATED THE FOLDER IN WHICH I SHOULD WRITE??? ##" << std::endl;
    std::cout << "Parameters:" << std::endl;
    std::cout << "  Output name: " << params_.output_name << std::endl;
    std::cout << "  Alpha stiff: " << params_.alpha_stiff << std::endl;
    std::cout << "  Fraction B cells: " << params_.fraction_B_cells << std::endl;
    std::cout << "  Switch rates: " << params_.switch_rate_B_to_C << " (B->C), " 
              << params_.switch_rate_C_to_B << " (C->B)" << std::endl;
    std::cout << "  Switch rates: " << params_.switch_rate_A_to_D << " (A->D), " 
              << params_.switch_rate_D_to_A << " (D->A)" << std::endl;
}

void IOManager::printTimeStep(int time_step, double cpu_time) const {
    std::cout << "############## time " << time_step << " #################" << std::endl;
    std::cout << "CPU Time = " << cpu_time << " seconds" << std::endl;
}

void IOManager::printCellCounts(const std::vector<Cell>& cells) const {
    int nA = 0, nB = 0, nC = 0, nD = 0;
    
    for (size_t i = 0; i < cells.size(); ++i) {
        switch (cells[i].getType()) {
            case TYPE_A_SENDER: nA++; break;
            case TYPE_B_RECEIVER: nB++; break;
            case TYPE_C_STICKY: nC++; break;
            case TYPE_D_RESPONSIVE: nD++; break;
        }
    }
    
    std::cout << "#CELLS: A=" << nA << ";  B=" << nB << "; C=" << nC << "; D=" << nD << std::endl;
}

std::string IOManager::generateFileName(const std::string& prefix, int time_step) const {
    std::stringstream ss;
    ss << "Alpha" << params_.alpha_stiff << "/" << prefix;
    if (time_step >= 0) {
        ss << time_step;
    }
    return ss.str();
}

void IOManager::createOutputDirectory() const {
    std::stringstream ss;
    ss << "Alpha" << params_.alpha_stiff;
    std::string dir_name = ss.str();
    
    // Create directory (Unix/Linux/Mac)
    mkdir(dir_name.c_str(), 0755);
} 