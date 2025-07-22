#include "Simulation.h"
#include "IO.h"
#include "Constants.h"
#include <iostream>
#include <cstdlib>
#include <ctime>

// Global debug flags (defined in Constants.h)
int g_debug_stop = 0;
int g_debug_speak = 0;

int main(int argc, char* argv[]) {
    std::cout << "## ATTENTION! HAVE YOU CREATED THE FOLDER IN WHICH I SHOULD WRITE??? ##" << std::endl;
    std::cout << "Type: 1. name output 2. cell adherence (stiff cells) 3. Fraction of B=receiver cells " 
              << "4. Switch Rate ON B->C; 5. Switch Rate OFF C->B; 6. cell adherence (stiffer cells) "
              << "7. cell adherence (stiffest cells) 8. Switch rate on A->D 9. Switch rate off D->A" << std::endl;
    
    // Initialize random seed
    srand(time(NULL));
    
    try {
        // Parse command line arguments
        SimulationParameters params = IOManager::parseCommandLine(argc, argv);
        
        // Create and run simulation
        MorphogenesisSimulation simulation(params);
        simulation.run();
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Unknown error occurred" << std::endl;
        return 1;
    }
    
    return 0;
} 