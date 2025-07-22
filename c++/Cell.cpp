#include "Cell.h"
#include <algorithm>
#include <cstdlib>

Cell::Cell() : id_(-1), type_(TYPE_A_SENDER), alpha_(ALPHA_STIFF_DEFAULT), 
               lambda_(LAMBDA_AREA_CONSTRAINT), area_(0), com_x_(0), com_y_(0),
               unwrapped_com_x_(0), unwrapped_com_y_(0), boundary_counts_(16, 0) {
}

Cell::Cell(int id, CellType type, double alpha) : id_(id), type_(type), alpha_(alpha),
                                                  lambda_(LAMBDA_AREA_CONSTRAINT), area_(0),
                                                  com_x_(0), com_y_(0), unwrapped_com_x_(0),
                                                  unwrapped_com_y_(0), boundary_counts_(16, 0) {
}

bool Cell::attemptGeneSwitch(double switch_probability) {
    double random_value = static_cast<double>(rand()) / RAND_MAX;
    return random_value < switch_probability;
} 