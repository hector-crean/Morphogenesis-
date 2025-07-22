#ifndef CELL_H
#define CELL_H

#include "Constants.h"
#include <vector>

class Cell {
public:
    // Constructors
    Cell();
    Cell(int id, CellType type, double alpha = ALPHA_STIFF_DEFAULT);
    
    // Getters
    int getId() const { return id_; }
    CellType getType() const { return type_; }
    double getAlpha() const { return alpha_; }
    double getLambda() const { return lambda_; }
    int getArea() const { return area_; }
    double getCOMX() const { return com_x_; }
    double getCOMY() const { return com_y_; }
    double getUnwrappedCOMX() const { return unwrapped_com_x_; }
    double getUnwrappedCOMY() const { return unwrapped_com_y_; }
    
    // Setters
    void setType(CellType type) { type_ = type; }
    void setAlpha(double alpha) { alpha_ = alpha; }
    void setLambda(double lambda) { lambda_ = lambda; }
    void setArea(int area) { area_ = area; }
    void setCOM(double x, double y) { com_x_ = x; com_y_ = y; }
    void setUnwrappedCOM(double x, double y) { unwrapped_com_x_ = x; unwrapped_com_y_ = y; }
    
    // Area operations
    void incrementArea() { area_++; }
    void decrementArea() { area_--; }
    
    // Boundary tracking
    void incrementBoundary(BoundaryType type) { boundary_counts_[type]++; }
    void resetBoundaries() { std::fill(boundary_counts_.begin(), boundary_counts_.end(), 0); }
    int getBoundaryCount(BoundaryType type) const { return boundary_counts_[type]; }
    
    // Gene switching probabilities
    bool attemptGeneSwitch(double switch_probability);
    
private:
    int id_;
    CellType type_;
    double alpha_;           // Stiffness parameter
    double lambda_;          // Area constraint parameter
    int area_;              // Current area
    double com_x_, com_y_;  // Center of mass
    double unwrapped_com_x_, unwrapped_com_y_;  // Unwrapped COM for PBC
    
    // Boundary interaction counts
    std::vector<int> boundary_counts_;
};

#endif // CELL_H 