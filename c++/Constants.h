#ifndef CONSTANTS_H
#define CONSTANTS_H

// Simulation Parameters
const int MAX_MONTE_CARLO_SWEEPS = 10000;
const int GRID_SIZE = 60;
const int TARGET_CELL_AREA = 16;
const int TOTAL_CELLS = int(GRID_SIZE * GRID_SIZE * 1.0 / TARGET_CELL_AREA);

// Physical Parameters
const double LAMBDA_AREA_CONSTRAINT = 1.0;
const double ALPHA_SOFT = 1.0;
const double ALPHA_STIFF_DEFAULT = 1.0;
const double PRESSURE = 0.0;

// Cell Types
enum CellType {
    TYPE_A_SENDER = 0,
    TYPE_B_RECEIVER = 1,
    TYPE_C_STICKY = 2,
    TYPE_D_RESPONSIVE = 3
};

// Boundary Interaction Types
enum BoundaryType {
    BOUNDARY_AA = 0, BOUNDARY_AB = 1, BOUNDARY_AC = 2, BOUNDARY_BA = 3,
    BOUNDARY_BB = 4, BOUNDARY_BC = 5, BOUNDARY_CA = 6, BOUNDARY_CB = 7,
    BOUNDARY_CC = 8, BOUNDARY_AD = 9, BOUNDARY_BD = 10, BOUNDARY_CD = 11,
    BOUNDARY_DD = 12, BOUNDARY_DA = 13, BOUNDARY_DB = 14, BOUNDARY_DC = 15
};

// Neighbor directions (8-connected)
const int NEIGHBOR_DIRECTIONS[8][2] = {
    {+1, 0}, {-1, 0}, {0, +1}, {0, -1},
    {+1, +1}, {+1, -1}, {-1, +1}, {-1, -1}
};

// Debug flags
extern int g_debug_stop;
extern int g_debug_speak;

#endif // CONSTANTS_H 