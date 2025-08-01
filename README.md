

# Morphogenesis

![Modelling Self-Assembly of Tissues in Silico](modelling_self-assembly_of_tissues_in_silico_poster.jpg)

**📄 [Read the Full Thesis: Cell_Sorting.pdf](Cell_Sorting.pdf)**

## Abstract

Genetically encoded algorithms control individual and collective cell behaviour. We aim to reconstruct these programs in silico to understand their logic. In vitro and in vivo studies have identified many genes that control cell-cell signalling and cell morphology. Most development systems, however, use contact cell-cell signalling interactions to induce morphological responses. We consequently explored whether simple synthetic circuits, in which morphological changes are driven by contact cell-cell signalling interactions, could suffice to generate self-organizing multicellular structures.

## Quick Use Guide

### 1. Compilation

`gwitch21.cp` is the source file for the standard Cellular Potts Model.

**Compile:**
```bash
c++ gwitch21.cp -o gs
```

### 2. Running the Simulation

**Usage:**
```bash
./gs outname(string) cell_adherence_alpha(0-infinite) Fraction_B_receiver_cells(0<f<1) switch_rate_ON_B->C(0-1) switch_Rate_OFF_C->B(0-1) cell_adherence_stiffer_cells(0-infinite) cell_adherence_stiffest_cells(0-infinite) switch_rate_ON_A->D(0-1) switch_rate_off_D->A(0-1)
```

**Example:**
```bash
./gs out 1 0.5 0.3 0.1 1.5 1.8 0.1 0.05
```

**Note:** In the simulation each cell type has a different adhesiveness, which affects the form of the adhesion matrix. In the accompanying code, these are referred to as:
- **α_soft** (B cell)
- **α_stiff** (A cell) 
- **α_stiffer** (D cell)
- **α_stiffest** (C cell)

Where the adhesiveness/stiffness hierarchy is: `α_soft < α_stiff < α_stiffer < α_stiffest`

Two additional control parameters (`k_on2` and `k_off2`) control the switch rate between type A and type D cells.

### 3. Output Files

**Lattice coordinates** (x,y), spin id and Alpha are reported in:
```
Alpha_1/out_t*
```

### 4. Center of Mass Data

**Center of mass positions** are reported in:
```
Alpha1/com_out
```

**Format:**
```
#t1
CellID x y 
..
CellID x y 

#t2
CellID x y 
..
CellID x y 
```

**Plotting examples:**
- First frame: `p "com_out" index 0 u 2:3 w p pt 7`
- Second frame: `p "com_out" index 1 u 2:3 w p pt 7`
- And so on...

### 5. Creating Animations

Use `load "Plot_movie.plt"` in the `MakeAnimation` folder to make animations of cell movements.

### 6. Mean Squared Displacement (MSD) Analysis

In the `Analysis` folder, use `ComputeMSD`:

**Usage:**
```bash
./ComputeMSD ../Alpha1.2/com_out 3000 10 225 60 out
```

This computes mean squared displacement of center of mass from `com_out`:
- For 275 frames
- With lagtime 10 Monte Carlo sweeps
- 225 cells on a 60x60 grid
- Writing to `MSD_out.dat`

**Example:**
```bash
./ComputeMSD ../Alpha1.4/com_out 270 10 225 60 out
```

### 7. Plotting MSD Results

**Average MSD across cells and time:**
```gnuplot
p "MSD_out.dat" i 0 u 1:2
```

**Individual cell MSD (e.g., cell 20):**
```gnuplot
p "MSD_out.dat" i 20 u 1:2
```

**All cells with average highlighted:**
```gnuplot
p for [i=1:225] "MSD_out.dat" index i u 1:2 w l lw 0.5 not, "MSD_out.dat" i 0 w l lw 5 lt -1 title "average MSD"
```

### 8. Center of Mass Output Format

`ComputeMSD` also creates `COM_out.dat` with the format:
```
#Cell id
time cell_id x_pbc y_pbc x_unwrapped y_unwrapped
...
```

The unwrapped coordinates are used to calculate MSD.

**Note:** `Ncells = 60×60/16 = 225` where 16 is the typical cell area.

## Background on Morphogenesis

Cells can differentiate into specific cell types that spatially self-organise into tissue architectures. This is achieved by cells communicating, and inducing the expression of proteins, in one-another. In embryology, there is specific interest in how regulatory programs encoded in the DNA of a single fertilized cell can integrate together into the gross organisation of a whole organism. It is still an open question as to how compact DNA programs encode algorithms that allow individual cells to construct complex macrostructures. 

The accompanying code involves understanding how genetically encoded algorithms control these individual and collective cell behaviours. Here we reconstruct in silico a 'toy model' of cell-cell interactions, attempting to find metrics that correlate with:
- Cell self-organisation into multi-layered structures
- Divergence of genotypically identical cells into distinct types
- Symmetry breaking

Computational modelling allows us to track the trajectories of these multi-cellular systems through **morphospace** (space of forms). Plotting the energy landscape corresponding to each point in the morphospace allows an understanding of stable and meta-stable states in the morphospace, and how tissues transition between them.

While this code is primarily focused on the biological realm, it is equally applicable in studying animal collective behaviours, such as flock formation and sociality of ants, which also arise in the absence of central control, and are generated by local interaction among individuals. Within the architectural sphere, morphogenetic programming finds expression in the construction of parametric and generative structures.

## Biological Implementation Details

Cell-to-cell communication networks comprise both intra- and intercellular processes. In order to simplify modelling, intracellular signal transduction networks are generally treated as "black boxes" with specified input-to-output response relationships. This abstracts molecular detail yet still captures essential dynamic properties of the system.

### Cell Signalling Mechanics

Cells have a cell-specific density of ligands/receptors on their surface. For two adjacent cells, the signalling strength between them is a function of the number of ligands/receptors in contact. Varying degrees of contact are correlated with differing signalling strengths. 

In the biological realm, furthermore, the response of a cell to this signalling is an integral of the signal strength over time, and so is not just dependent on the instantaneous ligand-receptor contact, but the history of this contact in the recent past. In addition, rather than certain ligand-receptor interactions leading precisely to morphological outputs, the same ligand-receptor bonding can lead to different outputs, only differentiated by the magnitude of the signal strength as well as the time over which the signalling occurs. 

For the sake of simplicity, our in silico model assumes a one-to-one correspondence between signalling and a morphological output.

### Implementation Strategy

Cell-cell contact dependent signalling can be implemented via a points-based system, where cells gain 'points' for having one of its lattice points being in contact with the lattice point of another cell. For this purpose, we can construct an array tracking the **heterotrophic boundary length** for each cell. This comprises a 6 by 225 array, with details of the number of contact points each cell has with other cell types at a given moment. This information can then be used to calculate the point at which a cell has had sufficient stimulation to express cadherin.

An additional concern arises when considering the two methods by which morphological changes are activated: **gene inhibition** and **gene activation**. Gene activation, for instance, has 6 main forms (in terms of the biomechanics), relating to how the RNA/DNA polymerase acts. This means that gene activation is often not permanent. If cell A can cause a change in cell B to state B', B' can transition back to B after a 'decay time', if it does not receive sufficient signalling through binding with cell A's receptors. The duration of gene activation can be altered by playing with stability of transcription factors. We must therefore also incorporate in our simulation a mechanism whereby cells can transition back to their original state.

### Simplified Gene Switch Implementation

While we did incorporate the functionality to control cell switching via heterotopic boundary lengths, it was decided that the relevant features could be elegantly modelled more simply. During the metropolis step of the CPM, we introduced a 'gene switch' step. If the spins of two adjoining lattice points are different (indicating different cell types), a random number between 0-1 is generated. If this is less than the `k_on` number, then the cell that the lattice point belongs to will change type. In this way, it simulates the fact that more cell contact with sender cells leads to larger ligand exposure, and so greater probability of a receiver cell changing type. Although it only incorporates information from one timestep, the speed of the simulation is fast.
