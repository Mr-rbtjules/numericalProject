# Solving PDEs Efficiently in C Using Numerical Methods

A geometric multigrid solver and multigrid-preconditioned conjugate gradient for the 2D Laplace equation on non-trivially shaped domains (rectangles with rectangular holes), implemented in C with CSR sparse storage.

---

## Overview

This project implements efficient iterative solvers for the 2D Laplace equation discretized with the standard 5-point finite difference stencil. The domain is a rectangle with a rectangular hole cut out — the hole can be interior or touching an edge, and all boundary cases are handled. The main solvers are a **geometric multigrid method** (V-, W-, and F-cycles) and a **preconditioned conjugate gradient** using multigrid as the preconditioner.

---

## The Problem

**PDE**: 2D Laplace equation with Dirichlet boundary conditions

```
-Δu = 0    on Ω
 u = g     on ∂Ω
```

where Ω is a rectangle with a rectangular hole, and g(x,y) = sin(√(x² + y²)) on all boundaries (outer edges and hole edges). Eight predefined domain configurations are available.

**Discretization**: Standard 5-point stencil on a uniform Cartesian grid with spacing h:

```
         -1/h²
-1/h²    4/h²   -1/h²
         -1/h²
```

The resulting sparse linear system Au = b is stored in **CSR format** (Compressed Sparse Row).

---

## Solvers

### Geometric Multigrid

The multigrid method supports three cycle types:

| Cycle | Description |
|-------|-------------|
| **V-cycle** | Single recursion down and up |
| **W-cycle** | Two recursive calls per level |
| **F-cycle** | V-cycle with additional sweeps going back down |

**Components:**
- **Pre-smoothing**: Forward Gauss-Seidel
- **Post-smoothing**: Backward Gauss-Seidel
- **Restriction**: Full-weighting operator (averaging with cardinal neighbors)
- **Prolongation**: Bilinear interpolation (4 cases: odd-odd, odd-even, even-odd, even-even)
- **Coarse solver**: UMFPACK direct solver, symmetric Gauss-Seidel, or Jacobi (compile-time selectable)

All multigrid levels are stored in **flat contiguous arrays** with offset indexing for cache efficiency.

### Multigrid-Preconditioned Conjugate Gradient

The PCG method uses 2 multigrid V-cycle iterations as the preconditioner at each CG step, combining the rapid error smoothing of multigrid with the optimal convergence of CG.

### Unpreconditioned Conjugate Gradient

A standalone CG implementation is also available for comparison.

---

## Convergence

Sample residual history (15 V-cycle iterations):

```
Iter    ||b - Au||₂
 1      569.8
 2       24.6
 3        0.66
 5        0.0019
 8        1.1e-7
11        1.8e-11
```

Approximately one order of magnitude reduction per iteration — characteristic of a well-functioning multigrid solver.

---

## Project Structure

```
numericalProject/
├── proto.h       # Macros, structs, function prototypes, compile-time config
├── main.c        # Entry point (calls mg_method or CGmethod)
├── method.c      # Multigrid cycles, smoothers, restriction, prolongation
├── prob.c        # Problem assembly: CSR matrix and RHS for each grid level
├── cg.c          # Preconditioned CG and plain CG
├── globVal.c     # Domain parsing, multi-level index offset computation
├── tools.c       # Linear algebra: dot product, vector ops, residual
├── plot.c        # Gnuplot visualization, cycle ASCII art
├── umfpack.c     # UMFPACK direct solver wrapper
├── time.c        # CPU timer
├── Makefile       # Build with SuiteSparse, METIS, OpenMP, BLAS, LAPACK
└── oldCode/      # Previous iterations
```

---

## Building

Requires SuiteSparse (UMFPACK), METIS, OpenMP, BLAS, and LAPACK. On macOS with Homebrew:

```bash
make
./main
```

Compile-time flags in `proto.h`:

| Flag | Options | Description |
|------|---------|-------------|
| `MODE` | 0 / 1 / 2 | Coarse solver: UMFPACK / sym. Gauss-Seidel / Jacobi |
| `BOUND` | BOUND0/1/2 | Boundary condition function |
| `DOMAIN` | DOMAIN1–8 | Domain geometry |
| `PLOT` | 0 / 1 | Enable Gnuplot visualization |
| `CHRONO` | 0 / 1 | Print timing |

---

## Visualization

- **3D surface plot** of the solution via Gnuplot (NaN masking for hole regions)
- **Convergence plot**: semilog residual norm vs. iteration (`iter.png`)
- **Cycle structure**: ASCII art of V/W/F-cycle printed to terminal
