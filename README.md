# 2D Field Solvers

FINITE DIFFERENCE · ELECTROSTATICS · C (2017), PYTHON (IN PROGRESS)

Finite-difference field solvers on a 2D grid. The same five-point stencil runs
all of them; they differ only in how many time derivatives they carry —
Laplace none, heat one, waves two.

**Status: the 2017 C is here, the Python rewrite is not yet.** The repository
was called `2d-heat-equation-solver` until October 2026, which is how it stayed
lost — the code solves electrostatics. That name was not wrong, only early:
steady-state heat and electrostatics are the same Laplace equation, and
relaxation is a diffusion process. `heat.py` will make it literal.

## What is here

```
t1/                 the 2017 Métodos Computacionales submission, unchanged
  placas.c          serial relaxation; "placas" = plates
  sinparalelisar/   the serial version that worked
  Entrega/          the MPI submission, which never compiled
docs/superpowers/   the design and plan for the rewrite
```

## The result

A parallel-plate capacitor, plates 2 m long and 1 m apart in a 5 m grounded
box, at ±50 V. Measured on a 251 × 251 grid:

| Parameter | Value | Unit | Source |
|---|---:|---|---|
| Field, gap centre | 100.2 | V/m | Measured |
| Field, ideal V/d | 100.0 | V/m | Calculated |
| C′, from the grid | 3.06 × 10⁻¹¹ | F/m | Measured |
| C′, textbook ε₀l/d | 1.77 × 10⁻¹¹ | F/m | Calculated |
| Edge excess | 72.6 | % | Calculated |

The textbook formula is wrong by 73 % for this geometry, because it assumes the
field outside the plates away. C′ is integrated over the whole box, whose
grounded walls truncate the fringing field, so it is quoted with its box.

## The 2017 code

The serial version works. The MPI version does not compile:

```
Entrega.c:160:12: error: expected ';' after expression
  160 |         V[i]=V1[i]
```

Seven defects in total, each with its fix, are written up in `legacy/BUGS.md`
once the restructure lands. The architecture of the parallel version is sound —
flat-array decomposition with halo rows and one `MPI_Allgather` per sweep, and
the Dirichlet condition re-stamped each sweep rather than masked per cell,
which is the right call for a decomposed domain.

## Next

`solvers/grid.py`, `laplace.py`, `plot.py` and their tests, then `heat.py` and
`wave.py`. Plan in `docs/superpowers/plans/`.

Written up at
[sergiosgutierrez.com/writing/the-capacitor-the-textbook-hides](https://sergiosgutierrez.com/writing/the-capacitor-the-textbook-hides).
