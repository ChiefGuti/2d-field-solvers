# 2D Field Solvers — Phase 1 (repo) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn a misnamed 2017 coursework repo into a documented NumPy electrostatics solver whose headline result is that the textbook parallel-plate formula is 73% wrong for its own geometry.

**Architecture:** One shared grid module holds the geometry, the Dirichlet mask and the 5-point Laplacian; `laplace.py` relaxes it to equilibrium with red-black SOR; `plot.py` renders any solved grid without knowing which physics produced it. The 2017 C moves to `legacy/` unchanged, with its seven defects and their fixes recorded in `legacy/BUGS.md`. Later phases add `heat.py` and `wave.py` against the same `grid.py` and `plot.py`.

**Tech Stack:** Python 3 with NumPy and Matplotlib only. No package, no build step, no new dependencies. `/usr/bin/python3` is the interpreter that has NumPy and Matplotlib on this machine — the Homebrew 3.13/3.14 do not.

**Spec:** `docs/superpowers/specs/2026-10-10-2d-field-solvers-design.md`

**Scope note:** The spec's phase 1 covers two independent subsystems. This plan is the **repo** only. The website (`src/lib/capacitor.ts`, `capacitor-field.tsx`, the WR-04 post, the `/research` connection) is a separate plan written after this one lands, because it is a different repo with its own test and review gate and can ship independently.

## Global Constraints

- **Never commit or push without Sergio's explicit approval in chat.** This overrides the commit steps below: every "Commit" step means *stage the change, show the diff, and wait*. No exceptions, any branch, any repo.
- **Never delete a file not named in this plan.** The two deletions here are pre-approved (`t1/Entrega/a.out`, `t1/sinparalelisar/output.txt`); anything else needs a fresh yes.
- SI units throughout: metres, volts, farads per metre. Every documented value carries a unit and a source.
- `EPS0 = 8.8541878128e-12` F/m, defined once in `grid.py`.
- Grid size `m` is **odd**, so the box has a true centre row and the geometry is exactly symmetric. This is a deliberate departure from the 2017 code's `-1` index offset (see Review Focus).
- Datasheet voice in all prose: outside reader, present tense, definite. Sections numbered `## 01`, figures captioned `FIG. n · what it is`, spec tables Parameter · Value · Unit · Source with numbers right-aligned. Reference: `../../hardware/docs/DESIGN.md` §3.1 and the WF-01 sheet.
- No new code-number series. The solvers are named in prose: `laplace`, `heat`, `wave`.
- `legacy/` sources are the 2017 submission, byte-for-byte. Fixes are recorded in `BUGS.md`, never applied to those files.

## Review Focus

Input classes the spec implies but no task's own tests would otherwise exercise. Each has its test added to the task that owns the code.

1. **Even `m`** — the 2017 bug class, reintroduced. With even `m` there is no centre row, the plates land half a cell off-centre and antisymmetry silently degrades from 0 to ~1e-1. Verified 2026-10-10 on the original output. `grid.py` must reject even `m` rather than return a subtly wrong grid. (Task 3)
2. **Plates wider than the box** — `plate_len >= L` indexes outside the array or silently clips, giving a capacitance for a geometry the caller did not ask for. Must raise. (Task 3)
3. **Non-convergence** — `max_sweeps` reached without hitting `tol` currently returns a half-relaxed grid that looks plausible and is wrong. Must be visible in the returned metrics, not silent. (Task 4)
4. **`gap` smaller than one cell** — `round(gap/(2h))` collapses to 0, both plates land on the same row and the field is undefined. Must raise. (Task 3)
5. **Plates close to the walls** — at `l/d = 4` in an `L = 5` box the grounded walls distort the result (measured: the C′ ratio flattens from 1.568 to 1.542, partly wall effect, not physics). The l/d sweep must not silently run into this; the figure caps where the box stops being big enough. (Task 6)

---

### Task 1: Repo hygiene and the legacy move

**Files:**
- Create: `.gitignore`
- Delete: `t1/Entrega/a.out`, `t1/sinparalelisar/output.txt`
- Move: `t1/` → `legacy/t1/`
- Restore: `legacy/t1/sinparalelisar/placas.c`, `plotter.py`, `V.png` to their 2017 state

**Interfaces:**
- Consumes: nothing.
- Produces: `legacy/t1/` holding the unmodified 2017 submission. No Python yet.

- [ ] **Step 1: Confirm what the working tree has changed since 2017**

```bash
cd /Users/sergiogutierrez/Documents/Projects/research/repos/2d-field-solvers
git status --short
git diff --stat
```

Expected: three modified files under `t1/sinparalelisar/` (`placas.c`, `plotter.py`, `V.png`) and untracked `docs/`, `output.txt`, `placas`.

- [ ] **Step 2: Save the 2026 patch before reverting it**

The patch is the raw material for `BUGS.md` in Task 2. Keep it outside the repo so it is not committed.

```bash
git diff t1/sinparalelisar/placas.c > /tmp/2026-fixes.patch
cat /tmp/2026-fixes.patch
```

Expected: the calloc, `V_new` pointer and tolerance hunks.

- [ ] **Step 3: Revert the 2017 sources to as-submitted**

```bash
git checkout -- t1/sinparalelisar/placas.c t1/sinparalelisar/plotter.py t1/sinparalelisar/V.png
git status --short
```

Expected: no `M` lines remain. Only untracked `docs/`, `output.txt`, `placas`.

- [ ] **Step 4: Delete the two approved files**

Both confirmed by Sergio on 2026-10-10. `a.out` is a 16 KB Mac binary committed in 2017 that no longer runs; `output.txt` is 608 KB of regenerable grid data. **Delete these two only.**

```bash
git rm t1/Entrega/a.out
rm -f t1/sinparalelisar/output.txt t1/sinparalelisar/placas
```

- [ ] **Step 5: Write `.gitignore`**

```gitignore
# Build products
a.out
*.o
placas
solvers/__pycache__/
__pycache__/

# Generated data and figures
output.txt
*.npy
out/
```

- [ ] **Step 6: Move the 2017 tree under legacy/**

```bash
git mv t1 legacy
ls legacy
```

Expected: `legacy/1placas.c  legacy/Entrega  legacy/Entrega.c  legacy/Untitled.ipynb  legacy/placas.c  legacy/sinparalelisar`

- [ ] **Step 7: Verify the legacy sources are byte-identical to 2017**

```bash
git diff HEAD --stat -- legacy
```

Expected: only renames (`R`), no content changes. If any line of `.c` or `.py` differs, Step 3 was incomplete — fix before continuing.

- [ ] **Step 8: Stage and show the diff. Do not commit.**

```bash
git add -A
git status --short
```

Then show Sergio the summary and wait for his go-ahead before `git commit -m "chore: move 2017 submission to legacy/, drop build products"`.

---

### Task 2: `legacy/BUGS.md`

**Files:**
- Create: `legacy/BUGS.md`

**Interfaces:**
- Consumes: `/tmp/2026-fixes.patch` from Task 1.
- Produces: the document §06 of the README links to.

- [ ] **Step 1: Write the file**

Four parts per bug — Symptom, Cause, Fix, Why it matters — per spec §5. The fixes are recorded here and deliberately **not** applied to the `legacy/` sources.

````markdown
# Bugs in the 2017 submission

Found on 2026-10-10 and 2026-10-11, nine years after the code was handed in.
The sources in this directory are unchanged: these are the fixes, recorded
beside the mistake rather than applied over it.

Finding them is most of what the original exercise was worth.

**There are two MPI files and they differ.** `../Entrega.c` has `h = 0.02`;
`Entrega/Entrega.c`, the submission, has `h = 5/512`. Line numbers below are
the submission's.

## The headline: it never compiled

```
Entrega.c:160:12: error: expected ';' after expression
  160 |         V[i]=V1[i]
```

The parallel half of the assignment — the actual deliverable — could not be
built by anyone, including its author. Everything below sits behind that.

## 1. Missing semicolon, Entrega/Entrega.c:160

**Symptom.** `mpicc` stops. No binary.

**Cause.** `V[i]=V1[i]` in the copy-back loop.

**Fix.** `V[i]=V1[i];`

**Why it matters.** Nothing else in this file was ever exercised.

## 2. Assignment where an increment was meant, Entrega/Entrega.c:162

**Symptom.** If it compiled, it would not stop.

**Cause.** `n=+1;`

**Fix.** `n+=1;`

**Why it matters.** `=+` is two tokens that read like one operator. `n = +1`
sets the counter to positive one every sweep, so `while(n<N)` with `N = 2*m*m`
never ends.

## 3. Integer division, Entrega/Entrega.c:10

**Symptom.** Every plate coordinate degenerates.

**Cause.** `double h = 5/512;`

**Fix.** `double h = 5.0/512.0;`

**Why it matters.** Both operands are `int`, so C does the division in integer
arithmetic and gets 0, *then* converts to `double`. The declared type on the
left does not reach into the expression on the right. With `h = 0`, `l/(h*2)`
divides by zero.

## 4. Dead branch for the last rank, Entrega/Entrega.c:132

**Symptom.** The final slab of the grid is never relaxed.

**Cause.** `if(world_rank==world_size)`

**Fix.** `if(world_rank==world_size-1)`

**Why it matters.** MPI ranks run 0 to `world_size-1`. The test is never true.

## 5. Always-true rank test, Entrega/Entrega.c:85

**Symptom.** Rank 0 relaxes its slab twice, the second time over the wrong range.

**Cause.** `if (!(world_rank==0 && world_rank==world_size))`

**Fix.** `if (world_rank!=0 && world_rank!=world_size-1)`

**Why it matters.** A rank cannot be both 0 and `world_size`, so the negation is
always true and every rank falls into the middle-rank branch as well as its own.
Like 4, it only misbehaves with more than two processors.

## 6. Uninitialised interior, sinparalelisar/placas.c

**Symptom.** Non-zero values at the box walls; interior values that change
between runs of the same binary.

**Cause.** `double *V = malloc(m*m*sizeof(double));`

**Fix.** `double *V = calloc(m*m, sizeof(double));`

**Why it matters.** The sweep loops run `i = 1` to `m-2`, so row 0, row `m-1`,
column 0 and column `m-1` are never written. Those cells are not scratch space —
they are the grounded box the capacitor sits in, and the solution depends on
them being zero. Relaxation faithfully propagates whatever the allocator held.

## 7. Aliased arrays, sinparalelisar/placas.c

**Symptom.** None visible. The answer is right.

**Cause.** `V_new = init(x0, x1, y0, y1, V);`

**Fix.** `V_new = init(x0, x1, y0, y1, V_new);`

**Why it matters.** `init` returns the pointer it is handed, so `V_new` was
reassigned to `V`. One array under two names, the freshly malloc'd buffer
leaked, and the copy-back loop assigned every cell to itself. That turns the
Jacobi iteration the two-array structure expresses into Gauss-Seidel, which
converges about twice as fast. The bug helped, which is why it was never found.

## Not a bug: re-stamping the plates

The MPI sweep has no plate guard, which looks like a defect next to the serial
version's `if (!(j >= x0 && ...))`. It is not. `init()` is called at the top of
every iteration:

```c
while(n<N)
{
	/*cada procesador restablece las condiciones de las varillas despues de cada iteracion*/
	V = init(x0,x1,y0,y1,V);
```

Let the averaging damage the plates, then write them back before anyone reads
them. Once the domain is cut into slabs, that is the better of the two designs:
no per-cell ownership test, same answer, and the comment says it was deliberate.

## Not a bug: the repository name

The repo was called `2d-heat-equation-solver` until 2026-10-10. The code solves
electrostatics. Both are true: steady-state heat and electrostatics are the same
Laplace equation, and Jacobi relaxation *is* a diffusion process. The name
described the method rather than the problem. It is now `2d-field-solvers`, and
`heat.py` will make the original name literal.
````

- [ ] **Step 2: Check every quoted line against the real file**

Line numbers and quoted code drift. Verify each:

```bash
cd /Users/sergiogutierrez/Documents/Projects/research/repos/2d-field-solvers
sed -n '23p;24p;27p' legacy/sinparalelisar/placas.c
sed -n '10p;85p;132p;160p;162p' legacy/Entrega/Entrega.c
```

Expected, in order: the `malloc` pair, `V_new = init(..., V)`, then `double h = 5/512;`, the always-true rank test, the dead last-rank branch, `V[i]=V1[i]`, and `n=+1;`. Correct the line numbers in `BUGS.md` to whatever these commands actually print — the move into `legacy/` does not shift them, but an accidental reformat would.

Then confirm the headline claim still holds, rather than trusting this plan for it:

```bash
mkdir -p /tmp/mpistub && cat > /tmp/mpistub/mpi.h <<'H'
#define MPI_COMM_WORLD 0
#define MPI_DOUBLE 0
int MPI_Init(void*,void*); int MPI_Comm_rank(int,int*); int MPI_Comm_size(int,int*);
int MPI_Allgather(void*,int,int,void*,int,int,int); int MPI_Finalize(void); int MPI_Barrier(int);
H
gcc -fsyntax-only -I/tmp/mpistub legacy/Entrega/Entrega.c
```

Expected: `error: expected ';' after expression` at the `V[i]=V1[i]` line. If it compiles clean, the headline of `BUGS.md` is wrong and must be rewritten before going further.

- [ ] **Step 3: Stage and show. Do not commit.**

```bash
git add legacy/BUGS.md && git diff --cached --stat
```

Wait for approval before `git commit -m "docs: record the seven 2017 defects and their fixes"`.

---

### Task 3: `solvers/grid.py`

**Files:**
- Create: `solvers/grid.py`
- Test: `solvers/test_solvers.py`

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `EPS0: float`
  - `geometry(m=251, L=5.0, plate_len=2.0, gap=1.0, v0=100.0) -> Grid`
  - `Grid` is a dataclass with fields `V: np.ndarray`, `fixed: np.ndarray` (bool), `h: float`, `m: int`, `L: float`, `plate_len: float`, `gap: float`, `v0: float`, `centre: int`, `y0: int`, `y1: int`, `x0: int`, `x1: int`.
  - `plate_len` and `gap` on the returned `Grid` are the **snapped** values — what the grid actually represents after rounding to whole cells — not what was requested.

- [ ] **Step 1: Write the failing tests**

Create `solvers/test_solvers.py`:

```python
"""Checks for the 2D field solvers. Run: /usr/bin/python3 solvers/test_solvers.py"""
import numpy as np
from grid import geometry, EPS0


def test_odd_m_has_a_true_centre():
    g = geometry(m=251)
    assert g.centre == 125
    assert g.y0 == g.centre - (g.centre - g.y0)
    # plates sit symmetrically about the centre row
    assert g.centre - g.y0 == g.y1 - g.centre


def test_even_m_is_rejected():
    # The 2017 bug class: with even m there is no centre row, the plates land
    # half a cell off-centre and antisymmetry silently degrades to ~1e-1.
    try:
        geometry(m=250)
    except ValueError as e:
        assert "odd" in str(e).lower()
    else:
        raise AssertionError("even m must raise ValueError")


def test_plates_wider_than_the_box_are_rejected():
    try:
        geometry(m=101, L=5.0, plate_len=5.0)
    except ValueError as e:
        assert "plate_len" in str(e)
    else:
        raise AssertionError("plate_len >= L must raise ValueError")


def test_gap_below_one_cell_is_rejected():
    # h = L/m = 5/101 ≈ 0.0495, so gap=0.01 rounds to 0 cells and both plates
    # would land on the same row.
    try:
        geometry(m=101, L=5.0, gap=0.01)
    except ValueError as e:
        assert "gap" in str(e)
    else:
        raise AssertionError("sub-cell gap must raise ValueError")


def test_plates_carry_the_voltage_and_are_fixed():
    g = geometry(m=251, v0=100.0)
    assert np.isclose(g.V[g.y0, g.x0:g.x1 + 1], 50.0).all()
    assert np.isclose(g.V[g.y1, g.x0:g.x1 + 1], -50.0).all()
    assert g.fixed[g.y0, g.x0:g.x1 + 1].all()
    assert g.fixed[g.y1, g.x0:g.x1 + 1].all()


def test_walls_are_fixed_at_zero():
    g = geometry(m=251)
    for edge in (g.fixed[0, :], g.fixed[-1, :], g.fixed[:, 0], g.fixed[:, -1]):
        assert edge.all()
    assert g.V[0, :].max() == 0.0 and g.V[-1, :].max() == 0.0


def test_snapped_dimensions_are_reported_not_requested():
    # h = 5/251 ≈ 0.01992, so a 2.0 m plate cannot land exactly on cells.
    g = geometry(m=251, plate_len=2.0, gap=1.0)
    assert g.plate_len != 2.0
    assert abs(g.plate_len - 2.0) < 2 * g.h
    assert np.isclose(g.gap, (g.y1 - g.y0) * g.h)


if __name__ == "__main__":
    import sys
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"ok  {fn.__name__}")
    print(f"\n{len(fns)} passed")
    sys.exit(0)
```

- [ ] **Step 2: Run them to verify they fail**

```bash
cd /Users/sergiogutierrez/Documents/Projects/research/repos/2d-field-solvers/solvers
/usr/bin/python3 test_solvers.py
```

Expected: `ModuleNotFoundError: No module named 'grid'`

- [ ] **Step 3: Write `solvers/grid.py`**

```python
"""
Shared geometry for the 2D field solvers: a square box on a square grid, with
two conducting plates held at ±v0/2 and grounded walls.

All three solvers (laplace, heat, wave) use this same grid and the same
Dirichlet mask; they differ only in how many time derivatives they carry.

SI units throughout: metres, volts.
"""
from dataclasses import dataclass

import numpy as np

EPS0 = 8.8541878128e-12  # vacuum permittivity, F/m


@dataclass
class Grid:
    """A box with its boundary conditions already written into V."""
    V: np.ndarray        # potential, volts; plates and walls already set
    fixed: np.ndarray    # True where V may never be updated
    h: float             # cell size, m
    m: int               # cells per side
    L: float             # box side, m
    plate_len: float     # snapped plate length, m
    gap: float           # snapped plate separation, m
    v0: float            # potential difference, V
    centre: int          # the centre row and column
    y0: int              # upper plate row
    y1: int              # lower plate row
    x0: int              # first plate column
    x1: int              # last plate column (inclusive)


def geometry(m=251, L=5.0, plate_len=2.0, gap=1.0, v0=100.0) -> Grid:
    """
    Build the parallel-plate geometry.

    `m` must be odd so the box has a true centre row. With even `m` the plates
    land half a cell off-centre and the solution loses its exact antisymmetry —
    this is the 2017 code's `-1` index offset, and it is rejected rather than
    reproduced.

    `plate_len` and `gap` are snapped to whole cells; the returned Grid reports
    the snapped values, which are what the solution actually describes.
    """
    if m % 2 == 0:
        raise ValueError(f"m must be odd so the box has a true centre row; got {m}")
    if plate_len >= L:
        raise ValueError(f"plate_len {plate_len} m does not fit in a box of {L} m")
    if gap >= L:
        raise ValueError(f"gap {gap} m does not fit in a box of {L} m")

    h = L / m
    centre = m // 2
    half_len = round(plate_len / (2 * h))
    half_gap = round(gap / (2 * h))
    if half_gap < 1:
        raise ValueError(
            f"gap {gap} m is below one cell ({h:.4f} m); both plates would "
            f"land on the same row. Raise gap or raise m."
        )
    if half_len < 1:
        raise ValueError(f"plate_len {plate_len} m is below one cell ({h:.4f} m)")

    x0, x1 = centre - half_len, centre + half_len
    y0, y1 = centre - half_gap, centre + half_gap
    if x0 < 1 or x1 > m - 2 or y0 < 1 or y1 > m - 2:
        raise ValueError("geometry does not fit inside the grounded walls")

    V = np.zeros((m, m))
    fixed = np.zeros((m, m), dtype=bool)
    V[y0, x0:x1 + 1] = v0 / 2
    V[y1, x0:x1 + 1] = -v0 / 2
    fixed[y0, x0:x1 + 1] = True
    fixed[y1, x0:x1 + 1] = True
    fixed[0, :] = fixed[-1, :] = fixed[:, 0] = fixed[:, -1] = True  # grounded box

    return Grid(
        V=V, fixed=fixed, h=h, m=m, L=L,
        plate_len=2 * half_len * h, gap=2 * half_gap * h, v0=v0,
        centre=centre, y0=y0, y1=y1, x0=x0, x1=x1,
    )
```

- [ ] **Step 4: Run the tests to verify they pass**

```bash
cd /Users/sergiogutierrez/Documents/Projects/research/repos/2d-field-solvers/solvers
/usr/bin/python3 test_solvers.py
```

Expected: `7 passed`

- [ ] **Step 5: Stage and show. Do not commit.**

```bash
git add solvers/grid.py solvers/test_solvers.py && git diff --cached --stat
```

Wait for approval before `git commit -m "feat: shared grid geometry with validated boundary conditions"`.

---

### Task 4: `solvers/laplace.py`

**Files:**
- Create: `solvers/laplace.py`
- Modify: `solvers/test_solvers.py` (append)

**Interfaces:**
- Consumes: `grid.geometry`, `grid.Grid`, `grid.EPS0`.
- Produces:
  - `solve(g: Grid, tol=1e-6, max_sweeps=100000, omega=None) -> tuple[np.ndarray, dict]`
  - The metrics dict has keys: `sweeps` (int), `maxdiff` (float), `converged` (bool), `omega` (float), `c_grid` (float, F/m), `c_ideal` (float, F/m), `edge_excess` (float, %), `e_centre` (float, V/m), `e_ideal` (float, V/m), `e_tip` (float, V/m).
  - `field(V, h) -> tuple[np.ndarray, np.ndarray, np.ndarray]` returning `(Ex, Ey, Emag)` in V/m.
  - `capacitance(V, h, v0) -> float` in F/m, by energy integration.
  - Passing `omega=1.0` gives plain Jacobi, which is how the tests check SOR did not change the physics.

- [ ] **Step 1: Write the failing tests**

Append to `solvers/test_solvers.py`, above the `__main__` block:

```python
from laplace import solve, field, capacitance


def test_field_in_the_gap_matches_v_over_d():
    g = geometry(m=251)
    V, mx = solve(g)
    assert mx["converged"]
    assert abs(mx["e_centre"] - mx["e_ideal"]) / mx["e_ideal"] < 0.005


def test_solution_is_exactly_antisymmetric():
    # Verified 2026-10-10: with odd m this is 0.0 to machine precision, because
    # the geometry and the iteration are both symmetric about the centre row.
    g = geometry(m=151)
    V, _ = solve(g)
    c = g.centre
    k = np.arange(1, 40)
    assert np.abs(V[c + k, :] + V[c - k, :]).max() < 1e-12
    assert np.abs(V[c, :]).max() < 1e-12


def test_grounded_walls_stay_at_zero():
    g = geometry(m=151)
    V, _ = solve(g)
    for edge in (V[0, :], V[-1, :], V[:, 0], V[:, -1]):
        assert np.abs(edge).max() == 0.0


def test_sor_and_jacobi_agree():
    # The check that the fast path did not change the answer.
    g1, g2 = geometry(m=81), geometry(m=81)
    v_sor, m_sor = solve(g1, tol=1e-8)
    v_jac, m_jac = solve(g2, tol=1e-8, omega=1.0)
    assert np.abs(v_sor - v_jac).max() < 1e-4
    assert m_sor["sweeps"] < m_jac["sweeps"] / 5   # measured: ~11x at m=96


def test_non_convergence_is_reported_not_hidden():
    g = geometry(m=151)
    V, mx = solve(g, max_sweeps=5)
    assert mx["converged"] is False
    assert mx["sweeps"] == 5
    assert mx["maxdiff"] > 1e-6


def test_edge_effects_raise_capacitance_above_the_ideal():
    g = geometry(m=251)
    V, mx = solve(g)
    assert mx["c_grid"] > mx["c_ideal"]
    assert 60.0 < mx["edge_excess"] < 85.0   # measured 2026-10-10: 72.6 %


def test_edge_excess_falls_as_the_plates_get_longer():
    # The textbook limit, earned. Measured ratios at m=151: 2.881, 2.155,
    # 1.711, 1.568 for l/d = 0.5, 1, 2, 3.
    ratios = []
    for plate_len in (0.5, 1.0, 2.0, 3.0):
        g = geometry(m=151, plate_len=plate_len, gap=1.0)
        _, mx = solve(g)
        ratios.append(mx["c_grid"] / mx["c_ideal"])
    assert all(b < a for a, b in zip(ratios, ratios[1:])), ratios
    assert ratios[0] > 2.5 and ratios[-1] < 1.7


def test_capacitance_is_grid_independent_under_refinement():
    # Measured 2026-10-10: m=201 → 3.0459e-11, m=401 → 3.0732e-11, 0.9 % apart.
    # Convergence is slow because the charge density diverges at the plate tips,
    # so the threshold is 2 %, not the 1 % the spec first assumed.
    _, coarse = solve(geometry(m=201))
    _, fine = solve(geometry(m=401))
    rel = abs(fine["c_grid"] - coarse["c_grid"]) / coarse["c_grid"]
    assert rel < 0.02, f"C' moved {100 * rel:.1f} % from m=201 to m=401"
```

- [ ] **Step 2: Run them to verify they fail**

```bash
/usr/bin/python3 solvers/test_solvers.py
```

Expected: `ModuleNotFoundError: No module named 'laplace'`

- [ ] **Step 3: Write `solvers/laplace.py`**

```python
"""
Laplace's equation on the shared grid: the electrostatic field of a
parallel-plate capacitor.

Every interior cell is replaced by the average of its four neighbours until the
grid stops changing. The plates and the walls are never updated, so they hold
their values and the solution grows out from them.

Solved by red-black successive over-relaxation. Measured 2026-10-10: at m=251
SOR converges in 689 sweeps where the 2017 Jacobi code took 32,813.
"""
import numpy as np

from grid import EPS0, Grid, geometry


def _optimal_omega(m: int) -> float:
    """SOR relaxation factor for an m x m Dirichlet grid. Grid-size dependent:
    computed, never hard-coded."""
    return 2.0 / (1.0 + np.sin(np.pi / m))


def solve(g: Grid, tol=1e-6, max_sweeps=100000, omega=None):
    """
    Relax to equilibrium. Returns (V, metrics).

    omega=None picks the optimal SOR factor; omega=1.0 is plain Jacobi, which
    the tests use to check SOR did not change the physics.

    `metrics["converged"]` is False when max_sweeps ran out. A half-relaxed grid
    looks plausible and is wrong, so callers must check it.
    """
    V = g.V.copy()
    free = ~g.fixed
    i, j = np.indices((g.m, g.m))
    colour = (i + j) % 2
    if omega is None:
        omega = _optimal_omega(g.m)

    avg = np.zeros_like(V)
    maxdiff = np.inf
    sweeps = 0
    for sweeps in range(1, max_sweeps + 1):
        maxdiff = 0.0
        for parity in (0, 1):
            # Recomputed per colour: that is what makes this Gauss-Seidel-like
            # rather than Jacobi. Computing it once for both colours would
            # silently drop back to Jacobi convergence.
            avg[1:-1, 1:-1] = 0.25 * (
                V[:-2, 1:-1] + V[2:, 1:-1] + V[1:-1, :-2] + V[1:-1, 2:]
            )
            sel = free & (colour == parity)
            delta = omega * (avg - V)
            V[sel] += delta[sel]
            maxdiff = max(maxdiff, np.abs(delta[sel]).max(initial=0.0))
        if maxdiff < tol:
            break

    return V, _metrics(V, g, sweeps, maxdiff, maxdiff < tol, omega)


def field(V, h):
    """Electric field from the potential: E = -grad V. Returns (Ex, Ey, |E|)."""
    Ey, Ex = np.gradient(-V, h)   # np.gradient is (row, col) = (y, x)
    return Ex, Ey, np.hypot(Ex, Ey)


def capacitance(V, h, v0):
    """
    Capacitance per unit length, F/m, by energy integration:
        W = (eps0/2) * integral |E|^2 dA,   C' = 2W / V^2

    The integral runs over the whole box, whose grounded walls truncate the
    fringing field, so C' depends on the box size. Always report L with it.
    """
    _, _, Emag = field(V, h)
    W = 0.5 * EPS0 * np.sum(Emag ** 2) * h * h
    return 2 * W / v0 ** 2


def _metrics(V, g, sweeps, maxdiff, converged, omega):
    _, _, Emag = field(V, g.h)
    c_grid = capacitance(V, g.h, g.v0)
    c_ideal = EPS0 * g.plate_len / g.gap
    return {
        "sweeps": sweeps,
        "maxdiff": float(maxdiff),
        "converged": bool(converged),
        "omega": float(omega),
        "c_grid": float(c_grid),
        "c_ideal": float(c_ideal),
        "edge_excess": float(100 * (c_grid / c_ideal - 1)),
        "e_centre": float(Emag[g.centre, g.centre]),
        "e_ideal": float(g.v0 / g.gap),
        "e_tip": float(Emag[g.y0, g.x0]),
    }


if __name__ == "__main__":
    g = geometry()
    V, mx = solve(g)
    print(f"grid        {g.m} x {g.m}, h = {g.h:.4f} m, box {g.L} m")
    print(f"plates      {g.plate_len:.2f} m long, {g.gap:.2f} m apart, {g.v0:.0f} V")
    print(f"converged   {mx['converged']} in {mx['sweeps']} sweeps "
          f"(omega = {mx['omega']:.4f}, maxdiff {mx['maxdiff']:.2e})")
    print(f"E centre    {mx['e_centre']:.2f} V/m   ideal V/d {mx['e_ideal']:.2f}")
    print(f"E tip       {mx['e_tip']:.1f} V/m")
    print(f"C' grid     {mx['c_grid']:.4e} F/m")
    print(f"C' ideal    {mx['c_ideal']:.4e} F/m")
    print(f"edge excess {mx['edge_excess']:.1f} %")
```

- [ ] **Step 4: Run the tests to verify they pass**

```bash
cd /Users/sergiogutierrez/Documents/Projects/research/repos/2d-field-solvers/solvers
/usr/bin/python3 test_solvers.py
```

Expected: `15 passed`. The two refinement tests are the slow ones (m=401 takes a few seconds).

- [ ] **Step 5: Run the solver and capture its numbers for the README**

```bash
/usr/bin/python3 solvers/laplace.py | tee /tmp/laplace-run.txt
```

Expected, within rounding of the prototype measured 2026-10-10: 689 sweeps, E centre ≈ 100.25 V/m, E tip ≈ 183.8 V/m, C' grid ≈ 3.0565e-11 F/m, edge excess ≈ 72.6 %.

**This output is the source for the README's spec table in Task 7.** If it disagrees with the spec's §3 table by more than rounding, the README quotes *this* run, not the spec — the spec's numbers came from the patched 2017 C and this is the independent reproduction that gates phase 1.

- [ ] **Step 6: Stage and show. Do not commit.**

```bash
git add solvers/laplace.py solvers/test_solvers.py && git diff --cached --stat
```

Wait for approval before `git commit -m "feat: red-black SOR Laplace solver with capacitance by energy integration"`.

---

### Task 5: `solvers/plot.py` — FIG. 1

**Files:**
- Create: `solvers/plot.py`

**Interfaces:**
- Consumes: `grid.Grid`, `laplace.field`.
- Produces: `two_panel(V, g, path, title=None) -> None`, writing a PNG. Takes an array and a Grid and nothing else, so `heat.py` and `wave.py` reuse it unchanged in later phases.

- [ ] **Step 1: Write `solvers/plot.py`**

This keeps the two-panel layout Sergio approved on 2026-10-10: potential with equipotentials and field lines on the left, `|E|` on the right.

```python
"""
The shared figure. Takes a solved grid and renders it; knows nothing about
which equation produced it, so all three solvers use it unchanged.

Figures sit on light stock in both themes, per the Datasheet rules in
../../hardware/docs/DESIGN.md section 3.1.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from laplace import field

STOCK = "#F1F2EF"   # --figure
INK = "#141518"     # --ink


def two_panel(V, g, path, title=None):
    """Potential with field lines, beside field magnitude. Writes a PNG."""
    _, _, Emag = field(V, g.h)
    x = np.linspace(0, g.L, g.m)
    xi, yi = np.meshgrid(x, x)
    Ey, Ex = np.gradient(-V, g.h)
    vmax = abs(g.v0) / 2

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 6.2), layout="constrained")
    fig.patch.set_facecolor(STOCK)

    im = ax.imshow(V, extent=[0, g.L, g.L, 0], cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.contour(xi, yi, V, levels=np.linspace(-vmax, vmax, 19)[1:-1],
               colors=INK, linewidths=0.4, alpha=0.5)
    ax.streamplot(xi, yi, Ex, Ey, color=INK, density=1.4, linewidth=0.6, arrowsize=0.8)
    fig.colorbar(im, ax=ax, label="V  [volt]")
    ax.set_title("Potential + field lines")

    # Clip the colour scale below the tip singularity, which would otherwise
    # take the whole range and leave the interior flat.
    im2 = ax2.imshow(Emag, extent=[0, g.L, g.L, 0], cmap="inferno",
                     vmin=0, vmax=2.5 * g.v0 / g.gap)
    fig.colorbar(im2, ax=ax2, label="|E|  [V/m]")
    ax2.set_title("Field magnitude — edge enhancement at the plate tips")

    for a in (ax, ax2):
        a.set_facecolor(STOCK)
        a.set_xlim(0, g.L)
        a.set_ylim(g.L, 0)
        a.set_xlabel("x [m]")
        a.set_ylabel("y [m]")

    fig.suptitle(title or
                 f"Parallel-plate capacitor, 2D Laplace "
                 f"(SOR, {g.m}x{g.m}, h={g.h:.4f} m, box {g.L:.0f} m)")
    fig.savefig(path, dpi=130, facecolor=STOCK)
    plt.close(fig)
```

- [ ] **Step 2: Render FIG. 1 and look at it**

```bash
cd /Users/sergiogutierrez/Documents/Projects/research/repos/2d-field-solvers
mkdir -p figures
/usr/bin/python3 -c "
import sys; sys.path.insert(0, 'solvers')
from grid import geometry
from laplace import solve
from plot import two_panel
g = geometry(); V, mx = solve(g)
two_panel(V, g, 'figures/fig-01-field.png')
print('converged', mx['converged'], mx['sweeps'], 'sweeps')
"
open figures/fig-01-field.png
```

Expected: uniform field and linear ramp between the plates, lines bowing out at the edges, two bright tips in the right panel. If the right panel is flat black with two dots, the `vmax` clip is wrong — raise it until the interior is visible.

- [ ] **Step 3: Stage and show. Do not commit.**

```bash
git add solvers/plot.py figures/fig-01-field.png && git diff --cached --stat
```

Wait for approval before `git commit -m "feat: shared two-panel field figure"`.

---

### Task 6: The remaining figures

**Files:**
- Create: `solvers/figures.py`

**Interfaces:**
- Consumes: `grid.geometry`, `laplace.solve`, `laplace.field`, `grid.EPS0`.
- Produces: `main()` writing `figures/fig-02-edge-excess.png`, `figures/fig-03-charge-density.png`, `figures/fig-04-convergence.png`.

- [ ] **Step 1: Write `solvers/figures.py`**

```python
"""
FIG. 2-4 for the README. Run: /usr/bin/python3 solvers/figures.py

FIG. 2 is the argument: the textbook formula is a limit, and this geometry is
nowhere near it.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from grid import EPS0, geometry
from laplace import field, solve
from plot import INK, STOCK

# At L=5 with gap=1, plates longer than about 3.5 m come close enough to the
# grounded walls that the walls, not the physics, flatten the curve. Measured
# 2026-10-10: the ratio moves only 1.568 -> 1.542 from l/d 3 to 4.
MAX_RATIO_FOR_BOX = 3.5


def fig_edge_excess(path, m=151, L=5.0, gap=1.0):
    lengths = [x for x in np.arange(0.5, 4.01, 0.25) if x <= MAX_RATIO_FOR_BOX]
    ld, ratio = [], []
    for plate_len in lengths:
        g = geometry(m=m, L=L, plate_len=float(plate_len), gap=gap)
        _, mx = solve(g)
        assert mx["converged"], f"did not converge at plate_len={plate_len}"
        ld.append(g.plate_len / g.gap)
        ratio.append(mx["c_grid"] / mx["c_ideal"])

    fig, ax = plt.subplots(figsize=(7, 4.4), layout="constrained")
    fig.patch.set_facecolor(STOCK); ax.set_facecolor(STOCK)
    ax.plot(ld, ratio, color=INK, marker="o", ms=3.5, lw=1.2)
    ax.axhline(1.0, color=INK, ls="--", lw=0.8, alpha=0.6)
    ax.annotate("textbook limit, C' = eps0 l / d", xy=(ld[-1], 1.0),
                xytext=(0, 8), textcoords="offset points", ha="right",
                fontsize=9, color=INK)
    ax.set_xlabel("plate length / separation,  l / d")
    ax.set_ylabel("C' / C'ideal")
    ax.set_title(f"Edge excess falls as the plates get long  (box {L:.0f} m, {m}x{m})")
    fig.savefig(path, dpi=130, facecolor=STOCK); plt.close(fig)
    return ld, ratio


def fig_charge_density(path, m=251):
    g = geometry(m=m)
    V, _ = solve(g)
    _, Ey, _ = field(V, g.h)
    # sigma = eps0 * E_n just below the upper plate's inner face
    cols = np.arange(g.x0, g.x1 + 1)
    sigma = EPS0 * Ey[g.y0 + 1, cols]
    pos = (cols - g.centre) * g.h

    fig, ax = plt.subplots(figsize=(7, 4.4), layout="constrained")
    fig.patch.set_facecolor(STOCK); ax.set_facecolor(STOCK)
    ax.plot(pos, np.abs(sigma), color=INK, lw=1.3)
    ax.set_xlabel("position along the plate from its centre  [m]")
    ax.set_ylabel("|sigma|  [C/m^2]")
    ax.set_title("Surface charge density piles up at the tips")
    fig.savefig(path, dpi=130, facecolor=STOCK); plt.close(fig)
    return pos, sigma


def fig_convergence(path, sizes=(51, 101, 151, 201, 251)):
    jac, sor = [], []
    for m in sizes:
        _, mj = solve(geometry(m=m), omega=1.0)
        _, ms = solve(geometry(m=m))
        jac.append(mj["sweeps"]); sor.append(ms["sweeps"])

    fig, ax = plt.subplots(figsize=(7, 4.4), layout="constrained")
    fig.patch.set_facecolor(STOCK); ax.set_facecolor(STOCK)
    ax.plot(sizes, jac, color=INK, marker="o", ms=3.5, lw=1.2, label="Jacobi")
    ax.plot(sizes, sor, color=INK, marker="s", ms=3.5, lw=1.2, ls="--", label="SOR")
    ax.set_yscale("log")
    ax.set_xlabel("grid size m"); ax.set_ylabel("sweeps to 1e-6")
    ax.set_title("Sweeps to converge: Jacobi O(m^2) against SOR O(m)")
    ax.legend(frameon=False)
    fig.savefig(path, dpi=130, facecolor=STOCK); plt.close(fig)
    return sizes, jac, sor


def main():
    ld, ratio = fig_edge_excess("figures/fig-02-edge-excess.png")
    print("FIG 2  l/d", [f"{v:.2f}" for v in ld])
    print("       C'/C'ideal", [f"{v:.3f}" for v in ratio])
    pos, sigma = fig_charge_density("figures/fig-03-charge-density.png")
    print(f"FIG 3  |sigma| centre {abs(sigma[len(sigma)//2]):.3e}  "
          f"tip {abs(sigma[0]):.3e} C/m^2  "
          f"ratio {abs(sigma[0] / sigma[len(sigma)//2]):.1f}x")
    sizes, jac, sor = fig_convergence("figures/fig-04-convergence.png")
    print("FIG 4  m", list(sizes)); print("       jacobi", jac); print("       sor", sor)


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it and record the printed numbers**

```bash
cd /Users/sergiogutierrez/Documents/Projects/research/repos/2d-field-solvers
/usr/bin/python3 solvers/figures.py | tee /tmp/figures-run.txt
```

Expected: the FIG. 2 ratios decrease monotonically from about 2.9 toward 1.5; the FIG. 3 tip/centre ratio is several times 1; FIG. 4's Jacobi column grows far faster than the SOR column. These printed numbers go in the README in Task 7.

- [ ] **Step 3: Look at all three figures**

```bash
open figures/fig-02-edge-excess.png figures/fig-03-charge-density.png figures/fig-04-convergence.png
```

Expected: FIG. 2 decays toward the dashed line at 1.0 without touching it. If it dips below 1.0, the box is too small for the longest plate — lower `MAX_RATIO_FOR_BOX` until it does not.

- [ ] **Step 4: Stage and show. Do not commit.**

```bash
git add solvers/figures.py figures/ && git diff --cached --stat
```

Wait for approval before `git commit -m "feat: edge-excess, charge-density and convergence figures"`.

---

### Task 7: The datasheet README

**Files:**
- Modify: `README.md` (full rewrite)

**Interfaces:**
- Consumes: `/tmp/laplace-run.txt` (Task 4) and `/tmp/figures-run.txt` (Task 6). **Every number in the README is copied from those two files.** No number is carried over from the spec, which measured the patched 2017 C rather than this solver.

- [ ] **Step 1: Re-read the house style before writing**

```bash
sed -n '/^### 3.1/,/^### 3.2/p' ../../../hardware/docs/DESIGN.md
sed -n '1,60p' ../../../website/src/content/prototypes/wf-01.mdx
```

The rules that matter: numbered sections `## 01 TITLE`; spec tables Parameter · Value · Unit · Source with numbers right-aligned; figures captioned `FIG. n · what it is`, numbered in reading order, caption ending with its source; present tense, definite, written for an outside reader; every value carries a unit; calculated values say `Calculated`, measured ones say `Measured`.

- [ ] **Step 2: Write `README.md`**

Structure, with the numbers filled in from the two run files:

```markdown
# 2D Field Solvers

FINITE DIFFERENCE · PYTHON + NUMPY · ELECTROSTATICS · REV 0.1

Finite-difference field solvers on a shared 2D grid. `laplace.py` solves the
electrostatic field of a parallel-plate capacitor; `heat.py` and `wave.py`
follow, over the same geometry and the same figure.

## 01 Description

[Two paragraphs. First: a conductor at fixed potential is a boundary condition
and the field follows from it; relaxation finds it. Second: the result is the
textbook picture with the textbook's lie removed — uniform field and linear
ramp between the plates, and at the tips a field <E_tip/E_centre> times higher
with the lines bowing out into space the formula says is empty.]

<FIG. 1>

### Measured

| Parameter | Value | Unit | Source |
|---|---:|---|---|
| Grid | 251 x 251 | cell | h = <h> m |
| Box | 5.00 | m | |
| Plate length | <plate_len> | m | Snapped to cells |
| Plate separation | <gap> | m | Snapped to cells |
| Potential difference | 100 | V | |
| Field, gap centre | <e_centre> | V/m | Measured |
| Field, ideal V/d | <e_ideal> | V/m | Calculated |
| Field, plate tip | <e_tip> | V/m | Measured |
| C', grid | <c_grid> | F/m | Measured |
| C', ideal eps0 l/d | <c_ideal> | F/m | Calculated |
| Edge excess | <edge_excess> | % | Calculated |
| Sweeps to 1e-6 | <sweeps> | — | Measured, SOR |

## 02 Method

[The 5-point stencil. The Dirichlet mask: plates and walls never updated.
Red-black SOR and why omega is computed from m rather than hard-coded. The
one-paragraph statement of the shared-stencil idea from spec section 2.1, so a
reader meets it here and not only in the spec.]

## 03 Results

[FIG. 2 and the l/d argument: the excess is not an error in the simulation, it
is what the formula leaves out, and it decays toward 1 as the plates get long.
Then FIG. 3, charge density piling up at the tips, which is *why* the tip field
is high. Then FIG. 4, the numerics, as supporting material.]

**Caveat.** C' is computed by integrating the field energy over the whole box,
whose grounded walls truncate the fringing field. C' therefore depends on the
box size, and every result above is quoted with its box. The trend in FIG. 2,
not the single number, is the result.

## 04 Running it

[Commands with the interpreter named, expected runtime from the measured runs.]

## 05 History

[The 2017 Metodos Computacionales assignment at Universidad de los Andes:
solve the capacitor, then parallelise it with MPI. legacy/ holds that code
unchanged. The repo was called 2d-heat-equation-solver until 2026-10-10, which
was not wrong so much as early — see section 06 and legacy/BUGS.md.]

## 06 Known bugs in the 2017 code

[The seven, one short paragraph each, linking to legacy/BUGS.md for the fixes.
State plainly that the serial version worked and the MPI submission could not
have.]
```

- [ ] **Step 3: Check every number against the run files**

```bash
cd /Users/sergiogutierrez/Documents/Projects/research/repos/2d-field-solvers
cat /tmp/laplace-run.txt /tmp/figures-run.txt
grep -nE '[0-9]' README.md | head -40
```

Every figure in the README must appear in one of the two run files or be a stated input (5.00 m box, 100 V). Anything else is invented — remove it or re-derive it.

- [ ] **Step 4: Check the voice**

Read the README once as an outside reader. It fails if it contains: "we", "I then", "as can be seen", "the assignment asked", a value without a unit, or a figure without a numbered caption and a source.

- [ ] **Step 5: Stage and show. Do not commit.**

```bash
git add README.md && git diff --cached --stat
```

Show Sergio the rendered README and wait before `git commit -m "docs: datasheet README"`.

- [ ] **Step 6: Final check of the whole repo**

```bash
/usr/bin/python3 solvers/test_solvers.py
git status --short
ls figures legacy solvers
```

Expected: all tests pass; nothing untracked but intended files; `legacy/` holds the 2017 tree plus `BUGS.md`.

---

## What this plan does not cover

- **The website.** `src/lib/capacitor.ts`, `src/components/sheet/capacitor-field.tsx`, the WR-04 post and the two one-line `/research` changes get their own plan, written once this one lands. That plan must also run `npm run check`, the public-copy grep gate from the website's `AGENTS.md`, and `phone-check` with the new post added to the page list.
- **`heat.py`** (spec phase 2) and **`wave.py`** (spec phase 3). Each gets its own design pass.
- **Pushing anything.** The repo is public and `main` is its default branch. Nothing is pushed without Sergio saying so.
