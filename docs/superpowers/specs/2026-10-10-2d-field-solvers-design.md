# 2D Field Solvers — design

Date: 2026-10-10
Status: approved in chat 2026-10-10, pending written review
Repo: `Projects/research/repos/electrostatics-2d-plates` (to be renamed `2d-field-solvers`)

## 1. Intent

Sergio wrote a parallel-plate capacitor solver as a Métodos Computacionales
assignment at Universidad de los Andes in 2017: relax Laplace's equation on a
2D grid in C, then parallelize it with MPI. He remembered it as "the electric
field between two conductive plates, and you could see the field around, inside
and around the surfaces". It was found on 2026-10-10 in the GitHub repo
`ChiefGuti/2d-heat-equation-solver` — misnamed, since the code solves
electrostatics, not heat.

What he asked for, in his words: documentation following his website, going
over the repo, bringing the plots to "more modern and useful data", and making
it dynamic if possible. On the figure's subject he said: "go for the 2D wave
equation but start with this and then do a second implementation. I like the
plots."

**Intended outcome.** One repo holding two finite-difference field solvers over
the same geometry — electrostatics now, the wave equation next — documented in
the Datasheet language his website uses, with an interactive figure published
as a writing post on sergiosgutierrez.com.

**Who it is for.** An outside reader arriving from his site who wants to see
that he can take a physics problem from boundary conditions to a number. The
2017 code is not that; it is an artifact that proves where the interest came
from.

**Success criteria.**
1. `python3 solvers/laplace.py` runs on a clean machine and reproduces the
   measured table in the README.
2. The README reads as a datasheet, not a lab report: definite, present tense,
   every value carrying a unit and a source.
3. The published figure answers a question the static PNG cannot — moving a
   slider changes the geometry and the reported capacitance follows.
4. The 2017 C code survives unchanged, with its bugs documented rather than
   quietly corrected.

**Explicitly not in scope.** 3D. Dielectric materials. Charge transport. A
general PDE framework. Anything on `/research`, which is the two-thesis
registry.

## 2. Decisions taken

| # | Decision | Rationale |
|---|---|---|
| D1 | **No new part-number series.** The post is `WR-04` | Sergio, 2026-10-10: "we should have less shorthands maybe so they are simpler to categorize". The site has six series already (SG, LAB, TH, WF, WR, PM) and a seventh earns nothing. The solvers are named in prose — `laplace`, `heat`, `wave` — not coded. WR-03 is taken by `msc-rerun-2026`, so the post is WR-04. |
| D2 | Rename repo to `2d-field-solvers` | Approved by Sergio 2026-10-10. `2d-heat-equation-solver` described one equation the repo did not yet contain; the new name covers all three. |
| D3 | Phase 1 subject = edge effects and capacitance | The honest "useful data". The 2017 code was one plot from it. Numerics story stays as a second figure, not the headline. |
| D4 | Python + NumPy for the solver | Matches his stack. The C stays in `legacy/`. |
| D5 | Physics duplicated in TypeScript for the web | Follows the site's existing pattern: `hom.ts` and `nlpn.ts` already restate thesis physics as illustrative models with tests. |
| D6 | Red-black SOR, not Jacobi | Measured 2026-10-10: at 96×96, SOR converges in 265 sweeps against Jacobi's 3,057 — 11.5×. This is what makes the browser figure interactive rather than animated-once. |
| D7 | Published as a `/writing` post, surfaced on `/research` | Sergio: "this is part of my research can I connect it there". `writing.ts` + MDX + `src/components/sheet/` already does everything needed, and `postsAbout()` already backlinks posts onto pages. One line on the hub connects it (§8). A new route or content type is unearned. |
| D8 | **Three solvers, one grid: Laplace, heat, wave** | Sergio: "maybe I want to continue to heat somehow?" The same 5-point stencil runs all three — they differ only in how many time derivatives they carry (§2.1). This makes the repo's subject the method, not one problem, and retires the question of whether the old name was wrong. |

### 2.1 One stencil, three equations

The repo's organising idea, and the reason the three solvers belong together
rather than in three repos. Every one of them evaluates the same quantity at
every interior cell — the Laplacian, as the four neighbours minus four times
the centre:

```
∇²u ≈ (u[i-1,j] + u[i+1,j] + u[i,j-1] + u[i,j+1] − 4u[i,j]) / h²
```

What differs is only how many time derivatives sit on the left:

| Solver | Equation | Order in t | What it does | Step rule |
|---|---|---|---|---|
| `laplace.py` | `∇²V = 0` | none | Settles to equilibrium; the steady state the other two end at | Relax: `u ← u + ω(avg − u)` |
| `heat.py` | `∂u/∂t = α∇²u` | first | Diffuses; sharp edges smear and never come back | `u += α·Δt·∇²u`, `Δt ≤ h²/4α` |
| `wave.py` | `∂²u/∂t² = c²∇²u` | second | Propagates and reflects; energy is conserved, so it rings | Leapfrog, `c·Δt/h ≤ 1/√2` |

Two consequences worth stating, because they are what make this worth building
as one thing:

1. **Laplace is the limit of heat.** Run `heat.py` long enough with fixed
   boundaries and it converges to exactly what `laplace.py` returns. That is a
   test, not a remark — the two solvers check each other, and it is the cleanest
   demonstration that Jacobi relaxation *is* a diffusion process in disguise.
   The 2017 assignment solved a capacitor with a heat method and nobody said so.
2. **The same geometry means the same figure.** A conductor at fixed voltage, a
   boundary at fixed temperature and a wall that reflects are all Dirichlet
   conditions on the same mask array. `plot.py` does not need to know which
   physics produced the grid it is handed.

The old repo name was not wrong so much as premature.

## 3. The headline result

Computed from the converged 250×250 field on 2026-10-10, by energy integration
`W = ½ε₀∫|E|²dA` and `C' = 2W/V²`:

| Parameter | Value | Unit | Source |
|---|---:|---|---|
| C′, grid | 3.058 × 10⁻¹¹ | F/m | Measured |
| C′, ideal ε₀l/d | 1.771 × 10⁻¹¹ | F/m | Calculated |
| Edge excess | 72.7 | % | Calculated |
| Field, gap centre | 99.8 | V/m | Measured |
| Field, ideal V/d | 100.0 | V/m | Calculated |
| Field, plate tip | 183.3 | V/m | Measured |

The textbook formula is wrong by 73% for this geometry, because the plates are
only twice as long as they are far apart. That is the post's argument: the
uniform-field picture is a limit, and this geometry is nowhere near it.

These numbers came from the **patched 2017 C code** (§5). They are the target
`laplace.py` must reproduce, not a result it may assume: phase 1 is not done
until the Python solver produces them independently and the README table is
rewritten from its output.

**Caveat to carry into the docs.** The energy integral runs over the whole box,
whose grounded walls truncate the fringing field, so C′ depends on the box size.
The docs must state the box dimensions with the result, and the figure must
sweep l/d to show the excess falling as the plates get long — approaching the
textbook limit rather than contradicting it.

## 4. Repo structure

```
README.md                  the datasheet (§5)
solvers/grid.py            the shared piece: geometry, the Dirichlet mask, laplacian()
solvers/laplace.py         phase 1: red-black SOR to equilibrium
solvers/heat.py            phase 2: explicit diffusion, converges to laplace.py
solvers/wave.py            phase 3: leapfrog, reflects off the conductors
solvers/plot.py            the two-panel figure, shared by all three
solvers/test_solvers.py    the check (§7)
docs/superpowers/specs/    this document
legacy/
  BUGS.md                  the seven 2017 defects, each with its fix
  t1/                      the 2017 C and MPI code, unchanged
```

Units: SI throughout, metres and volts. `solvers/` is a plain directory of
scripts, not a package — there is one import boundary (`plot.py` takes an array)
and no more.

**Deleted** (confirmed by Sergio 2026-10-10): `t1/Entrega/a.out`, a compiled
binary committed in 2017 that no longer runs on his machine, and
`t1/sinparalelisar/output.txt`, 608 KB of generated grid data. A `.gitignore`
covers `a.out`, `*.o`, compiled solver binaries, `output.txt` and `__pycache__`
so neither returns. `a.out` remains in git history; removing it from there
needs a force-push and is not worth it for 16 KB.

### Interfaces

`laplace.py` exposes one function:

```python
def solve(m=250, L=5.0, plate_len=2.0, gap=1.0, v0=100.0, tol=1e-6)
    -> tuple[np.ndarray, dict]
```

Returns the potential grid and a metrics dict (`sweeps`, `maxdiff`, `c_grid`,
`c_ideal`, `edge_excess`, `e_centre`, `e_tip`). `plot.py` takes that pair and
nothing else, so the wave solver can reuse it by returning the same shape.

## 5. Documentation

The README becomes a datasheet, following `hardware/docs/DESIGN.md` §3.1 and
the WF-01 sheet:

- Header: title, fields joined by ` · ` in upper case, `REV`.
- Numbered sections `## 01 Description`, `## 02 Method`, …
- Spec tables: Parameter · Value · Unit · Source, numbers right-aligned.
- Figures captioned `FIG. n · what it is`, numbered in reading order, ending
  with their source.
- Voice: outside reader, present tense, definite. No "we then observe that".
- Every value carries a unit; calculated values say so.

Section plan:

| § | Content |
|---|---|
| 01 | Description — what the field looks like and the 73% result |
| 02 | Method — relaxation, the boundary conditions, red-black SOR, why ω matters |
| 03 | Results — the measured table, the two-panel figure, the l/d sweep |
| 04 | Running it — commands, arguments, runtime |
| 05 | History — the 2017 assignment, and that `legacy/` is kept unchanged |
| 06 | Known bugs in the 2017 code — the seven below; the MPI submission never compiled |

`legacy/BUGS.md` records what was found on 2026-10-10. **Each entry carries its
fix, not just its diagnosis** (Sergio, 2026-10-10: the documentation should
maintain the solution and the bug fix as part of the documentation). The
format per bug is four parts:

1. **Symptom** — what the output did.
2. **Cause** — the line, quoted.
3. **Fix** — the corrected line, as a diff hunk.
4. **Why it matters** — the physics or C behaviour it turns on.

So the patch is not thrown away when the C moves to `legacy/`: the source files
stay exactly as submitted in 2017, and the fixes live in `BUGS.md` as a
readable record beside them. A reader sees both the mistake and the correction
without the artifact being quietly rewritten. This supersedes the earlier
decision to simply revert the patch.

Seven, found 2026-10-10 and 2026-10-11. **The headline is that the parallel
submission never compiled**, which is a stronger and simpler statement than the
earlier draft of this spec made.

Note there are **two** MPI files and they are not the same. `t1/Entrega.c` has
`h = 0.02`; `t1/Entrega/Entrega.c`, the submission, has `h = 5/512`. An earlier
draft of this spec conflated them.

In the MPI submission (`Entrega/Entrega.c`), verified by syntax-checking against
a stubbed `mpi.h`:

| # | Line | Code | Effect |
|---|---|---|---|
| 1 | 160 | `V[i]=V1[i]` | No semicolon. Hard compile error — the deliverable never built |
| 2 | 162 | `n=+1;` | Parses as `n = +1`, not `n += 1`. The counter sticks at 1 and `while(n<N)` never terminates |
| 3 | 10 | `double h = 5/512;` | Both operands `int`, so the division happens in integer arithmetic and yields 0 before the conversion to `double`. Every plate coordinate then divides by zero |
| 4 | 132 | `if(world_rank==world_size)` | Ranks run 0..`world_size-1`, so never true. The last rank's branch is dead and its slab is never relaxed |
| 5 | 85 | `if (!(world_rank==0 && world_rank==world_size))` | A rank cannot be both, so always true. Rank 0 runs its own branch and then this one, over the wrong range |

In the serial version (`sinparalelisar/placas.c`), which did work:

| # | Code | Effect |
|---|---|---|
| 6 | `malloc` instead of `calloc` | The sweep never writes row/col 0 and m−1. Those cells are the grounded box and the solution depends on them being zero; `malloc` leaves whatever the allocator last held, and relaxation propagates it inward |
| 7 | `V_new = init(..., V)` | `init` returns the pointer it is given, so both names aliased one array. The copy-back loop assigned each cell to itself, silently turning the intended Jacobi iteration into Gauss–Seidel — which converges about twice as fast, so the bug helped |

**Not a bug, and `BUGS.md` must say so.** The MPI version has no plate guard in
its sweep, which looks like a defect until you see that `init()` is called at
the top of every iteration, with the comment *"cada procesador restablece las
condiciones de las varillas despues de cada iteracion"*. Re-stamping the
Dirichlet condition each sweep instead of masking cells is the better design
once the domain is cut into slabs: no per-cell ownership test, same answer. An
earlier pass of this work nearly recorded it as bug 8.

The working tree was patched on 2026-10-10 to produce a figure (calloc, the
pointer, and a convergence tolerance replacing the fixed 125,000 sweeps). The
`legacy/` sources are restored to the 2017 submission, and that patch is
preserved as the Fix sections of `BUGS.md` — recorded, not discarded.

§06 of the README summarises these for a reader who never opens `legacy/`, and
links to it.

## 6. The figures

The two-panel layout Sergio approved is kept: potential with equipotentials and
field-line streamplot on the left, `|E|` magnitude on the right showing tip
enhancement. What changes is what surrounds it.

- **FIG. 1** the two panels, as now, with the geometry stated in the caption.
- **FIG. 2** C′/C′_ideal against l/d, from about 0.5 to 20, showing the excess
  decaying toward 1 — the textbook limit earned rather than assumed.
- **FIG. 3** surface charge density σ = ε₀E_n along a plate, which diverges at
  the tips. This is why the tip field is 1.8× the interior field.
- **FIG. 4** (numerics, kept as supporting material) sweeps to converge against
  grid size for Jacobi and SOR, showing the O(m²) against O(m) separation.

Colours follow §3: figures on light stock in both themes, `--spi` blue for
focus, `--est` amber for calculated-not-measured, nothing decorative.

## 7. Testing

One file, `solvers/test_solvers.py`, plain asserts, runnable with `python3`:

- The gap-centre field equals V/d within 0.5%.
- The potential is antisymmetric about the midplane within 1e-9.
- Grounded walls stay at exactly 0.
- SOR and Jacobi converge to the same field within tolerance — the check that
  the fast path did not change the physics.
- `c_grid > c_ideal` always, and the ratio decreases monotonically as l/d grows.
- Halving h moves C′ by less than 1% — grid independence.

The last two are the ones that would actually catch a wrong answer; the rest
are cheap.

## 8. Website

No new dependencies. Follows the existing figure pattern exactly.

```
src/lib/capacitor.ts        pure model: solveSOR(), capacitance(), fieldAt()
src/lib/capacitor.test.ts   same assertions as §7, ported
src/components/sheet/capacitor-field.tsx   client component
src/content/writing.ts      the WR-04 row
src/content/writing/the-capacitor-the-textbook-hides.mdx   the post
public/writing/the-capacitor-the-textbook-hides/           static fallback figures
src/app/research/page.tsx   one line, below
```

### Connecting it to /research

Verified in the code on 2026-10-10, not assumed. `postsAbout(href)` in
`writing.ts:85` returns every post whose `related` array points at `href`, and
`/research` already calls it — but only per thesis
(`research/page.tsx:54`, `postsAbout('/research/${t.slug}')`), so a post about
the capacitor would not appear on the hub on its own.

Two changes, both one line:

1. The post's `related` includes `{ href: "/research", label: "SG-04 · Research" }`.
2. `/research`'s `SeeAlso` gains `...postsAbout("/research")` beside the
   per-thesis spread it already has.

That is the whole connection. No new route, no new content type, no new code
series, and the mechanism is the one the site already uses for theses,
datasheets and prototypes. Any later research-code post joins the hub by
pointing at it.

The post itself links back to the repo, and `person.ts` `software` gains the
repo row so it also appears on `/projects`.

`capacitor-field.tsx` draws plain SVG with Tailwind, ink `#141518` / muted
`#5e625c` / rule `#C6C9C4`, sliders in the style of `hom-dip.tsx`. Grid 96×96,
red-black SOR over a flat `Float64Array`. Measured 2026-10-10: 265 sweeps,
~2.4M cell updates, a few milliseconds — so the figure **re-solves on every
slider change**, not just once. A "Show the numbers" table sits under it per
§3.1, reporting C′, C′_ideal and the excess for the current geometry.

Accessibility and fallback: sliders are real `<input type="range">` with
labels; the component renders a static converged field on first paint so a
no-JS or headless reader still sees FIG. 1. `phone-check` must pass, and its
headless Chrome has no WebGL — irrelevant here, since this is SVG, not three.js.

Before any commit touching the website: `npm run check`, the public-copy grep
gate from `AGENTS.md`, and the phone-check page list extended with the new post.

## 9. Phasing

**Phase 1 — Laplace.** Repo restructure and rename, `legacy/` with BUGS.md,
`grid.py`, `laplace.py`, `plot.py`, tests, the datasheet README, then the
website model, component, the WR-04 post and the `/research` connection.

**Phase 2 — heat.** `∂u/∂t = α∇²u` on the same grid. Small on its own, because
`grid.py` and `plot.py` already exist, and it buys the test that matters: run it
to steady state and it must reproduce `laplace.py` cell for cell. It also makes
the repo's old name true and gives the browser figure its first animation, since
diffusion is one array update per frame.

**Phase 3 — wave.** `∂²u/∂t² = c²∇²u`, leapfrog, conductors as reflecting
boundaries. The one Sergio asked for first. Cheap per frame, so it animates
natively. Its CFL condition (`c·Δt/h ≤ 1/√2`) is the thing most likely to go
wrong and gets its own test: step past the limit and the energy must blow up,
step under it and the total energy must stay flat.

Phases 2 and 3 are specified at this level deliberately. Each gets its own
design pass once phase 1 is standing.

## 10. Risks

| Risk | Handling |
|---|---|
| C′ depends on box size, so the 72.7% is geometry-specific | State the box with every result; FIG. 2 sweeps l/d so the trend, not the single number, carries the argument |
| Duplicated physics in Python and TypeScript drifts | Same assertions in both test files; the post states the web model is illustrative, as `nlpn.ts` does |
| SOR's ω is grid-size dependent | `ω = 2/(1+sin(π/m))` computed from m, not hard-coded; the Jacobi-agreement test catches a bad ω |
| The site is live and a push to main deploys | Nothing pushed without explicit approval; `npm run check` and the public-copy gate before any commit |

## 11. Open questions for Sergio

Resolved 2026-10-10: no new code series (D1), rename approved (D2), heat added
as phase 2 (D8), connected to `/research` through `postsAbout` (§8).

One left:

1. **Repo visibility.** It is public today, and the plan puts the seven 2017
   bugs in a public `BUGS.md`. Finding them is the best part of the work, but
   it is a public page in his name saying the submission was broken. If he does
   not want to decide, everything builds locally and nothing is pushed —
   `legacy/BUGS.md` is the only piece that would need to change, and it changes
   by deletion.
