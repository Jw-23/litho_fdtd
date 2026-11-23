# litho_fdtd

`litho_fdtd` is a C++23 three-dimensional finite-difference time-domain (FDTD) solver for
photolithography mask spectra. It uses the geometry types from `litho_geometry`, adds a full
Yee-grid Maxwell solver, and returns complex diffraction orders suitable for imaging or
source-mask optimization models.

The public namespace is `litho::fdtd`. The library follows the same conventions as
`litho_geometry`: value types, free functions, `std::expected` error returns, no explicit
`throw`, and English source comments.

## Features

- 3D Cartesian Yee grid with second-order interior stencils.
- Convolutional PML (CPML) on any non-periodic axis.
- Periodic and Bloch-periodic lateral boundaries for mask periods and oblique illumination.
- Pulsed plane-wave illumination initialized in a uniform layer, with correct `E`/`H` phase and
  envelope-derivative terms for both s and p polarization.
- Isotropic materials from a complex refractive index at a reference wavelength, represented as
  a lossy dielectric or a single Drude pole for negative real permittivity.
- DFT plane monitors and a mask-spectrum post-processor: complex orders, actual-field Poynting
  efficiencies, and incident-normalized spectra.
- 3D mask geometry through `litho_geometry::common::extrusion` and `box3` solids, voxelized one
  z slice at a time.
- Structure-of-arrays storage, fused H/E scheduler phases, and a chunk scheduler from the sibling
  `parallel` project; field kernels use at most four workers.

## Build

Requirements: CMake >= 3.25, Ninja (recommended), a C++23 compiler with `std::expected`,
Boost headers, and the sibling `parallel` scheduler library. `litho_fdtd` adds `../parallel`
automatically when the `parallel` target is not already defined.

```sh
cmake --workflow --preset release
./build/release/tests/litho_fdtd_tests
./build/release/examples/litho_fdtd_mask_example
```

Or step by step:

```sh
cmake --preset release
cmake --build --preset release
ctest --preset release
```

`litho_fdtd` adds the sibling `litho_geometry` project automatically when
`litho_geometry::litho_geometry` is not already defined. If you vendor both projects, use
`add_subdirectory(path/to/litho_fdtd)` and link `litho_fdtd::litho_fdtd`.

To build the throughput benchmark:

```sh
cmake --preset release -DLITHO_FDTD_BUILD_BENCHMARKS=ON
cmake --build --preset release
./build/release/benchmarks/litho_fdtd_benchmark [workers]
```

The benchmark reports million cell-updates per second for a 120×120×240 vacuum grid; the optional
`workers` argument follows `run_options::workers`. An optional second argument selects `monitors`
(three four-frequency monitors), `patterned` (a lossy magnetic dielectric pattern), or `oblique`
(the same pattern with complex Bloch fields):

```sh
./build/release/benchmarks/litho_fdtd_benchmark 4 patterned
./build/release/benchmarks/litho_fdtd_benchmark 4 oblique
# Empty trace directory, then the scheduler grain in rows:
./build/release/benchmarks/litho_fdtd_benchmark 4 patterned "" 128
```

## Quick start: a 193 nm line/space mask

```cpp
#include <litho_fdtd/fdtd.hpp>

using namespace litho::fdtd;
namespace common = litho::geometry::common;

mask_problem p;
p.grid.cells = {80, 4, 120};
p.grid.spacing = {5e-9, 5e-9, 5e-9};

// add() returns std::expected<size_t>; the quick start assumes success.
(void)p.materials.add(material{"air", 1.0, 0.0, 1.0, std::nullopt});
(void)p.materials.add(*from_refractive_index("fused silica", {1.56, 0.0}, 193e-9));
(void)p.materials.add(*from_refractive_index("MoSi", {2.5, 0.6}, 193e-9));

p.mask.background = 0;
p.mask.layers.push_back(planar_layer{0.0, 300e-9, 1}); // silica

auto line = common::make_extrusion(
    *common::rectangle(common::box<double, 2>{{0.0, 0.0}, {200e-9, 20e-9}}),
    300e-9, 360e-9);
p.mask.solids.push_back(patterned_solid{*line, 2});

p.source.z = 120e-9;         // initial pulse centre inside the silica
p.source.wavelength = 193e-9;
p.source.polarization = 0;   // 0 = s (TE), pi/2 = p (TM)
p.source.pulse_fwhm = 0.5;   // optical periods
p.run.duration = 20;         // optical periods
p.run.workers = 0;           // 0 = up to four workers (no CPU pinning)
p.run.chunk_rows = 128;      // scheduler grain in rows

p.monitor_z = 450e-9;        // in air, above the absorber
p.spectrum.min_order_x = -6;
p.spectrum.max_order_x = 6;
p.spectrum.background_index = 1.0;

auto result = compute_mask_spectrum(p);
if (!result) { /* result.error(), message(result.error()) */ }

for (const auto& order : result->spectrum.orders) {
    // order.field is {Ex, Ey, Ez}; order.efficiency is set when a reference run exists.
    // result->normalized contains orders divided by the reference (0,0) field.
}
```

The example program [`examples/mask_spectrum.cpp`](examples/mask_spectrum.cpp) prints the
orders and writes `litho_fdtd_mask_spectrum.csv`.

## Geometry

The 3D definitions live in `litho_geometry/common/solid.hpp` and are re-exported through
`<litho_geometry/common.hpp>`:

- `common::extrusion<T>`: an xy polygon extruded from `z_min` to `z_max`, with an optional
  uniform `top_scale` for tapered sidewalls.
- `common::box3<T>` / `box3d`: an axis-aligned 3D box.
- `common::solid3<T>` / `solid3d`: `std::variant<extrusion<T>, box3<T>>`.
- `common::section_at(solid, z)`: exact 2D cross-section, used by the voxelizer.
- `common::bounds(solid)`: 3D bounding box.

The solver rasterizes one z slice at a time; it never materializes a full voxel volume in the
geometry layer. Solids override planar layers, and later solids override earlier ones.

## API overview

| Header | Contents |
| --- | --- |
| `<litho_fdtd/grid.hpp>` | `grid`, `axis`, `component`, Yee offsets |
| `<litho_fdtd/material.hpp>` | `material`, `material_library`, `from_refractive_index`, `from_permittivity` |
| `<litho_fdtd/structure.hpp>` | `planar_layer`, `patterned_solid`, `structure`, `voxelize` |
| `<litho_fdtd/boundary.hpp>` | `boundaries`, `boundary_kind`, `cpml_options` |
| `<litho_fdtd/source.hpp>` | `plane_wave_source` |
| `<litho_fdtd/monitor.hpp>` | `monitor`, `monitor_result`, DFT plane accumulation |
| `<litho_fdtd/spectrum.hpp>` | `spectrum_options`, `diffraction_order`, `compute_spectrum`, `normalize_spectrum` |
| `<litho_fdtd/solver.hpp>` | `problem`, `run_options`, `run_result`, `solve` |
| `<litho_fdtd/mask.hpp>` | `mask_problem`, `mask_spectrum`, `compute_mask_spectrum` |
| `<litho_fdtd/fdtd.hpp>` | umbrella header |

All fallible entry points return `litho::fdtd::result<T>` (`std::expected<T, errc>`). The
implementation catches allocation failures at public boundaries and translates them to
`errc::allocation_failure`; user code never needs exceptions.

## Units and conventions

- SI units internally: `E` in V/m, `H` in A/m, lengths in metres, time in seconds.
- Time dependence is `exp(-i omega t)`. A complex refractive index is entered as
  `n = n' + i n''` with `n'' >= 0` for a passive material.
- `from_refractive_index` builds either:
  - a non-dispersive lossy dielectric when `Re(n^2) > 0`, with
    `epsilon_r = Re(n^2)` and `sigma = omega * eps0 * Im(n^2)`; or
  - a single Drude pole when `Re(n^2) < 0`, fitted to match `n^2` exactly at the reference
    wavelength. This keeps the time-domain update passive and stable while representing
    metals and negative-permittivity absorbers.
- Angles are radians. `theta` is the incidence angle from `+z`; `phi` is the azimuth in the xy
  plane. `polarization = 0` is s (TE), `pi/2` is p (TM). `from_below` selects +z or -z
  propagation.
- `plane_wave_source::z` is the centre of the initial pulse, not a current sheet. The pulse is
  a one-way wave packet, so there is no backward-propagating source wave.

## Numerics

### Yee grid

All six components share one array shape and are staggered on the lower edges/faces of each
cell:

| Component | Position |
| --- | --- |
| `Ex` | `(i+1/2, j, k)` |
| `Ey` | `(i, j+1/2, k)` |
| `Ez` | `(i, j, k+1/2)` |
| `Hx` | `(i, j+1/2, k+1/2)` |
| `Hy` | `(i+1/2, j, k+1/2)` |
| `Hz` | `(i+1/2, j+1/2, k)` |

E updates use backward differences of H; H updates use forward differences of E. At non-periodic
outer boundaries the outermost stencil is one-sided and CPML absorbs the field before it reaches
the boundary.

### Time step

`dt = courant / (c_max * sqrt(1/dx^2 + 1/dy^2 + 1/dz^2))`, where `c_max` is the largest
material phase velocity estimated from the minimum `epsilon_inf * mu_r`. The default Courant
number is 0.95. Drude materials are additionally checked against `omega_p * dt <= 1`.

### CPML

Each non-periodic axis gets a polynomial-graded CPML with independent integer and half-integer
coefficient tables, matching the staggered field locations. The default parameters are
`layers = 10`, `reflection = 1e-6`, `order = 3`, `kappa_max = 5`, and
`alpha_max = 0.05 * sigma_max`. For optical masks, use enough layers that the PML thickness is
roughly half a wavelength or more; 20–40 layers at 5–10 nm spacing is typical. The PML is
implementation-transparent: it is not a scattering object and needs no special handling in the
geometry.

### Periodic and Bloch boundaries

Axes with `boundary_kind::periodic` wrap the update stencil. For an oblique incidence angle the
solver automatically switches to complex analytic fields and applies the Bloch phase
`exp(i k_parallel L)` at the wrap. Real fields are used for normal incidence or when the lateral
boundaries are CPML.

### Source normalization

The pulse is initialized as a Gaussian plane-wave packet at `t = 0` (E) and `t = -dt/2` (H). The
H field includes the envelope-derivative term, so the initial condition satisfies Maxwell's
equations for a pulsed plane wave rather than only for a monochromatic wave. Monitor DFT sums are
divided by the analytic Gaussian envelope integral, so a unit-amplitude plane wave produces unit
complex spectrum coefficients. A monitor time window can reject late CPML or chamber echoes; in
that case the normalization is no longer the full pulse integral, so compare full and reference
runs with the same window.

### Mask spectrum

For a z-normal monitor, the Fourier order convention is

```text
kx = kx_incident + 2*pi*m/period_x
ky = ky_incident + 2*pi*n/period_y
kz = sqrt((n_bg*k0)^2 - kx^2 - ky^2)
```

The post-processor records all six DFT fields on the plane. Order powers are formed from the
actual Fourier coefficients of E and H (`S_z = 0.5 Re(Ex Hy* - Ey Hx*)`), including the half-cell
spatial phase between the E plane and the staggered H samples. Evanescent orders report zero
propagating power. Set `propagation_sign = -1` for a reflection monitor; the default `+1` is for
transmission.

`compute_spectrum` returns complex orders and, when `incident_power` is supplied, efficiencies.
`normalize_spectrum(full, reference)` divides every order by the dominant component of the
reference (0,0) order, preserving relative phase. `compute_mask_spectrum` performs a full run and
an unpatterned reference-stack run automatically when `use_reference = true`; it returns raw
orders, normalized orders, incident power, and total transmitted efficiency.

## Accuracy guidance

- Resolve the shortest wavelength in the highest-index material with at least 20 cells per
  wavelength; 30–40 is preferable for phase-sensitive mask spectra.
- Make the CPML thick enough in wavelengths. A ten-cell PML can reflect several percent at
  optical wavelengths; 20–40 layers at λ/20 spacing is a safer default.
- Place the initial pulse in a uniform, lossless layer and keep the pulse envelope away from the
  mask and the PML. The solver validates a ±4σ_z slab around the source and rejects a source that
  crosses a material boundary or enters CPML.
- Material interfaces are rasterized at staggered field locations. A planar interface between two
  z samples has an effective half-cell offset; align critical interfaces to the desired half-cell
  plane or refine the grid. Subpixel interface smoothing is not implemented.
- Use a reference run for efficiencies. The high-level mask API does this when
  `use_reference = true`; the reference removes all `solids`, so put the substrate and other
  unpatterned layers in `mask.layers`.
- The material model is exact at the reference wavelength. Use a narrow pulse around that
  wavelength for quantitative spectra; for broadband work, run separate simulations per
  wavelength or fit multiple poles.

## Validation

`tests/` contains:

- 3D extrusion/solid slicing, bounds, and voxelization.
- Complex-index conversion, Drude round-trip, and active-material rejection.
- Vacuum plane-wave amplitude, impedance, and energy.
- Normal-incidence Fresnel reflection against `(1-n)/(1+n)`.
- Oblique s- and p-polarized Bloch-periodic plane waves, including transverse `E · k`.
- Lossless dielectric-grating energy conservation (`T + R ≈ 1`) with actual-field order powers.
- The high-level mask-spectrum API against the low-level grating result.
- A Drude absorber run with bounded, finite output.

The test suite passes in Release in about 4–6 seconds on the development machine; Debug is slower
because the update kernels are not optimized.

## Performance and memory

The solver stores six field arrays, up to twelve CPML auxiliary arrays (only for CPML axes), and
the diagonal material coefficients for the three electric and three magnetic components. For a
complex (Bloch) run, field and CPML arrays use `std::complex<double>`.

Update kernels are dispatched through `parallel::for_each_chunk` from the sibling `parallel`
work-stealing scheduler. The scheduler splits the z×y rows into cache-sized chunks; the calling
thread participates, and an idle team is reused between updates from the same thread.
`run_options::workers = 0` selects up to four workers, `1` runs inline, and explicit values
are capped at four. This limits thread count without pinning CPUs; macOS still chooses performance
or efficiency cores according to system load. `run_options::chunk_rows`
sets the scheduler grain in rows (default 128; `0` uses the framework default). Monitor accumulation
and material setup remain single-threaded because they are small compared with the field updates.
No OpenMP runtime or build option is used.

H and E updates are fused into two scheduler phases per time step: all three H components share one
parallel region, then all three E components share one. CPML auxiliary variables are read and
written only inside the PML layers; outside, the derivative is the plain Yee stencil.

Two costs limit four-worker scaling: shared cache/memory traffic, and boundary/CPML/material
branches that obstruct SIMD in the general stencil. Earlier tracing already showed about 95.6%
worker utilization; increasing the grain mainly reduced a small amount of waiting. Compiler
remarks confirmed that the original six component loops were not vectorized.

Material coefficients now use lossless row packing: an x-uniform row stores one double, a patterned
row retains every value, and mixed layouts are decoded once per row. The vacuum benchmark's nine
coefficient arrays shrink from 248.8 MB to 2.07 MB. Wrapped cells and x CPML edges retain the general
kernel; y/z CPML state and row coefficients are selected before entering the contiguous bulk loop,
allowing SIMD for real fields. Drude electric updates retain the general formula. The public API
and floating-point accuracy settings are unchanged.

Apple M4, AppleClang 21, Release, 120×120×240, 300 steps, 128-row grain, tracing off.
Each configuration ran three times with old/new versions alternating sequentially; median loop times:

| Workload | Workers | Before | After | Speed relative to old version |
| --- | ---: | ---: | ---: | ---: |
| Vacuum, real fields | 1 | 15.006 s | 6.128 s | 2.45× |
| Vacuum, real fields | 4 | 4.909 s | 2.625 s | 1.87× |
| Lossy magnetic dielectric pattern, real fields | 4 | 4.863 s | 2.660 s | 1.83× |
| Same pattern, oblique complex fields | 4 | 8.010 s | 8.070 s | 0.99× (no stable gain) |

These reruns were contended: a separate FDTD test/benchmark workload and desktop applications were
active on the same machine. They include extra scheduling and bandwidth competition and do not
represent idle-machine peak throughput. Four-worker vacuum loop ranges were 4.631–5.600 s before
and 2.204–3.078 s after; median whole-solve times were 5.037 → 2.791 s. Earlier lighter-load runs
reproduced the original ~7.9 s single-worker and ~3.0 s four-worker results, matching the reported
~2.6× thread speedup. No stable improvement was measured for oblique complex fields in this round.

Compute thread speedup against the single-worker time of the same version. Faster single-core
SIMD can reduce that ratio while still substantially reducing elapsed time. The remaining scaling
limit is consistent with shared data bandwidth; hardware counters have not separately quantified
DRAM and cache effects. Numerical regressions cover real/complex fields, layered/patterned media,
conductivity, magnetism, all-axis/mixed CPML, Drude and narrow grids, including pre-optimization
numerical signatures and identical serial/parallel outputs.

## Limitations

- Spectrum extraction is implemented for z-normal monitors only; field-only monitors may use x
  or y normals.
- One Drude pole, no multi-pole Lorentz/Drude fits, and no anisotropic or magnetic loss.
- No MPI domain decomposition; the shared-memory `parallel` scheduler only.
- No subpixel interface smoothing and no conformal FDTD.
- The high-level reference run removes all solids; patterned substrates should be described as
  planar layers plus separate 3D solids.

## Related documents

- [Chinese documentation](README.zh-CN.md)
- [`litho_geometry` geometry library](../litho_geometry/README.md)

## API checks and tracing

`validate(problem)` checks configuration before voxelization and field allocation.
Source uniformity and material-dependent stability still require `solve`.
When `run.steps` is set, `duration` is ignored, including validation.
`monitor_at(grid, axis, position)` constructs an x/y/z monitor at a physical SI
coordinate and rejects off-grid positions (tolerance: 1e-6 cells).
`component_values` returns an empty span for invalid indices or truncated buffers.
`run_result::workers` reports the selected worker count, including serial automatic
execution on small grids.

The build reuses `litho::trace` or adds the sibling `../litho_tracing` collector
without building its tools, examples, or tests. Override its source path with
`LITHO_FDTD_TRACING_SOURCE_DIR`, or call `find_package(litho_trace CONFIG REQUIRED)`
in the parent project first.

```cpp
litho::Options options;
options.directory = "runs/fdtd-001"; // must be a fresh directory
options.build_id = "fdtd-release";
options.source_revision = "your-git-revision";
auto opened = litho::Session::open(options);
if (!opened) { /* opened.error().describe() */ }
auto& trace = opened.value();
p.run.trace.session = &trace; // works with problem and mask_problem
p.run.trace.context = trace.stream("mask-job", {}, LT_SOURCE);
p.run.trace.step_interval = 1;   // sample every step; 0 disables kernel spans
p.run.trace.detail_steps = 1;    // also emit per-chunk worker spans for step 0
{
    auto result = compute_mask_spectrum(p);
    if (!result) { /* message(result.error()) */ }
}
trace.close();
auto stats = trace.stats(); // inspect dropped_events and error
```

The caller owns the Session, which must outlive the solve. A null Session disables
instrumentation. Collection is best effort: synchronous tracing exceptions do not
change numerical results. Session statistics expose collector drops and I/O errors.

Events carry source locations and explicit parent contexts. Category `fdtd` covers
configuration, preparation, completion, errors, and stage spans: mask spectrum,
reference/patterned runs, solve, voxelization, coefficients, initialization, time
loop, and postprocessing. Category `fdtd.step` records each block and its progress.
Category `fdtd.kernel` records one span per update component
(`update_h`, `update_e`) and monitor accumulation on **the first step of each
block**. Every such span carries a `summary` instant with:

| Field | Meaning |
| --- | --- |
| `wall_ns` | phase wall time on the calling thread |
| `main_work_ns` | time the calling thread spent inside its own chunks |
| `wait_ns` | `wall_ns - main_work_ns` |
| `work_ns` | sum of all chunk durations; only when `detail_steps` applies |
| `chunks` | executed chunk count |
| `utilization` | `work_ns / (wall_ns * workers)`; only for detailed steps |
| `main_fraction` | `main_work_ns / wall_ns` |

`detail_steps` (default 1) additionally records child spans in category
`fdtd.chunk` for the first N sampled steps, so the timeline shows each worker's
active intervals and the barrier tail. A high `utilization` with a poor speedup
means per-worker throughput dropped (memory/cache contention); a low
`utilization` means workers are idle or imbalanced. Filter categories with
`Options::categories`.

`elapsed_seconds` measures the time loop including any enabled instrumentation;
the outer solve span and benchmark wall time also cover setup and postprocessing.

```sh
./build/release/benchmarks/litho_fdtd_benchmark 4 runs/fdtd-w4
./build/release/benchmarks/litho_fdtd_benchmark 4 monitors runs/fdtd-benchmark-monitors
../litho_tracing/build/litho_analyze runs/fdtd-w4 summary
../litho_tracing/build/litho_analyze runs/fdtd-w4 serve web/dist 8090
../litho_tracing/build/litho_analyze runs/fdtd-w4 export runs/fdtd-w4/perfetto.json
./build/release/benchmarks/litho_fdtd_spectrum_benchmark
```

Build the analyzer separately in `litho_tracing`. The monitor workload records four
frequencies on three z planes. The spectrum benchmark extracts 169 orders from a
120×120 plane, repeats five times, and prints a checksum. Separable x/y phases and
factored stagger offsets reduce complex exponentials from `5*nx*ny` to `nx+ny+2`
per order. Monitor accumulation now traverses each channel contiguously and avoids
per-step phase-array allocations.

Local AppleClang 21 Release comparison (identical workloads run sequentially,
tracing disabled):

| Workload | Before | After |
| --- | ---: | ---: |
| 120×120 spectrum, 169 orders, average of five extractions | 81.023 ms | 13.703 ms |
| 120×120×240, 300 steps, 4 workers, three four-frequency monitors, median of three loops | 4.553 s | 4.254 s |

These correspond to about 5.9× faster spectrum extraction and 6.6% less time in
the monitored loop, rather than a universal field-kernel speedup. The field-loop
entry predates the H/E fusion and chunk experiments below.

The following historical H/E fusion measurements predate the row packing and SIMD optimization above.
The fused four-worker loop was measured against the pre-fusion trace on the same
120×120×240, 300-step vacuum grid with no monitors (`step_interval = 1`, `detail_steps = 1`):

| Configuration (4 workers) | Loop time | Step-0 phase work | Worker utilization | Worker idle |
| --- | ---: | ---: | ---: | ---: |
| Before: 6 phases, full-grid CPML, 64-row grain | 4.63 s | 57.1 ms | 0.805 | 19.5% |
| Fused H/E, boundary-only CPML, 64-row grain | 3.27 s | 40.9 ms | 0.924 | 7.6% |
| Same, 128-row grain (default) | 3.28 s | 44.0 ms | 0.956 | 4.4% |
| Same, 256-row grain | 3.29 s | 43.5 ms | 0.969 | 3.1% |
| Same, 512-row grain | 3.36 s | 42.5 ms | 0.976 | 2.4% |

The serial loop dropped from 12.11 s before fusion to about 8.0–8.3 s with one worker. The four-worker
step-0 speedup is `4 * 0.956 / 1.56 = 2.45`, matching the measured ~2.5×: fusion and boundary-only
CPML reduced absolute work, but concurrent workers still inflate total work by 1.56× because of
memory/cache bandwidth and per-core throughput. Larger grains reduce worker idle (19.5% → 4.4% →
2.4%) but do not improve wall time once the tail becomes coarse; 128 rows is the default balance.
The main thread itself waits only a few percent; the earlier waiting was on workers at phase
boundaries. Traced runs wrote 3,622–5,430 events with zero drops. These numbers are
hardware-specific.

For future batch workloads, reference-spectrum reuse and prepared grid/material
coefficients are useful candidates. The current high-level API still voxelizes
and runs the reference for each request; no cache-lifetime API has been added.
