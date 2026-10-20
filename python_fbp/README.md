# Python fan-beam reference implementation

This directory ports all four scripts in `../li_fbp/`, which remains unchanged.
`fbp.py` shares the ray integration, cubic interpolation, raw I/O, filtering, and
backprojection code between the linear and curved detector geometries.

## Run

Use Python 3.10 or newer and run these commands from the repository root, inside
your Python virtual environment:

```sh
python -m pip install -r python_fbp/requirements.txt
python -m python_fbp run --detector linear --plot
python -m python_fbp run --detector curved --plot
```

The default input is `li_fbp/without_lesion.raw`. Outputs are separated by detector:

- `python_fbp/results/<detector>/projections.raw`
- `python_fbp/results/<detector>/img.raw`
- `python_fbp/results/<detector>/run.json` (parameters, backend, and array layouts)
- `images/python_fbp/<detector>/reconstruction.png` when `--plot` is supplied

A repeated run overwrites that detector's outputs. Use `--output-dir` to keep
individual experiments. Both the CLI and raw writers reject output paths under
`li_fbp/`, including paths resolving there through symlinks.

To run the two steps independently:

```sh
python -m python_fbp forward --detector linear
python -m python_fbp reconstruct --detector linear --plot
```

`--input PATH` selects a big-endian raw image for `forward`/`run`, or a projection
raw file for `reconstruct`. The latter must match the selected detector geometry
and original sampling; raw files contain no geometry metadata. The default
`reconstruct` input is the projection file in its output directory.
`--chunk-rays 64` bounds forward-projector memory use without changing ray
positions or sampling. CPU execution uses NumPy by default.

For an NVIDIA CUDA machine, install the [CuPy distribution matching its CUDA
installation](https://docs.cupy.dev/en/stable/install.html), then pass
`--backend cupy`. The implementation uses the same formulas with CuPy arrays.
This backend has an optional regression test but has **not been executed on the
current Mac**, which has no supported CUDA GPU.

## Mapping and preserved parameters

| MATLAB script | Python operation |
| --- | --- |
| `c1_forward_projection_v6_linear_detector.m` | `forward_project(image, "linear")` |
| `c1_forward_projection_v6_curved_detector.m` | `forward_project(image, "curved")` |
| `c2_fbp_v6_linear_detector.m` | `reconstruct(projections, "linear")` |
| `c2_fbp_v6_curved_detector.m` | `reconstruct(projections, "curved")` |

The production defaults are unchanged: 512 × 512 images, FOV 400 mm,
SID 750 mm, SDD 1200 mm, 2480 detector elements, detector spacing 0.35 mm,
and 2000 views on the original grid from zero through `2*pi - 2*pi/2000`.
Forward ray sampling retains `ceil((sid+sdd)/(dx/2))` midpoint samples per ray.
The `Geometry` API supports explicit overrides for synthetic regression tests;
the CLI uses the original defaults.

```python
from python_fbp import read_image, forward_project, reconstruct, write_image

image = read_image("li_fbp/without_lesion.raw")
projections = forward_project(image, "curved")
reconstruction = reconstruct(projections, "curved")
write_image("python_fbp/results/curved/img.raw", reconstruction)
```

## Numerical and layout details

- Image array row/column indices match MATLAB's in-memory matrix, including the
  transpose used by the reference display. Input files are big-endian float32
  read into double arrays, matching MATLAB's `fread(..., 'single')` default.
- Projection arrays have shape `(views, detectors)`. Their little-endian float32
  row order matches the reference forward projector's per-view writes.
- Reconstruction files are little-endian float32 in MATLAB column order. The
  display uses `image.T` with the reference intensity range [0.015, 0.035].
- Cubic interpolation implements cubic convolution with quadratic ghost samples
  at the image boundary and zero extrapolation. It was checked against MATLAB
  `interp2`; it is not substituted with a cubic B-spline.
- Forward geometry follows the reference's single-precision arithmetic, while
  interpolated samples retain the input image's double precision. Reconstruction
  matrices and detector axes use single precision, as in the reference.
- The curved reconstructor's implicit workspace dependency on lowercase `sdd`
  is made explicit through `Geometry.sdd`, preserving the forward script's
  single-precision angular spacing.
- The backprojectors intentionally retain the original unnormalized interpolation
  weight (`position - detector_coordinate`, without dividing by spacing) and the
  first-bin fallback for invalid detector indices. These behaviors are part of
  this port; no algorithmic correction or attenuation clipping is introduced.
- Ramp filtering uses the symmetric kernel's equivalent full convolution and
  central slice, avoiding the reference's repeated ramp matrix allocation.
  Floating-point reduction order and CPU/GPU math libraries can differ.

See [MATLAB interp2](https://www.mathworks.com/help/matlab/ref/interp2.html) for the
reference interpolation method and [NumPy reshape order](https://numpy.org/doc/stable/reference/generated/numpy.reshape.html)
for the explicit column-order convention.

## Verification

Run the regression suite from the repository root:

```sh
python -m unittest discover -s python_fbp/tests -v
```

The included fixture was generated by MATLAB R2023b Update 11 on Apple silicon,
using temporary in-memory copies of the reference scripts with GPU/display
plumbing replaced by CPU equivalents. `li_fbp` was never edited.
The complete synthetic pipeline uses a 12 × 12 image, 16 detector elements,
9 views, and 35 mm pitch solely for testing; production defaults are unchanged.
The fixture also includes four projections of the included 512 × 512 input at
original zero-based view indices `[0, 250, 500, 750]`, with every production
geometry and sampling value retained.

| Comparison against MATLAB CPU | Linear relative L2 error | Curved relative L2 error |
| --- | ---: | ---: |
| Synthetic forward projection | 9.02e-8 | 2.59e-7 |
| Synthetic ramp filtering | 2.75e-7 | 3.19e-7 |
| Synthetic backprojection | 5.56e-7 | 6.22e-8 |
| Four production-size forward views | 1.54e-7 | 9.12e-7 |

The suite also checks cubic edges/extrapolation, file byte order and matrix
orientation, chunk independence, streaming output, legacy boundary behavior,
CLI integration, and protected output paths. The optional CUDA test skips when
CuPy or a supported device is unavailable. Source/input SHA-256 values in the
fixture bind the comparisons to the preserved MATLAB files and image.

These checks establish agreement with CPU evaluations of the source formulas.
They do **not** establish bitwise GPU equivalence, full 2000-view reconstruction
agreement, or clinical correctness. A [synthetic comparison figure](../images/python_fbp_validation/matlab_comparison.png)
shows both reconstructions and their differences; it is a test visualization,
not a patient reconstruction.

To regenerate fixtures, run `python_fbp/tests/generate_matlab_fixtures.m` in MATLAB,
then package the raw files from the output directory printed by that script:

```sh
python python_fbp/tests/pack_matlab_fixtures.py /path/to/cbct_python_validation
```

The generator locates the repository from its own path. Optional environment
variables `CBCT_REPOSITORY` and `CBCT_FIXTURE_OUTPUT` override the repository and
fixture output paths. Keep the MATLAB version/provenance metadata in the packer
consistent with the environment used to regenerate the fixture.
