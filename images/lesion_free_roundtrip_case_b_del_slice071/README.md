# Lesion-free forward projection / filtered backprojection verification

The curved-detector MATLAB reference completed a full round trip on Case B DEL,
slice 71 of 142. The reconstruction preserves the visible image structure, with
fine-grained differences and substantial boundary errors. This demonstrates
execution and measures fidelity for this one input; it does not establish an
accuracy threshold, validate the linear-detector code, or validate GPU execution.

## Results

Open [the before/after comparison and error maps](02_before_after_and_errors.png).
Both image panels use the same display window, 7000–9000 stored pixel units.
Difference panels show the complete error range. Display windowing does not
change the matrices or metrics.

All errors below are in **stored pixel units, not HU**.

| Region | Pixels | MAE | RMSE | Mean signed error | Maximum absolute error |
| --- | ---: | ---: | ---: | ---: | ---: |
| All pixels (primary result) | 262,144 | 78.5983 | 289.5290 | -19.3770 | 5220.7271 |
| Nonzero ground truth | 255,616 | 76.3833 | 289.8749 | -22.5563 | 5220.7271 |
| Interior diagnostic | 244,596 | 53.4378 | 67.2978 | -0.3698 | 325.7778 |

Signed error is reconstruction minus ground truth. MAE is mean absolute error;
RMSE is the square root of mean squared error. Full-image median absolute error
is 46.3438, and the 95th percentile is 147.7674.

The maximum error occurs at display row 512, column 256: ground truth 7271,
reconstruction 2050.27294921875. All of the largest 2622 absolute errors (the
largest 1%, rounded up) are within five pixels of an image edge or stored zero
padding. The interior diagnostic excludes both these boundary regions and zero
pixels; it is not a body or tissue segmentation. It supplements the full-image
result rather than replacing it. The reconstruction range is -444.4174 to
9808.0400; values were not clipped for measurement.

## Input and encoding limitation

Source, relative to the repository root:
`lesion_insertion_project/cases without lesion/case b/ABD PEL AX 3 DEL_40cm_512x512x142.raw`.
The lesion-free designation comes from the source folder. Slice 71 is the chosen
middle slice, using one-based indexing, without lesion insertion or registration.

Owner-provided metadata: 512 × 512 × 142, reconstruction FOV 40 cm, 3 mm slice
thickness, delayed contrast phase, float32 big-endian. However, this file is
74,448,896 bytes: exactly two bytes per voxel. Float32 at the supplied dimensions
would require 148,897,792 bytes. The run therefore uses **inferred big-endian
uint16** data, which produces a coherent image. Big-endian signed and unsigned
16-bit interpretations agree on this slice (range 0–9818). This does not confirm
the original export encoding. HU scaling/intercept and attenuation calibration
are unknown; none were invented or applied. A confirmed export could require
repeating this test.

The slice starts at byte offset `70 * 512 * 512 * 2`. MATLAB read it using
`fread(fid,[512 512],'uint16=>double')` with big-endian byte order. The extracted
slice was written to a temporary big-endian float32 file for the reference
loader; equality after conversion was checked exactly. This tests the numerical
projection/reconstruction of stored values, not calibrated physical attenuation.

## Reference implementation and execution

Used `li_fbp/c1_forward_projection_v6_curved_detector.m` and
`li_fbp/c2_fbp_v6_curved_detector.m` at the source revision recorded in
[05_metrics.json](05_metrics.json). Source and input SHA-256 hashes are recorded
there. The originals were not edited.

MATLAB R2023b Update 11 ran on CPU because no MATLAB-supported GPU was available.
Temporary in-memory copies removed `gpuDevice(1)`, replaced `gpuArray`/`gather`
wrappers with their CPU values, replaced GPU array allocations with CPU
allocations, redirected file paths, and suppressed reference plotting while
reducing progress logging. The reference projection, filtering, interpolation,
and backprojection formulas were retained, including their existing boundary
behavior. GPU numerical equivalence was not tested.

Retained parameters: 512 × 512 pixels; FOV 400 mm; pixel spacing 0.78125 mm;
source-to-isocenter distance 750 mm; source-to-detector distance 1200 mm;
2480 curved-detector elements at 0.35 mm pitch; 2000 views. The original sampling
rule was retained. The supplied 3 mm slice thickness is metadata only for this
2D calculation. Forward projection took 94.62 seconds; filtering/backprojection
took 6.23 seconds. No noise was deliberately added.

No fitted gain/offset, registration, independent image normalization, or
post-reconstruction denoising was used. Both images are transposed identically
for display, following the reference plotting convention.

## Files and checks

- [01_ground_truth.png](01_ground_truth.png): input image and source description.
- [02_before_after_and_errors.png](02_before_after_and_errors.png): matched-window
  before/after images, signed differences, and absolute differences.
- [03_roundtrip_arrays.mat](03_roundtrip_arrays.mat): ground truth, reconstruction,
  signed/absolute/squared errors, 2000 × 2480 projections and filtered projections,
  diagnostic masks, geometry, execution metadata, and metrics. MATLAB v7.3 format.
- [04_pixel_errors.csv](04_pixel_errors.csv): all 262,144 pixels, with one-based
  display row/column, ground truth, reconstruction, signed error, and absolute error.
- [05_metrics.json](05_metrics.json): complete numerical summaries, percentile
  definitions, mask definitions, geometry, encoding caveat, and provenance.

MAT arrays retain reference indexing. CSV display row `r`, column `c` corresponds
to MAT array element `(c,r)` and to the figures. CSV rows are ordered by display
column, then display row; use the coordinate columns rather than assuming raster
order. Error values remain continuous floating-point values, without quantization.

Checks passed: expected dimensions; finite image/projection values; exact input
conversion; signed/absolute/squared error identities; reloaded MAT consistency;
CSV row count and agreement with computed errors to within 1e-7; and JSON RMSE
agreement to within 1e-9. Figures were visually inspected. All 74 original tracked
source/research files and untracked dataset files in the pre-run manifest remain
byte-for-byte unchanged. No lesion data, kernels, or reference code were modified.
