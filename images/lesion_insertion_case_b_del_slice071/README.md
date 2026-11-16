# Paired image-domain and projection-domain lesion insertion

This experiment uses the same healthy Case B DEL slice 71/142 as the preceding
round-trip test and the same real segmented DEL lesion in both methods. The owner
selected **replacement**, rather than adding absolute lesion intensities to
healthy values. The original data and `li_fbp/` scripts were preserved.

## What was compared

Let H be the healthy image, M the embedded nonzero lesion mask, L the absolute
lesion values, A the curved-detector forward projector, and R the reference
filtered backprojection. Define delta = M .* (L - H) elementwise, zero outside the patch.

1. **Direct image-domain insertion:** I = H + delta. Every masked pixel is
   replaced exactly by its corresponding lesion value; all other pixels are unchanged.
2. **Projection-domain insertion:** P = R(A(H) + A(delta)). The negative/positive
   replacement difference is projected, avoiding double-counting healthy tissue.
3. **Matched reconstruction control:** C = R(A(H + delta)), computed independently.
   Comparing P with C measures whether insertion before or after projection gives
   the same reconstructed result. Comparing P with I also includes reconstruction
   error, since I has not passed through projection and reconstruction.
4. **Incremental lesion error:** [P - R(A(H))] - delta. Subtracting the healthy
   reconstruction isolates the reconstruction's effect on the inserted change.

No noise was added, and no registration, fitted rescaling, intensity offset,
smoothing, mask interpolation, or image resizing was applied.

## Numerical results

All intensity errors are in **stored pixel units, not HU**. Signed differences are
projection-domain minus image-domain. The matched comparison uses C as the image-domain output.

| Comparison | Region | Pixels | MAE | RMSE | Bias | Max absolute |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Projection insertion vs direct image replacement | all_pixels | 262,144 | 78.598224 | 289.528907 | -19.377942 | 5221.080078 |
| Projection insertion vs direct image replacement | lesion | 246 | 51.827990 | 65.082317 | -5.609923 | 176.747070 |
| Projection insertion vs direct image replacement | interior | 244,596 | 53.437697 | 67.297616 | -0.370856 | 325.438477 |
| Projection insertion vs reconstructed image replacement | all_pixels | 262,144 | 0.002708 | 0.003573 | +0.000001 | 0.051758 |
| Projection insertion vs reconstructed image replacement | lesion | 246 | 0.012455 | 0.015598 | -0.000931 | 0.051758 |
| Projection insertion vs reconstructed image replacement | interior | 244,596 | 0.002734 | 0.003612 | +0.000000 | 0.051758 |
| Inserted change vs intended replacement change | all_pixels | 262,144 | 0.202288 | 0.342349 | -0.000934 | 30.853516 |
| Inserted change vs intended replacement change | lesion | 246 | 3.306893 | 5.482363 | +0.116511 | 30.853516 |
| Inserted change vs intended replacement change | interior | 244,596 | 0.196229 | 0.339343 | -0.001087 | 30.853516 |

Projection linearity residual, A(H) + A(delta) - A(H + delta):
relative L2 = 8.87467468e-09 and maximum absolute
residual = 0.25 in projection units
(stored pixel units times mm).
The matched reconstructed images have relative L2 difference
4.77524877e-07.

These are numerical consistency checks for this one image, lesion, and CPU run.
They are not validation of physical lesion contrast, realistic acquisition noise,
or diagnostic performance. Agreement after matched reconstruction does not mean
that the reconstructed lesion exactly equals its source patch.

## Input, mask, and placement

Healthy source: `/Users/ncao/Desktop/cbct_analysis/lesion_insertion_project/cases without lesion/case b/ABD PEL AX 3 DEL_40cm_512x512x142.raw`.
Healthy image and forward projections are reused exactly from
`../lesion_free_roundtrip_case_b_del_slice071/03_roundtrip_arrays.mat`.

Lesion source: `lesion_insertion_project/lesion modesl/case a/S13_I137_DEL_19x18.raw`.
SHA-256: `bbe8e1200bde37827d6fcbfff7372fe984914ad2ed20910bcf3f24be5e62bac6`.

The 684-byte patch is interpreted as big-endian uint16 with reference matrix size
19 x 18. Its 246 nonzero pixels define the mask; zero pixels are background and
never overwrite the healthy image. Lesion values range from
8185 to 8315.
The patch occupies **display rows 251-268 and columns 166-184**, one-based, with no
resampling or rotation. Placement is illustrative, not anatomically targeted.
Lesion pixel spacing is unknown: a 1:1 pixel mapping is used on the healthy
0.78125 mm grid, so physical lesion-size fidelity is not established.

Both raw encodings remain inferred. The healthy file contains two bytes per
voxel despite the reported float32 metadata; its earlier analysis used big-endian
uint16. No HU scaling or attenuation calibration is available. Replacement assumes
the lesion and healthy stored values have compatible scaling, which is not
independently established from export metadata.

The interior region is more than five pixels from both image boundaries and
stored zero padding. The JSON also reports the area outside the lesion and a
three-pixel square ring around it. These are geometric masks, not anatomical segmentations.

## Geometry, execution, and checks

Used the preserved `li_fbp/c1_forward_projection_v6_curved_detector.m` and
`li_fbp/c2_fbp_v6_curved_detector.m` formulas on CPU. Temporary copies replace GPU
allocation/wrappers, redirect paths, and suppress plotting/logging; source files
are not edited. The FBP script's inherited lowercase `sdd` is supplied as the
same single-precision 1200 used by its forward script. Projection storage,
addition, and FBP retain the reference single precision. This checks CPU
execution; GPU equivalence was not tested.

Geometry remains 512 x 512, FOV 400 mm, SID 750 mm, SDD 1200 mm, 2480 curved
detector elements with 0.35 mm pitch, and 2000 views, with the original ray
sampling, filtering and interpolation. Supplied slice thickness 3 mm is metadata
only in this 2D test. Source hashes and the preceding validation provenance are
in `05_metrics.json`.

Timings: delta forward projection 101.16 seconds;
combined-image forward projection 109.45 seconds;
three filtering/backprojection runs 19.65 seconds.

Checks passed: exact replacement; exact preservation outside the mask; exact
float32 conversion of the integer input/delta images; finite 512 x 512 outputs
and 2000 x 2480 projections; projection linearity relative L2 below 1e-6; matched
reconstruction relative L2 below 1e-5 and maximum absolute difference below 0.25;
and reconstruction additivity relative L2 below 1e-5. These thresholds were set
before execution for numerical consistency, not clinical acceptance.
Exported arrays were reloaded, and all 262,144 CSV rows were checked for subtraction
consistency to within 1e-7. Figures were visually inspected.

## Results and reproduction

- [01_insertion_comparison.png](01_insertion_comparison.png): healthy image,
  both insertion outputs, and zoomed comparison including the matched reconstruction.
  Shared display window 7000-9000; orange rectangles mark the patch bounds.
- [02_method_differences.png](02_method_differences.png): signed subtraction maps;
  direct-comparison panels saturate at +/-200 for visibility, while the matched
  reconstruction difference uses its full symmetric range. Gray is zero.
- [03_insertion_arrays.mat](03_insertion_arrays.mat): all inputs, mask, replacement
  delta, projections, reconstructed images, errors, metadata, and metrics.
- [04_pixel_differences.csv](04_pixel_differences.csv): coordinates, healthy values,
  mask, delta, both method outputs, matched control and three differences for every pixel.
- [05_metrics.json](05_metrics.json): complete statistics and definitions.
- [validate_lesion_insertion.m](validate_lesion_insertion.m): executable experiment
  and artifact exporter.

MAT arrays retain reference indexing; figures and CSV use the transpose, so
CSV display row r, column c corresponds to MAT element (c,r). CSV rows are ordered
by display column, then row. Continuous numerical differences are not quantized.

From this folder in MATLAB:

```matlab
validate_lesion_insertion('compute');
% Then export one artifact at a time:
validate_lesion_insertion('comparison');
validate_lesion_insertion('difference');
validate_lesion_insertion('arrays');
validate_lesion_insertion('pixels');
validate_lesion_insertion('metrics');
```

The default cache is a subdirectory of MATLAB's tempdir. Each call accepts the
same optional second argument to use another cache directory. This run used
`/private/tmp/cbct_lesion_insertion_validation`. Requires the prior healthy
validation arrays, original reference scripts, raw lesion file, and Image
Processing Toolbox for the geometric interior mask.

Preservation audit: all 81 pre-existing tracked and dataset files in the pre-run
SHA-256 manifest remain byte-for-byte unchanged, including the golden MATLAB
reference, kernels, original raw inputs, and preceding verification results.
