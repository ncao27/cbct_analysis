# Case B registration: delayed phase onto portal venous phase

The existing `demons_iodine_registration.m` ran successfully with its default
settings on a matched real 2D phase pair and a local 3D slab. It improved several
structural-alignment measures and passed a known-motion control. No functional
refactoring was needed, and the registration function remains byte-for-byte
unchanged. This is bounded validation, not proof of complete anatomical accuracy.

## Reference phase and tested data

**PV is fixed; DEL moves.** Portal venous is an established reference choice for
multiphase abdominal CT, rather than a universal requirement for every task.
A [published multiphase liver CT study](https://link.springer.com/article/10.1007/s00330-008-0996-8)
registered its other phases to a portal venous reference.
[Siemens' liver-analysis workflow](https://academy.siemens-healthineers.com/_/en-us/syngo-via-ct-liver-analysis/)
also uses a selected liver series, often venous or the most contrasted series,
as the reference. The intended analysis determines the target coordinate system.
The first argument of our function is the fixed reference, even though the
original pre/post-iodine example calls it `preImage`.

The available Case B inputs are **both post-contrast**, not a true pre-injection
and post-injection pair:

- Fixed: `/Users/ncao/Desktop/cbct_analysis/lesion_insertion_project/cases without lesion/case b/ABD PEL AX 3 PV_40cm_512x512x171.raw` (171 slices).
- Moving: `/Users/ncao/Desktop/cbct_analysis/lesion_insertion_project/cases without lesion/case b/ABD PEL AX 3 DEL_40cm_512x512x142.raw` (142 slices).

The native 2D test matches **PV slice 91 to DEL slice 71**. It searches all PV
slices using bone-focused cross-correlation with a bounded translation search.
PV 91 was the strongest candidate; adjacent candidates scored lower. The images
were also inspected visually. This is image-content correspondence, not a
verified DICOM slice-position match. The initial correlation peak had no
in-plane shift at the coarse search resolution, so no separate rigid step was
added for this already closely aligned pair.

The 3D test uses **PV 83-99 and DEL 63-79**, resampled in-plane to 256 x 256 while
retaining all 17 slice planes. It tests 3D code on a local overlap only. All raw
inputs remain untouched; no reconstruction geometry or detector parameters changed.

## Quantitative results

For the 2D test, Dice is an overlap score (1 means perfect overlap). Masks are
simple stored-intensity proxies, not expert organ or bone segmentations.

| Measure | Before | After default demons |
| --- | ---: | ---: |
| Body-mask Dice | 0.990984 | 0.998571 |
| Bright-structure Dice | 0.809401 | 0.849695 |
| Symmetric source-edge distance, pixels | 1.721577 | 1.406336 |
| Edges re-detected after interpolation, pixels | 1.721577 | 1.938095 |
| Phase intensity disagreement, RMSE | 171.358691 | 147.058834 |
| Phase intensity disagreement, MAE | 65.840872 | 49.429722 |
| Optimized-feature RMSE | 0.031046 | 0.020628 |

The phase-intensity measurements use 84,215 pixels on the same
fixed body-interior/common-support mask before and after. They are **stored pixel
units, not HU, and not ground-truth registration errors**: PV and DEL genuinely
have different enhancement, bowel contents/motion, and potentially through-plane
anatomy. Optimized-feature RMSE is not an independent accuracy measure.

The edge measurements are deliberately both reported. Propagating the original
moving edges with the deformation improved their geometric agreement with fixed
edges. Re-detecting edges in the interpolated output worsened the score. Edge
detection can change with interpolation and contrast, but this does not prove
that every local structure is correctly aligned. Visible residual mismatches
remain, especially internally; neither score replaces landmark validation.

A bounded alternative used the existing `RegistrationPair` option with clipped
7400-8800 feature copies and `FieldSmoothing=2`. It yielded propagated-edge distance
1.484320, re-detected-edge distance
1.644261, and phase RMSE
153.978169. It did not improve all measures. The default remains
the primary result; the alternative arrays and metrics are retained for review.

## Deformation checks and controls

- 2D minimum finite-difference Jacobian: **0.335112**; nonpositive
  fraction **0.000000** (no detected folds).
- 2D valid interpolation support: **99.416%**;
  maximum displacement **8.493 pixels**.
- 3D minimum Jacobian: **0.323781**; nonpositive fraction
  **0.000000**; valid support **92.241%**.
- 3D phase RMSE: **120.948 -> 85.083**
  on 302,194 common-support voxels with fixed stored values
  above 7700, excluding two planes at either end. This ROI differs from the 2D
  ROI, so the two RMSE values should not be compared directly.
- Known-motion control: the real PV image was synthetically translated by
  (+4 columns, -3 rows). Median backward-field endpoint error was
  **0.033471 pixels**, with 95th percentile
  **0.113703 pixels**. Intensity RMSE fell from
  **206.618 to 1.411**. This control has known truth;
  the actual patient phase pair does not.
- Identity input returned exactly the same values and a zero field. The final
  real-phase output was verified to equal a linear warp of original DEL values
  on valid samples; unsupported samples are NaN. Feature normalization was not
  applied to the returned intensities.

All smoothing/displacements use the tested grid's pixel/voxel units. The 3D grid
is not asserted to be physically isotropic, and positive sampled Jacobians do
not establish clinically valid anatomy or a continuous diffeomorphism.

## Limits and reproducibility

Both volumes contain two bytes per voxel, conflicting with the reported
float32 export setting. They are interpreted as inferred big-endian uint16,
consistent with the previous experiments. HU scaling/offset remains unknown.
The supplied 40 cm in-plane FOV implies 0.78125 mm native pixel spacing; 3 mm was
reported as slice thickness, **not confirmed inter-slice spacing**. DICOM origins,
direction cosines, and slice positions are absent. Consequently this is **not a
validated full-volume registration in patient coordinates**. Do not resize the
171-slice and 142-slice volumes to equal depth or assume equal indices match.
Full-volume use needs geometry-aware alignment and independent anatomical checks.

Parameters remain the existing defaults: iterations [100 50 25], field smoothing
1.5 voxels, local normalization sigma 5, noise floor 0.03. No iodine map was
registered to a conventional CT image. Original source/intensity values were
preserved apart from the interpolation inherently required by registration.
[MATLAB's demons documentation](https://www.mathworks.com/help/images/ref/imregdemons.html)
describes the backward sampling field and its pixel/voxel component convention.

Arrays in this folder are in **display orientation**, already transposed relative
to raw-file MATLAB slices and the earlier reconstruction MAT files. D(:,:,1) is
column displacement and D(:,:,2) is row displacement. For 3D, component 3 is slice
displacement. Do not apply these fields to untransposed reconstruction arrays.

From this folder in MATLAB:

```matlab
validate_case_b_registration('compute');
validate_case_b_registration('comparison');
validate_case_b_registration('differences');
validate_case_b_registration('slab');
validate_case_b_registration('arrays');
validate_case_b_registration('metrics');
```

Requires Image Processing Toolbox. Default cache is beneath MATLAB tempdir;
this run supplied `/private/tmp/cbct_registration_case_b` as the optional second
argument to every call. No-input calls to the original demons function generate
its separate synthetic demo; this validation does not overwrite that demo.

## Files

- [01_registration_comparison.png](01_registration_comparison.png): fixed PV,
  original DEL, registered DEL, matched-window overlays and displacement magnitude.
- [02_phase_subtractions.png](02_phase_subtractions.png): before/after DEL minus
  PV, shown at +/-200 stored units on identical valid support. Bright/dark
  residuals include genuine contrast differences and remaining misalignment.
- [03_slab_and_motion_control.png](03_slab_and_motion_control.png): local 3D center
  plane and the real-image known-shift control.
- [04_registration_arrays.mat](04_registration_arrays.mat): images, displacement
  fields, valid masks, structural masks/edges, slice-match scores, diagnostic
  trial, 3D slab, known-motion control, metadata and metrics.
- [05_metrics.json](05_metrics.json): complete settings, numerical results and
  explicit scope/encoding limitations.
- [validate_case_b_registration.m](validate_case_b_registration.m): reproducible
  computation, checks, and per-file export stages.

All 91 pre-existing tracked/data files in the SHA-256 baseline remain
byte-for-byte unchanged, including `demons_iodine_registration.m`, `li_fbp/`,
C++/CUDA kernels, original scans, and previous results. New validation files are
committed separately in accordance with repository practice.

Source hashes:

- Registration function: `ed375c02e5421b7694a5e90ce861361f83c24f14840f086e93ade6ab64b4cf8d`.
- PV raw volume: `2d183adeddd2e1c4ae43350612038ab1260001012a42242ba7284fe92363b709`.
- DEL raw volume: `135f1aec17ccc5c762a7e393768683a35f7e81f9e57d0411632a7af21dd513e6`.
