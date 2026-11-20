# Regenerated reconstruction and lesion-insertion images

The excess speckle is substantially reduced using normalized detector interpolation and explicit correction for an assumed stored-value baseline. No image smoothing was applied. These are corrected experimental previews, not calibrated HU/attenuation reconstructions. The correction uses the original image's known support; its accuracy is not an independent assessment from projections alone.

- [Healthy reconstruction: original, regenerated, and subtraction](01_reconstruction.png)
- [Image-domain and projection-domain lesion insertion](02_lesion_insertion.png)
- [Pixel differences and matched-method comparison](03_pixel_differences.png)
- [Numerical images, masks, and pixel error maps](04_arrays.mat)
- [Metrics and method metadata](05_metrics.json)
- [Reproduction script](regenerate_previews.m)

The healthy input remains Case B DEL slice 71. Both insertion methods use the original Case A `S13_I137_DEL_19x18.raw` lesion, the same 246-pixel mask, replacement values, and placement. The original lesion patch is also included in `04_arrays.mat`; its unchanged display preview is [here](../display_window_minus200_to300/03_inserted_lesion.png).

## What changed

The separate MATLAB workflow reads the curved-detector reference scripts and adapts copies in memory for CPU execution. It normalizes the detector interpolation fraction by `d_gamma` and zeros contributions from invalid detector indices. Reference files, geometry, detector parameters, and ray sampling are unchanged.

For the large stored-value baseline, let h be the healthy image, delta the original lesion replacement difference, A the original forward projector, and Rc the corrected backprojector. The workflow uses:

```text
p = 8192 * (h ~= 0)
healthy reconstruction = Rc(A(h) - A(p)) + p
projection-domain insertion = Rc(A(h) - A(p) + A(delta)) + p
matched image-domain control = Rc(A(h + delta) - A(p)) + p
```

The same p is used for both methods, and the code verifies that insertion does not change its support. A(p) is computed afresh; A(h), A(delta), and A(h+delta) come from the preceding saved validation. The baseline p is restored exactly after reconstructing the residual. This prevents the reconstruction error of that large known baseline from dominating the images. It does not establish that 8192 is the true export offset or that the resulting values represent physical attenuation.

The display remains stored value minus 8192 with limits [-200, 300]. All error metrics use unshifted stored units. Subtraction panels identify their signs and scales; finer error scales deliberately expose the remaining errors. No values are zeroed or clipped in the saved numerical arrays for presentation.

## Measured results

RMSE values are in stored pixel units, not HU.

| Comparison | Original RMSE | Regenerated RMSE |
| --- | ---: | ---: |
| Healthy reconstruction, interior | 67.2978 | 6.3099 |
| Projection insertion vs direct image insertion, interior | 67.2976 | 6.3115 |
| Projection insertion vs direct image insertion, lesion mask | 65.0823 | 7.0840 |
| Healthy reconstruction, all pixels | 289.5290 | 36.3980 |
| Projection insertion vs direct image insertion, all pixels | 289.5289 | 36.3982 |

The interior mask is the original 244,596-pixel mask, excluding pixels within five pixels of the image edge or zero padding. It includes background and is not tissue-only. The all-pixel results retain boundary errors: the regenerated healthy image has maximum absolute error 489.6180 stored units across the full image, versus 29.9796 in the interior. The figures should not be interpreted as error-free reconstructions.

After applying the same corrected reconstruction to both insertion methods, their full-image difference has RMSE **0.00155361** and maximum absolute difference **0.0174561**. The incremental lesion error, `(projection insertion - healthy reconstruction) - delta`, has lesion-mask RMSE **2.7262**.

## Validation and preservation

MATLAB R2023b executed the workflow on CPU. The uncorrected reference reconstruction first reproduced the saved original exactly. Assertions checked unchanged image-domain replacement outside the lesion, identical paired support, finite outputs, reduced interior errors, matched-method agreement, and numerical additivity. The regenerated figures were inspected visually, including lesion detail. Exported arrays were read back; copied result files were checked byte-for-byte. The GPU path was not tested.

All pre-existing code, kernels, inputs, and results, including `li_fbp/`, were preserved byte-for-byte. This directory contains the new workflow and outputs; the preceding reconstruction, insertion, diagnosis, display, and registration results remain unchanged. No registration image was regenerated because this correction applies to the forward/backprojection experiment.

## Reproduce

From the repository root in MATLAB:

```matlab
addpath('images/speckle_reduced_case_b_del_slice071');
regenerate_previews('compute');
regenerate_previews('reconstruction');
regenerate_previews('insertion');
regenerate_previews('differences');
regenerate_previews('arrays');
regenerate_previews('metrics');
```

Computation writes intermediate files to `fullfile(tempdir,'cbct_speckle_reduced')`. Each export stage writes one result file. Optional second and third arguments select a computation cache and export destination. The saved MAT arrays use the reference matrix orientation; transpose them to match these figures. `healthy_error` is reconstruction minus input, while the subtraction panel in `01_reconstruction.png` shows its negative, input minus reconstruction.
