# Reconstruction speckle diagnosis

Case B DEL slice 71, investigated 2026-11-19. The excess fine bright/dark texture is reproducible reconstruction error. The lesion-insertion pipeline inherits the same healthy reconstruction artifact. No random noise was added by either validation pipeline.

[Controlled comparison](01_speckle_diagnosis.png)

All errors below use stored pixel units, not confirmed HU. The original input encoding remains inferred big-endian uint16; export metadata and HU/attenuation calibration remain unavailable. Image panels use stored value minus 8192 with limits [-200, 300], solely as the previously authorized display convention. Error panels use reconstruction minus input with limits [-200, 200].

## Confirmed interpolation defect

In `li_fbp/c2_fbp_v6_curved_detector.m`, line 109 computes:

```matlab
w = gamma0 - gamma(k);
```

Lines 113-114 use this as a dimensionless interpolation fraction. It should be divided by the detector angular spacing:

```matlab
w = (gamma0 - gamma(k)) ./ d_gamma;
```

The spacing is approximately 0.000291667 radians. The existing weight therefore stays near zero instead of spanning 0 to 1 between detector bins, almost always selecting the lower bin. Correcting only this expression in an in-memory CPU copy reduces interior RMSE from 67.2978 to 46.4115. The original file was not changed.

## Large stored-value baseline reproduces the artifact

The published test forward-projects raw stored values directly, with no conversion to physical attenuation. The forward projector uses cubic interpolation with zero extrapolation at line 87 and midpoint ray integration at lines 47-49 and 91. Thus the large baseline inside the image support ends abruptly at its boundary. The backprojector applies the unwindowed ramp at lines 51-73.

To isolate the baseline contribution, define h as the healthy input, A as the original forward projector, R as the original filtered backprojector, and p = 8192 * (h ~= 0). Forward-project and reconstruct p with the same geometry and sampling. The error R(A(p)) - p has correlation 0.998117 with the original healthy reconstruction error inside the original interior ROI. This establishes that a constant-valued support image reproduces almost the same error pattern, without anatomical texture.

A second diagnostic subtracts the baseline's projections, reconstructs the residual, and restores the known baseline exactly:

```text
bypass = R(A(h) - A(p)) + p
```

Its interior RMSE is 9.5668. With normalized detector interpolation it is 6.3099. These experiments support a large baseline interacting with discretized projection/reconstruction as the dominant source of the observed pattern. Sampling/interpolation errors at the abrupt support boundary, emphasized by ramp filtering, are a plausible mechanism; this experiment does not separately apportion the remaining error among ray integration, detector sampling, filtering, and backprojection.

**The bypass uses the original image support and restores p exactly. It is diagnostic evidence, not a validated physical reconstruction, a recommended unconditional subtraction of 8192, or proof that 8192 is the true HU offset.** Display-only subtraction cannot remove artifacts generated earlier during reconstruction. No image denoising was applied.

## Numerical controls

All metrics in this table use the original saved interior mask (244,596 pixels). It excludes zero pixels and pixels within five pixels of the image edge or zero padding; it includes nonzero air/background, so it is not an anatomical tissue-only ROI.

| Experiment | Interior MAE | Interior RMSE | High-frequency error RMS |
| --- | ---: | ---: | ---: |
| Published reconstruction, reproduced exactly | 53.4378 | 67.2978 | 59.2177 |
| Normalize detector interpolation only | 36.9690 | 46.4115 | 41.1306 |
| Higher-precision filtering/backprojection control | 53.4368 | 67.2958 | 59.2164 |
| Higher precision plus normalized interpolation | 36.9690 | 46.4115 | 41.1306 |
| Reconstruct support-only p, compared with p | 60.1878 | 75.8016 | 66.7291 |
| Known-baseline bypass | 7.5709 | 9.5668 | 8.2979 |
| Known-baseline bypass plus normalized interpolation | 5.0217 | 6.3099 | 5.5314 |

High-frequency error means e - imgaussfilt(e, 1), measured only for analysis; this operation was not applied to reconstructed images. The precision control promotes projection/filter/backprojection arrays and accumulators, but retains original float32 projections and some single-precision geometry/kernel intermediates. It does not establish the behavior of a fully double-precision forward-and-back pipeline.

The original baseline was reproduced bit-for-bit before running variants. The RMS residual in the decomposition R(A(h)) - R(A(p)) - R(A(h)-A(p)) is 0.0188663 inside the ROI, supporting approximate numerical linearity at this scale.

## Why projection-domain insertion shows the same flakes

`images/lesion_insertion_case_b_del_slice071/validate_lesion_insertion.m`, lines 66-70, reconstructs the healthy projections plus the projected lesion replacement difference. With linear projection/reconstruction, the healthy reconstruction error remains a common component.

Outside the lesion, within the same interior ROI, the healthy error and projection-domain insertion error have correlation **0.9999906153**. Their difference has RMS **0.2915655** stored units. This explains why the pattern appears in both outputs; it does not establish absence of lesion-specific errors.

## Additional observations and scope

The invalid detector-index handling in backprojection lines 105-114 sets invalid indices to 1 and weights to zero, then still contributes detector bin 1. Affected contributions should be masked. There are 2,056 pixels with at least one invalid ray, but none lie in the measured interior ROI, so this does not explain the interior speckle.

The earlier display changed from a 2,000-unit window to a 500-unit window, making the same intensity differences four times more prominent in the unsaturated grayscale range. This reveals existing errors; it does not create them.

Tests ran in MATLAB R2023b on CPU copies of the curved-detector reference scripts. Only paths, GPU wrappers, plotting/logging, and the explicitly described experimental variants differed. Detector geometry, sampling, and reference files were preserved. The in-memory backprojection inherited sdd = single(1200), matching the original validation's forward-projector workspace. The native GPU path was not validated on this machine.

Inputs were the saved arrays in `images/lesion_free_roundtrip_case_b_del_slice071/03_roundtrip_arrays.mat` and `images/lesion_insertion_case_b_del_slice071/03_insertion_arrays.mat`. Temporary execution files and detailed numerical arrays were placed in `/private/tmp/cbct_speckle_diagnosis` and are not durable repository artifacts. Only this report and its comparison figure are added to the repository. No production correction has been applied to li_fbp or other existing code.
