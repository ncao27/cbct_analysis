# Direct projection and FBP: research validation

Use this folder for the current validation. For the lesion-insertion research, start with the corrected original-bandwidth output shown in the middle column of Figure 02 and Figure 05; Figure 01 instead highlights the smoother bandwidth comparison. The earlier support-mask/8192-baseline workaround is excluded. No input pixel is selectively offset, rescaled, replaced, or restored to improve the reconstruction. The previously authorized lesion replacement remains the same independent experiment.

## Results to open

- [01: Healthy input, original reference, and updated reconstruction](01_reconstruction.png)
- [02: Paired image-domain and projection-domain lesion insertion](02_lesion_insertion.png)
- [03: Complete pixel-error maps and lesion differences](03_pixel_differences.png)
- [04: Independent analytic phantom validation](04_phantom_validation.png)
- [05: Reconstruction-bandwidth comparison](05_bandwidth_comparison.png)
- [06: Metrics, numerical checks, and assumptions](06_metrics.json)
- [07: Numerical images, projections, error maps, and reference-source snapshots](07_arrays.mat)
- [08: File hashes and provenance](08_provenance.json)

All patient-image units are **stored pixel units, not confirmed HU or physical attenuation**. The input is Case B DEL slice 71; the lesion is Case A `S13_I137_DEL_19x18.raw`, with the same 246-pixel replacement mask and location as before. Raw encoding is still inferred big-endian uint16 because file length conflicts with the reported float32 metadata. The comparison tests a numerical round trip, not a calibrated acquisition model or clinical accuracy.

## Retained model and explicit changes

The method remains the professor's sampled-ray cubic interpolation followed by equiangular fan-beam FBP. It is not Siddon. The original `li_fbp/` files and kernels were not modified. SID 750 mm, SDD 1200 mm, 2480 detector elements, 0.35 mm detector pitch, 2000 views, 400 mm reconstruction FOV, and the 512 x 512 grid remain fixed.

1. **Detector interpolation and filter domain:** divide the angular displacement by the angular bin spacing to obtain a fractional interpolation weight. Evaluate the full filtered convolution wherever the output image needs it, including beyond the original raw-detector interval. Zero extension for the convolution adds no measured detectors and preserves all 2480 physical samples and their coordinates. Setting filtered values outside the raw interval to zero was incorrect: filtered projections can be nonzero there even when raw projections are zero. Geometric validity guards remain; these are not masks applied to input intensities.
2. **Ray integration:** find the ray's entry/exit at the rectangular interpolation domain and integrate inside it. The same midpoint rule and cubic interpolation are retained. Step refinement is a numerical convergence test; scanner geometry and sampling locations on the detector stay fixed.
3. **Exterior boundary model:** use zero-valued samples outside the original array so the cubic interpolant has a defined extension. Four exterior rows/columns are used; six give equivalent projections to the stated tolerance. Original samples and their physical coordinates are checked for exact preservation. This is an explicit change from the reference's abrupt zero extrapolation at the outermost sample centers. It introduces no anatomical segmentation and does not subtract a baseline.
4. **Reconstruction bandwidth:** retain the fan-beam ramp-filter form, but set its angular cutoff to `SID/(2*dx) = 480 cycles/radian`, corresponding to `1/(2*dx) = 0.64 cycles/mm` at isocenter. This value comes from the output pixel spacing, not a fit to the target image. The original detector-Nyquist bandwidth is retained as a control. Bandwidth changes spatial resolution and ringing; it is not a claim of unchanged sharpness everywhere.
5. **Arithmetic:** use double precision throughout the new numerical calculation. No post-reconstruction image denoising, ground-truth error-map subtraction, fitted gain, or computational intensity offset is applied.

The forward model is the line integral of a cubic-convolution interpolant of the unchanged image samples with zero extension. For healthy image h, replacement difference delta, projector A, and the declared FBP B:

```text
healthy reconstruction       = B(A(h))
projection-domain insertion  = B(A(h) + A(delta))
image-domain control         = B(A(h + delta))
```

No original-image support mask enters these formulas. The image-domain control is forward-projected independently rather than constructed by adding the other two projection arrays.

The last point of the original interpolation grid is not the same as an image cell's outer edge. The zero-extension model explicitly allows the cubic basis to extend outside the original sample-center rectangle. It preserves the reconstruction grid but changes the continuous boundary representation. It is a numerical boundary condition, not measured information about anatomy outside the scan.

## Accuracy and remaining limitations

RMSE in stored units:

| Healthy-image method | All 262,144 pixels | Original interior ROI | Central radius 150 mm |
| --- | ---: | ---: | ---: |
| Original professor reference | 289.52895 | 67.29781 | 64.93803 |
| FBP corrections with original projections | 257.21642 | 46.41103 | 44.89129 |
| Boundary/integration corrections, original bandwidth | 21.97551 | 15.28276 | 17.19178 |
| Boundary/integration corrections, grid bandwidth | 64.85631 | 11.01359 | 4.56963 |

The original interior ROI contains 244,596 nonzero pixels more than five pixels from image edges/zero padding; it includes background. The central ROI is a fixed geometric circle, not a tissue segmentation. These masks measure error only. The reconstruction is never replaced or masked using them. Full-image results remain the primary error accounting.

The updated healthy reconstruction has full-image maximum absolute error **1219.59931** and interior maximum **140.10858**. Significant boundary error remains. A visually cleaner image is not a zero-error reconstruction.

For projection insertion versus direct image-domain replacement, updated RMSE is **64.85957** over the whole image, **11.03417** in the interior, and **15.33304** within the lesion. The incremental lesion error, `[projection insertion - healthy reconstruction] - delta`, has lesion RMSE **16.23050**.

After both insertion paths undergo the same reconstruction, their full-image RMSE is **3.9589917e-12**, maximum absolute difference **5.2750693e-11**. This establishes numerical agreement of the simulated methods, not clinical validity.

## Bandwidth tradeoff: do not treat smoothness as lesion accuracy

With the numerical boundary and FBP corrections, the original-bandwidth projection insertion has lesion-mask RMSE **20.26008**. Its incremental lesion RMSE is **2.72620**, compared with **16.23050** for the grid-bandwidth result. The narrower bandwidth improves background artifacts and total lesion error here, but increases the error of the inserted lesion signal itself, particularly near its hard boundary.

The corrected original bandwidth also has lower whole-image healthy RMSE (**21.97551**) than the grid-bandwidth result (**64.85631**). The narrower filter reduces central-image artifacts while increasing boundary error. For this lesion experiment, the corrected original bandwidth is the more faithful starting point; the smoother result is a declared bandwidth comparison.

Both outputs are retained. Figure 02 shows direct image replacement, original-bandwidth projection insertion, and grid-bandwidth projection insertion. `unwindowed_healthy` and `unwindowed_insertion` retain the original ramp bandwidth after the numerical corrections; `healthy_reconstruction` and `projection_domain` use the grid-derived bandwidth. Figure 01's updated panel uses the latter and must not be interpreted as an unchanged-resolution improvement. There is no claim that one filter is best for every research endpoint.

## Independent and convergence checks

- Analytic off-center Gaussian: numerical forward-projection relative L2 error **5.964899e-08**; FBP of exact projections relative image L2 error **1.0570452e-05**. Acceptance limits were 1e-4 and 0.005, respectively. No gain was fitted.
- Analytic unit disk: interior reconstruction RMSE **0.0045643599**, excluding a three-pixel edge band; acceptance limit 0.01 of unit amplitude. Full-image RMSE **0.015766147** includes edge ringing/partial-volume effects.
- Unit square: original midpoint path-length RMS error **0.094937566 mm** on the tested views, versus **1.4898663e-14 mm** for integration with exact entry/exit bounds. The latter was required to have maximum error below 1e-9 mm.
- Maximum integration step: **0.2403846154 mm** versus **0.1201923077 mm** over all 2000 views. Projection relative L2 change **7.9717426e-07**; reconstructed full-image RMS change **0.017899509** stored units. Limits were 1e-5 and 0.1.
- An additional step halving on 17 evenly distributed views changed projections by **1.6478877e-07** relative L2 (limit 1e-5).
- Two versus four exterior zero rows differed by **0.00026865527** relative L2, failing the padding criterion. Four versus six differed by **1.9049644e-08**, passing 1e-5. The failed control is retained rather than hidden.
- Independently projected image-domain and projection-domain inputs agree to relative L2 **3.8646924e-17** (limit 1e-12). Backprojection additivity maximum discrepancy is **7.2306835e-11** stored units (limit 1e-6).

The unit-disk right-edge 10-90% transition width is **1.07931 mm** at the original bandwidth and **1.06695 mm** at the grid-derived bandwidth. This is one sampled edge-profile measurement, not a complete MTF characterization. Fan-beam spatial bandwidth varies with position; the stated cutoff is defined at isocenter.

## Research basis

The equiangular weighting, ramp form, and full-scan factor follow the formulation reviewed in [Fessler, Analytical Tomographic Image Reconstruction Methods, section 3.9.2.1](https://web.eecs.umich.edu/~fessler/book/c-tomo.pdf), particularly equations 3.9.29-3.9.32. The implementation uses the band-limited angular ramp with cutoff b, multiplied by the same fan-beam correction `(gamma/sin(gamma))^2`; its zero-lag value is b^2. This changes the reconstruction bandwidth while preserving the fan geometry.

[Kak and Slaney, chapter 5](https://engineering.purdue.edu/~malcolm/pct/CTI_Ch05.pdf) explains how sharp boundaries, finite projection sampling, and reconstruction sampling produce aliasing patterns. That motivates explicit convergence and bandwidth controls rather than an image-derived correction. [MATLAB interp2 documentation](https://www.mathworks.com/help/matlab/ref/interp2.html) specifies cubic convolution and the treatment of samples outside the interpolation domain; that distinction motivates the boundary tests. These references support the numerical approach, not patient-specific calibration claims.

## Display, preservation, and reproduction

The uniform shift by 8192 occurs **only in plotting**, with limits [-200, 300]. It is the previously requested display convention, not a numerical reconstruction offset or confirmed HU conversion. Error arrays and metrics use original stored units; displayed color limits do not clip saved arrays. Matrix arrays follow the reference orientation; transpose them to match the figures.

MATLAB R2023b ran on CPU. GPU behavior was not validated. The original reference output had already been reproduced bit-for-bit using the control workflow. All pre-existing input, kernel, reference-code, and numerical-result files were preserved. The earlier baseline-based results remain historical diagnostics and must not be used as direct projection/FBP validation.

From the repository root, use the research workflow below. `validate_direct_reference.m` is a separate reference-only control; if exporting that control, supply a distinct destination so its older filenames do not overwrite these results.

```matlab
addpath('images/direct_reference_case_b_del_slice071');
validate_reference_numerics('compute');
for stage = {'reconstruction','insertion','differences','phantoms','bandwidth','arrays','metrics'}
    validate_reference_numerics(stage{1});
end
```

The optional second argument selects the computation cache; the optional third selects an export destination. The default cache is `fullfile(tempdir,'cbct_reference_numerics')`. The saved MAT file includes original reference-source snapshots, quantitative arrays, projections, test results, and analytic phantoms. The metrics JSON and provenance manifest record the reported run; the README is a static report of that run.
