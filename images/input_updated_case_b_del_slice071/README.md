# Input and updated reconstruction comparisons

Each figure contains only the input image on the left and an updated reconstruction on the right. These are new layouts of the saved, validated results; the reconstruction was not recomputed or altered. The earlier artifact-heavy reconstruction and historical baseline-corrected previews are not included.

| Image | Reference ramp bandwidth | Grid-derived ramp bandwidth |
| --- | --- | --- |
| Healthy | [01: Input and updated](01_healthy_reference_bandwidth.png) | [02: Input and updated](02_healthy_grid_bandwidth.png) |
| Lesion inserted | [03: Input and updated](03_lesion_reference_bandwidth.png) | [04: Input and updated](04_lesion_grid_bandwidth.png) |

Both bandwidths use the numerical corrections from the [direct reconstruction validation](../direct_reference_case_b_del_slice071/README.md). Here, **reference bandwidth** means the corrected reconstruction retaining the professor's ramp bandwidth, not the old uncorrected reconstruction. The grid-derived bandwidth is the smoother comparison. Their accuracy and resolution tradeoffs remain as documented in that validation.

The healthy input is Case B DEL slice 71. In the lesion figures, the left panel contains direct image-domain replacement and the right panel contains projection-domain insertion followed by reconstruction, using the same lesion and placement. Orange boxes mark the lesion patch.

All panels use the same display convention: stored value minus 8192, window [-200, 300]. This is an assumed display offset, not a confirmed HU conversion. No new smoothing, pixel correction, intensity fitting, or baseline subtraction was performed.

The source is [07_arrays.mat](../direct_reference_case_b_del_slice071/07_arrays.mat). Healthy figures use `healthy` paired with `unwindowed_healthy` or `healthy_reconstruction`; lesion figures use `image_domain` paired with `unwindowed_insertion` or `projection_domain`. The original numerical arrays and existing figures are preserved.

To reproduce from the repository root in MATLAB:

```matlab
addpath('images/input_updated_case_b_del_slice071');
for stage = {'healthy_reference','healthy_grid','lesion_reference','lesion_grid'}
    render_input_updated(stage{1});
end
```

An optional second argument selects another export folder. MATLAB verified two image panels per figure, matching display scales, and unchanged displayed array values apart from the documented transpose and uniform display shift. All four PNGs were visually inspected before saving.
