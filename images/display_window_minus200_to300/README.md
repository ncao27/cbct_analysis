# Display previews: -200 to 300

These are new display versions of the saved reconstruction, lesion-insertion,
and registration results. They use the requested window **[-200, 300]**
(width 500, midpoint 50), with the same mapping for every compared CT image:

```text
displayed value = stored pixel value - 8192
grayscale = clip((displayed value + 200) / 500, 0, 1)
```

**8192 is an assumed display offset, not a confirmed HU conversion.** The owner
confirmed that the raw export's conversion settings are unavailable and asked
us to proceed with the best available display method. These figures explicitly
label the assumption. The raw files and saved arrays have no verified HU slope
or intercept. A numerical shift alone cannot establish Hounsfield units.
[DICOM defines the calibrated mapping using rescale slope and intercept](https://dicom.nema.org/medical/dicom/2026a/output/chtml/part03/sect_C.8.2.html).

For this assumption, the display window corresponds to stored values
**7992 through 8492**. Values outside it saturate to black or white for display;
the saved numerical arrays are not clipped or shifted. Zero-valued padding is
black, not interpreted as water. Unsupported registration samples remain NaN.
The narrower window makes fine tissue differences and reconstruction error more
visible than the previous 7000-9000 stored-value window.

## Figures

- [01_reconstruction.png](01_reconstruction.png): healthy image, reconstruction,
  and their direct subtraction.
- [02_lesion_insertion.png](02_lesion_insertion.png): both insertion methods,
  matched reconstruction, and lesion-region zooms.
- [03_inserted_lesion.png](03_inserted_lesion.png): the exact Case A
  `S13_I137_DEL_19x18.raw` lesion patch used in the experiment.
- [04_registration.png](04_registration.png): PV/DEL before and after
  registration, overlays using the same window, and displacement magnitude.
- [05_registration_controls.png](05_registration_controls.png): the local 3D
  slab and known-motion control using the same display mapping.

Subtractions retain their zero-centered difference scale; displacement retains
pixel units. The common offset cancels in subtraction, so MAE, RMSE and bias do
not change when that same offset is applied to both compared arrays. Existing
metrics were not recomputed or relabeled as HU. Relative errors normalized by
absolute image values depend on the offset and must not be interpreted as
HU-relative errors.

Original figures, MAT files, CSVs, metrics, raw inputs, and reconstruction and
registration code are retained unchanged. The original validation exporters
still reproduce their original stored-value figures; use this renderer for
the new display previews. This does not change any projection, reconstruction,
insertion, registration, or DICOM-derived parameter.

## Reproduce

Run from this folder in MATLAB:

```matlab
render_display_previews('check');
render_display_previews('reconstruction');
render_display_previews('insertion');
render_display_previews('lesion');
render_display_previews('registration');
render_display_previews('controls');
```

The optional second argument changes the assumed display offset; it remains
labeled as an assumption. Confirmed calibration should be implemented using
the actual slope and intercept rather than relabeling these previews.

Validation checked mapping endpoints, axes limits, displayed matrix values,
NaN support, and subtraction invariance. All five PNGs were visually inspected.
A SHA-256 audit verified that all 98 pre-existing tracked/data files remained
byte-for-byte unchanged.
