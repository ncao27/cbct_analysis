# Repository cleanup record

Pre-cleanup commit: `17ed43cfae77771205945e97ec47ac6d57acb50a`.

## Scope

Removed 17 superseded or unused files: ten MATLAB implementation scripts,
six editor autosaves, and one empty script. Preserved the independent lesion
and visualization helpers. Moved four native kernel sources and three Windows
MEX binaries without editing their contents. All files under `li_fbp/`, research
datasets, figures, videos, and supporting documents remain at their original
paths with identical contents. No retained numerical code or parameter changed.

## Removed files

- `cbct_pipeline_accelerated.m`
- `fbct_pipeline.m`
- `+siddon/siddoncpp_two.m`
- `+siddon/siddoncpp_three.m`
- `+siddon/siddoncu.m`
- `+recon/fdk.m`
- `+recon/fbp.m`
- `past_matlab_iterations (unused)/cbct_head_phantom_analysis1.m`
- `past_matlab_iterations (unused)/cbct_pipeline_standard.m`
- `past_matlab_iterations (unused)/faster_recon_practice.m`
- `cbct_head_phanton_analysis.asv`
- `cbct_pipeline_accelerated.asv`
- `fbct_pipeline.asv`
- `+viz/double.asv`
- `past_matlab_iterations (unused)/cbct_head_phantom_analysis.asv`
- `past_matlab_iterations (unused)/faster_recon_practice.asv`
- `untitled.m`

Every removed file was tracked at the pre-cleanup commit. Inspect any of them
with `git show '17ed43cfae77771205945e97ec47ac6d57acb50a:ORIGINAL_PATH'`.

## Relocated files

| Original path | New path | SHA-256 (unchanged) |
| --- | --- | --- |
| `siddon_kernel.cpp` | `kernels/cpp/siddon_kernel.cpp` | `438e92b4f1cae9c29d0903ea03f71755ff0dbc89102a42bc80fd5e531355b009` |
| `siddon_kernel_two.cpp` | `kernels/cpp/siddon_kernel_two.cpp` | `a4c7a1b7d61f4d5ecdd55320d6f3a0ca838f7f85425d6cc60d7c08025a4ae81d` |
| `siddon_kernel.cu` | `kernels/cuda/siddon_kernel.cu` | `873a462ad9ea7666cdc6ee1a8ef197880c3a8f927229753934aa276634a58c4d` |
| `fdk_kernel.cu` | `kernels/cuda/fdk_kernel.cu` | `c341b327c0572150d936a244da484ad944c2f9d2aa5fee280e8947e27491c139` |
| `siddon_kernel.mexw64` | `kernels/prebuilt/windows-x64/siddon_kernel.mexw64` | `66ef485e83000cf5dcc95bd0de7d7c98812f19c89fe32d46d055d04d7845197c` |
| `siddon_kernel_two.mexw64` | `kernels/prebuilt/windows-x64/siddon_kernel_two.mexw64` | `b4146c17fea82cbaea0c47451da28ab6ae321bed69f7280d7d91c0cd452db03f` |
| `fdk_kernel.mexw64` | `kernels/prebuilt/windows-x64/fdk_kernel.mexw64` | `7127e939d4a435a86ffdf029b4f287c550abdee5737c0e22661ef48520ef5ccb` |

## Verification

- Before each commit, compare retained files and every relocated file against
  their original contents, and check that only the intended file is staged.
- MATLAB R2023b Update 11 on Apple silicon: retained `lesion.sphere` and both
  visualization helpers passed runtime checks. A 20 × 20 × 20 mask with voxel
  spacing [1, 2, 4] and radius 4 had the correct center and included/excluded
  physical boundaries in all three axes. Its center row contained ones at
  positions 6 through 14 and zeros elsewhere. Both plotting helpers displayed
  the supplied matrix unchanged; `viz.double` retained square image pixels.
- The original big-endian `li_fbp/without_lesion.raw` loaded as a finite
  512 × 512 image, with range [-0.0021000004, 0.063360006].
- Full reference projection/reconstruction and CUDA execution were **not run**:
  this Mac reports zero supported GPUs. The existing Windows MEX binaries were
  not rebuilt or executed. Byte equality establishes preservation, not numerical
  correctness of those experiments or the reference algorithms.

The existing personal introduction in `readme.md` is preserved. The README now
identifies the conference workflow, its workspace/file dependencies, and the
new kernel locations. Ignore rules cover autosaves, macOS metadata, new MEX
builds, and the reference scripts' generated outputs while retaining the input
image and original Windows binaries. `AGENTS.md` carries the owner's working
practices into future repository sessions.
