# A Cone-Beam Computed Tomography (CBCT) Analysis Pipeline

## 🌐 About Me
Hi! I'm Nathan Cao, a current junior at Rice University. I am an electrical and computer enginering major, and outside of academics I pursue medical imaging research at MD Anderson Cancer Center under the guidance of Dr. Ke Li. Of course, research isn't my only extracurricular endeavor. I am a part-time triathlete (hoping to compete at the upcoming Ironman 70.3 Waco), ex-competitive badminton player, and a coffee addict (love a good flat white). Post-graduation I hope to pursue a Ph.D. in medical imaging and one day develop low-radiation imaging modalities suitable for underserved regions. 

---

## 🧠 Problem Statement
CT is one of the most widely used diagnostic and intervention imaging tools. From trauma to brain surgeons, CT is the modality that most rely on. Why? It's simple. It's inexpensive (making it a crowd pleaser among the general population), and fast (so anyone experiencing trauma can get scanned in a timely manner). The project first started when I commenced my research journey with Dr. Li. I had partaken in imaging research before, but because my primary duties involved building simulations, I never got a taste of learning about the field from a classical approach. Dr. Li very generously devoted his time every week to guiding me up the steep learning curve of everything CT, and this project is the result. I hope whoever stumbles across this (very heavily commented) repository uses it however they wish, whether to support their own project ambitions or as a guide. 

---

## Reference workflow

`li_fbp/` contains the forward projection and filtered backprojection code used
for the conference abstract. It is the golden reference for this photon-counting
CT project. Both detector geometries are retained unchanged:

| Detector | Forward projection | Reconstruction |
| --- | --- | --- |
| Curved | `c1_forward_projection_v6_curved_detector.m` | `c2_fbp_v6_curved_detector.m` |
| Linear | `c1_forward_projection_v6_linear_detector.m` | `c2_fbp_v6_linear_detector.m` |

On the Windows MATLAB setup, make `li_fbp` the current folder and run the forward
projector followed by the reconstruction for the **same detector geometry**.
Keep the same MATLAB workspace between the scripts: the curved reconstruction
uses `sdd` defined by its forward projector. These scripts use GPU arrays and
require MATLAB with Parallel Computing Toolbox and a supported GPU.
The forward projectors retain their original Windows-style input/output paths.

The included `without_lesion.raw` is read as a big-endian 512 × 512 single-precision
image. Each forward projector writes `projections.raw`; each reconstruction reads
that file and writes `img.raw`. Both detector pairs share those output names, so
copy results to a separate experiment folder before running the other geometry.
Save exported figures under a subdirectory of `images/`.

All geometry, sampling, attenuation, and DICOM-derived parameters are preserved.
The reference scripts are standalone; they do not call the experimental native
kernels or the MATLAB helpers below.

## Repository layout

| Path | Contents |
| --- | --- |
| `li_fbp/` | Unmodified conference reference scripts and input image |
| `kernels/cpp/` | Original 2D and 3D C++ Siddon kernel sources |
| `kernels/cuda/` | Original CUDA Siddon and FDK kernel sources |
| `kernels/prebuilt/windows-x64/` | Preserved Windows MEX binaries |
| `+lesion/` | Reusable spherical lesion mask helper |
| `+viz/` | Reusable image and sinogram display helpers |
| `ct_volume/` | Original MATLAB phantom datasets |
| `images/` | Existing research figures and subdirectories for new figures |
| `reconstruction videos/` | Existing reconstruction recordings |
| `helpful_documents/` | Existing MEX/CUDA setup notes and screenshots |

Add the repository root to the MATLAB path to call `lesion.sphere`, `viz.double`,
and `viz.sinogram`. These helpers are retained for future lesion insertion work;
they are not an integrated lesion insertion pipeline in the reference workflow.

## Native kernels

The C++ and CUDA experiments are preserved byte-for-byte, separately from the
conference reference. The old MATLAB pipelines and wrappers have been removed.
The `.mexw64` files are preserved Windows binaries, not macOS builds; their build
settings and correspondence to the source variants have not been revalidated.
The C++ and CUDA Siddon sources share the `siddon_kernel` name, so keep their
compiled outputs in separate directories when revisiting these experiments.

## Cleanup and recovery

Superseded MATLAB pipelines, unused wrappers, editor autosaves, and the empty
`untitled.m` were removed in individual commits. Research data, figures, videos,
and supporting documents were retained. See [the cleanup record](docs/cleanup.md)
for the file inventory and verification results.

The complete pre-cleanup repository is available at commit
`17ed43cfae77771205945e97ec47ac6d57acb50a`. To inspect a removed file without
changing the working tree, use:

```sh
git show '17ed43cfae77771205945e97ec47ac6d57acb50a:fbct_pipeline.m'
```

Use the same command with any original path listed in the cleanup record.
