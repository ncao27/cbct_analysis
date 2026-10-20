# Repository working practices

These instructions carry forward the owner's playbook at
`~/Desktop/agent_playbooks/cbct_analysis.md`.

- Treat `li_fbp/` as the golden reference: its curved and linear detector
  forward/backprojection scripts were used for the conference abstract.
  Preserve this directory unless the owner explicitly authorizes a change.
- Preserve all C++ and CUDA kernels, now under `kernels/`, and the original
  research inputs and results.
- Do not change geometry, detector, attenuation, sampling, or DICOM-derived
  parameter values as part of cleanup or refactoring.
- Read and review the complete intended change before editing. Implement one
  file change at a time and commit that file before changing the next file.
  Treat a file relocation as one change and commit both paths together.
- Before committing, run the affected code and verify its outputs visually or
  mathematically. Save new figures in a subdirectory of `images/`. For moves
  and deletions, also verify retained code/data and moved files byte-for-byte.
  If the hardware or inputs prevent a reconstruction run, state the limitation;
  do not claim numerical reconstruction validation from static checks.
- The research goal is forward projection of healthy and lesion images, followed
  by image-domain or projection-domain lesion insertion. Lesion location is not
  currently the priority. The reference scripts have not yet been integrated
  with the retained lesion/visualization helpers.
- Removed experiments remain available in Git history; see `docs/cleanup.md`.
