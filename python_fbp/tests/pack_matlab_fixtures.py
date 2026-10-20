"""Package the raw outputs of generate_matlab_fixtures.m as a regression fixture."""
from pathlib import Path
import argparse
import hashlib
import json

import numpy as np


def pack(directory, destination, repository):
    directory, repository = Path(directory), Path(repository)
    payload = {"image": np.fromfile(directory / "fixture_input.raw", ">f4").reshape((12, 12), order="F").astype(np.float64)}
    for detector in ("linear", "curved"):
        for name, shape, order in (("projections", (9, 16), "C"),
                                   ("filtered", (9, 16), "F"), ("img", (12, 12), "F")):
            payload[f"{detector}_{name}"] = np.fromfile(directory / f"{detector}_{name}.raw", "<f4").reshape(shape, order=order)
    for name, shape in (("v", (5, 6)), ("xq", (7, 12)), ("yq", (7, 12)), ("z", (7, 12))):
        payload[f"interp_{name}"] = np.fromfile(directory / f"interp_{name}.raw", "<f8").reshape(shape, order="F")
    payload["production_view_indices"] = np.array([0, 250, 500, 750])
    for detector in ("linear", "curved"):
        payload[f"{detector}_production_projections"] = np.fromfile(directory / f"{detector}_full_selected.raw", "<f4").reshape(4, 2480)
    payload["reference_image_sha256"] = np.array(hashlib.sha256((repository / "li_fbp/without_lesion.raw").read_bytes()).hexdigest())
    sources = {}
    for path in sorted((repository / "li_fbp").glob("*.m")):
        sources[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    if len(sources) != 4:
        raise ValueError("Expected all four preserved MATLAB reference scripts")
    payload["source_sha256"] = np.array(json.dumps(sources, sort_keys=True))
    payload["provenance"] = np.array("MATLAB R2023b Update 11, Apple silicon CPU. GPU/display plumbing replaced in memory; synthetic size=12, detectors=16, views=9, du=35; other parameters unchanged. Also includes four production-size projections at original view indices [0,250,500,750] for each detector.")
    for key, value in payload.items():
        if value.dtype.kind in "fc" and not np.isfinite(value).all():
            raise ValueError(f"Non-finite fixture: {key}")
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(destination, **payload)
    print(f"Packed {len(payload)} arrays/metadata entries into {destination}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--repository", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("matlab_cpu_reference.npz"))
    args = parser.parse_args()
    pack(args.directory, args.output, args.repository)
