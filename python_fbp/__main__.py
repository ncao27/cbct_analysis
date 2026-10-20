"""Run either reference detector geometry without changing li_fbp."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from .fbp import (
    Geometry, REFERENCE_DIR, forward_project, read_image, read_projections,
    reconstruct, write_image,
)

ROOT = Path(__file__).resolve().parents[1]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("forward", "reconstruct", "run"))
    parser.add_argument("--detector", required=True, choices=("linear", "curved"))
    parser.add_argument("--backend", default="numpy", choices=("numpy", "cupy"))
    parser.add_argument("--input", type=Path, help="Image raw for forward/run; projection raw for reconstruct")
    parser.add_argument("--output-dir", type=Path, help="Default: python_fbp/results/<detector>")
    parser.add_argument("--chunk-rays", type=int, default=64, help="Forward ray batch size; does not change geometry")
    parser.add_argument("--plot", action="store_true", help="Save the MATLAB-style reconstruction display under images/")
    args = parser.parse_args(argv)
    if args.chunk_rays < 1:
        parser.error("--chunk-rays must be positive")
    if args.command == "forward" and args.plot:
        parser.error("--plot requires reconstruction (use run or reconstruct)")
    output = (args.output_dir or ROOT / "python_fbp" / "results" / args.detector).resolve()
    if output == REFERENCE_DIR or REFERENCE_DIR in output.parents:
        parser.error("--output-dir must be outside the preserved li_fbp directory")
    geometry = Geometry()
    projections_path = output / "projections.raw"
    input_path = args.input or (projections_path if args.command == "reconstruct" else REFERENCE_DIR / "without_lesion.raw")

    def progress(done, total):
        if done == 1 or done % 50 == 0 or done == total:
            print(f"View {done}/{total}", flush=True)

    try:
        if args.command in ("forward", "run"):
            image = read_image(input_path, geometry)
            projections = forward_project(image, args.detector, geometry,
                                          backend=args.backend, chunk_rays=args.chunk_rays,
                                          output=projections_path, progress=progress)
            print(f"Projections: {projections_path}")
        else:
            projections = read_projections(input_path, geometry)
        if args.command in ("reconstruct", "run"):
            image = reconstruct(projections, args.detector, geometry, backend=args.backend, progress=progress)
            image_path = output / "img.raw"
            write_image(image_path, image)
            print(f"Reconstruction: {image_path}")
            if args.plot:
                import matplotlib
                matplotlib.use("Agg")
                import matplotlib.pyplot as plt
                figure, axes = plt.subplots()
                artist = axes.imshow(image.T, cmap="gray", vmin=0.015, vmax=0.035,
                                     origin="upper", interpolation="nearest")
                axes.set_title(f"FBP Reconstruction ({args.detector}, Python)")
                figure.colorbar(artist, ax=axes)
                directory = ROOT / "images" / "python_fbp" / args.detector
                directory.mkdir(parents=True, exist_ok=True)
                filename = directory / "reconstruction.png"
                figure.savefig(filename, dpi=160, bbox_inches="tight")
                plt.close(figure)
                print(f"Figure: {filename}")
        output.mkdir(parents=True, exist_ok=True)
        manifest = dict(command=args.command, detector=args.detector, backend=args.backend,
                        input=str(input_path.resolve()), geometry=vars(geometry),
                        projection_shape=[geometry.views, geometry.detectors],
                        projection_layout="views x detectors, little-endian float32, C order",
                        image_layout="MATLAB matrix, little-endian float32, Fortran order")
        (output / "run.json").write_text(json.dumps(manifest, indent=2) + "\n")
    except (OSError, ValueError, RuntimeError, ImportError) as error:
        parser.exit(1, f"Error: {error}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
