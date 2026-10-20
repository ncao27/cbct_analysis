"""Numerical regression against MATLAB CPU outputs, plus raw-layout/CLI checks."""
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import asdict
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from python_fbp import (
    Geometry, backproject, cubic_interpolate, forward_project, ramp_filter,
    read_image, read_projections, reconstruct, write_image, write_projections,
)
from python_fbp.fbp import REFERENCE_DIR, array_module

FIXTURE = Path(__file__).with_name("matlab_cpu_reference.npz")
SMALL = Geometry(image_size=12, detectors=16, views=9, du=35)


class ReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = np.load(FIXTURE, allow_pickle=False)

    @classmethod
    def tearDownClass(cls):
        cls.fixture.close()

    def test_default_geometry(self):
        self.assertEqual(asdict(Geometry()), dict(image_size=512, detectors=2480,
                         views=2000, fov=400, sid=750, sdd=1200, du=0.35))
        self.assertEqual(Geometry().angles[0], 0)
        self.assertLess(Geometry().angles[-1], 2 * np.pi)

    def test_reference_source_hashes(self):
        for filename, digest in json.loads(str(self.fixture["source_sha256"])).items():
            self.assertEqual(hashlib.sha256((REFERENCE_DIR / filename).read_bytes()).hexdigest(), digest)

    def test_production_size_selected_views(self):
        image_path = REFERENCE_DIR / "without_lesion.raw"
        self.assertEqual(hashlib.sha256(image_path.read_bytes()).hexdigest(), str(self.fixture["reference_image_sha256"]))
        image = read_image(image_path)
        for detector in ("linear", "curved"):
            with self.subTest(detector=detector):
                actual = forward_project(image, detector, view_indices=self.fixture["production_view_indices"])
                expected = self.fixture[f"{detector}_production_projections"]
                relative = np.linalg.norm(actual - expected) / np.linalg.norm(expected)
                self.assertLess(relative, 2e-6)
                np.testing.assert_allclose(actual, expected, atol=1e-4, rtol=2e-5)

    def test_cubic_edges_and_extrapolation(self):
        actual = cubic_interpolate(self.fixture["interp_v"], self.fixture["interp_xq"], self.fixture["interp_yq"])
        np.testing.assert_allclose(actual, self.fixture["interp_z"], atol=1e-14, rtol=1e-14)

    def test_forward_projection(self):
        for detector in ("linear", "curved"):
            with self.subTest(detector=detector):
                actual = forward_project(self.fixture["image"], detector, SMALL, chunk_rays=5)
                np.testing.assert_allclose(actual, self.fixture[f"{detector}_projections"], atol=2e-5, rtol=2e-6)
                self.assertEqual(actual.dtype, np.float32)

    def test_ramp_filter(self):
        for detector in ("linear", "curved"):
            with self.subTest(detector=detector):
                actual = ramp_filter(self.fixture[f"{detector}_projections"], detector, SMALL)
                np.testing.assert_allclose(actual, self.fixture[f"{detector}_filtered"], atol=1e-5, rtol=1e-5)

    def test_backprojection(self):
        for detector in ("linear", "curved"):
            with self.subTest(detector=detector):
                actual = backproject(self.fixture[f"{detector}_filtered"], detector, SMALL)
                np.testing.assert_allclose(actual, self.fixture[f"{detector}_img"], atol=2e-6, rtol=2e-5)

    def test_complete_pipeline(self):
        for detector in ("linear", "curved"):
            with self.subTest(detector=detector):
                projections = forward_project(self.fixture["image"], detector, SMALL)
                actual = reconstruct(projections, detector, SMALL)
                np.testing.assert_allclose(actual, self.fixture[f"{detector}_img"], atol=3e-6, rtol=3e-5)

    def test_chunking_selected_views_and_streaming(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "projection.raw"
            for detector in ("linear", "curved"):
                selected = forward_project(self.fixture["image"], detector, SMALL,
                                           chunk_rays=3, view_indices=[0, 3], output=path)
                alternate = forward_project(self.fixture["image"], detector, SMALL,
                                            chunk_rays=13, view_indices=[0, 3])
                np.testing.assert_array_equal(selected, alternate)
                np.testing.assert_array_equal(np.fromfile(path, "<f4").reshape(2, 16), selected)
                np.testing.assert_allclose(selected, self.fixture[f"{detector}_projections"][[0, 3]], atol=2e-5, rtol=2e-6)

    def test_invalid_bins_preserve_first_sample(self):
        g = Geometry(image_size=3, detectors=2, views=1)
        img = backproject(np.array([[2, 9]], dtype=np.float32), "linear", g)
        x = -g.fov / 3
        expected = g.sid / (g.sid - x) ** 2 * 2 * np.pi
        self.assertAlmostEqual(float(img[0, 0]), expected, places=8)

    def test_interpolation_weight_is_not_normalized(self):
        g = Geometry(image_size=3, detectors=4, views=1, fov=3, du=2)
        img = backproject(np.array([[1, 3, 7, 11]], dtype=np.float32), "linear", g)
        self.assertAlmostEqual(float(img[1, 1]), 7 * np.pi / g.sid, places=8)

    def test_raw_order_and_endianness(self):
        g = Geometry(image_size=3, detectors=3, views=2)
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            path = directory / "image.raw"
            np.arange(9, dtype=">f4").tofile(path)
            expected = np.array([[0, 3, 6], [1, 4, 7], [2, 5, 8]], dtype=np.float64)
            np.testing.assert_array_equal(read_image(path, g), expected)
            write_image(path, expected)
            np.testing.assert_array_equal(np.fromfile(path, "<f4"), np.arange(9))
            projections = np.array([[1, 2, 3], [4, 5, 6]], dtype=np.float32)
            write_projections(path, projections)
            self.assertEqual(path.read_bytes(), np.arange(1, 7, dtype="<f4").tobytes())
            np.testing.assert_array_equal(read_projections(path, g), projections)
            with self.assertRaises(ValueError):
                read_image(path, g)

    def test_invalid_input_and_protected_output(self):
        with self.assertRaises(ValueError):
            write_image(REFERENCE_DIR / "must_not_exist.raw", np.zeros((3, 3)))
        with self.assertRaises(ValueError):
            forward_project(self.fixture["image"], "unknown", SMALL)
        for indices in ([], [-1], [SMALL.views], [0.5]):
            with self.assertRaises(ValueError):
                forward_project(self.fixture["image"], "linear", SMALL, view_indices=indices)
        with self.assertRaises(ValueError):
            ramp_filter(np.zeros((1, 1)), "linear", SMALL)
        with self.assertRaises(ValueError):
            reconstruct(np.full((9, 16), np.nan), "linear", SMALL)

    def test_cli_pipeline_and_protected_directory(self):
        import python_fbp.__main__ as cli
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / "input.raw"
            self.fixture["image"].astype(">f4").ravel(order="F").tofile(source)
            with patch.object(cli, "Geometry", return_value=SMALL), redirect_stdout(io.StringIO()):
                result = cli.main(["run", "--detector", "curved", "--input", str(source),
                                   "--output-dir", str(directory / "output")])
            self.assertEqual(result, 0)
            output = np.fromfile(directory / "output/img.raw", "<f4").reshape(12, 12, order="F")
            np.testing.assert_allclose(output, self.fixture["curved_img"], atol=3e-6, rtol=3e-5)
            manifest = json.loads((directory / "output/run.json").read_text())
            self.assertEqual(manifest["geometry"]["du"], 35)
            with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as context:
                cli.main(["forward", "--detector", "linear", "--output-dir", str(REFERENCE_DIR)])
            self.assertEqual(context.exception.code, 2)

    def test_optional_cuda_backend(self):
        try:
            array_module("cupy")
        except RuntimeError as error:
            self.skipTest(str(error))
        for detector in ("linear", "curved"):
            actual = forward_project(self.fixture["image"], detector, SMALL, backend="cupy")
            np.testing.assert_allclose(actual, self.fixture[f"{detector}_projections"], atol=3e-5, rtol=2e-5)
            img = reconstruct(actual, detector, SMALL, backend="cupy")
            np.testing.assert_allclose(img, self.fixture[f"{detector}_img"], atol=1e-5, rtol=1e-4)


if __name__ == "__main__":
    unittest.main()
