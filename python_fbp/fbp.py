"""Array-based ports of the four li_fbp MATLAB reference scripts.

Public matrices use (view, detector) for projections and the same row/column
indices as MATLAB for images. Raw image I/O explicitly uses MATLAB column order.
The backprojection's unusual weights and invalid-bin fallback are intentional.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Literal

import numpy as np

Detector = Literal["linear", "curved"]
Progress = Callable[[int, int], None]
REFERENCE_DIR = (Path(__file__).resolve().parents[1] / "li_fbp").resolve()


@dataclass(frozen=True)
class Geometry:
    """Defaults reproduce li_fbp; overrides are useful for synthetic tests."""

    image_size: int = 512
    detectors: int = 2480
    views: int = 2000
    fov: float = 400.0
    sid: float = 750.0
    sdd: float = 1200.0
    du: float = 0.35

    def __post_init__(self):
        for name, minimum in (("image_size", 3), ("detectors", 2), ("views", 1)):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
                raise ValueError(f"{name} must be an integer >= {minimum}")
        for name in ("fov", "sid", "sdd", "du"):
            if not np.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive and finite")

    @property
    def dx(self):
        return self.fov / self.image_size

    @property
    def angles(self):
        return np.arange(self.views, dtype=np.float64) * (2 * np.pi / self.views)


def array_module(backend="numpy"):
    """Choose CPU explicitly, or require a working CuPy/CUDA installation."""
    if backend == "numpy":
        return np
    if backend != "cupy":
        raise ValueError("backend must be 'numpy' or 'cupy'")
    try:
        import cupy as cp
        if cp.cuda.runtime.getDeviceCount() < 1:
            raise RuntimeError("No CUDA device is available")
    except (ImportError, RuntimeError) as error:
        raise RuntimeError("The cupy backend requires CuPy and a supported CUDA GPU") from error
    return cp


def _host(array, xp):
    return np.asarray(array) if xp is np else xp.asnumpy(array)


def _detector(detector):
    if detector not in ("linear", "curved"):
        raise ValueError("detector must be 'linear' or 'curved'")


def _output_path(path):
    path = Path(path).resolve()
    if path == REFERENCE_DIR or REFERENCE_DIR in path.parents:
        raise ValueError("Python outputs must be outside the preserved li_fbp directory")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _read_raw(path, shape, dtype, order):
    path = Path(path)
    expected = int(np.prod(shape)) * np.dtype(dtype).itemsize
    if path.stat().st_size != expected:
        raise ValueError(f"{path}: expected {expected} bytes for shape {shape}")
    data = np.fromfile(path, dtype=dtype).reshape(shape, order=order)
    if not np.isfinite(data).all():
        raise ValueError(f"{path}: input contains non-finite values")
    return data


def read_image(path, geometry=Geometry()):
    """Read big-endian float32, column order, as double (MATLAB fread default)."""
    return _read_raw(path, (geometry.image_size,) * 2, ">f4", "F").astype(np.float64)


def read_projections(path, geometry=Geometry()):
    """Read little-endian float32 into (views, detectors)."""
    return _read_raw(path, (geometry.views, geometry.detectors), "<f4", "C").astype(np.float32)


def write_projections(path, projections):
    """Write each view's detector samples consecutively, matching MATLAB fwrite."""
    np.asarray(projections, dtype="<f4").tofile(_output_path(path))


def write_image(path, image):
    """Write little-endian float32 in MATLAB column order."""
    np.asarray(image, dtype="<f4").ravel(order="F").tofile(_output_path(path))


def _pad_cubic(image, xp):
    # Quadratic ghost samples reproduce interp2's cubic boundary behavior.
    padded = xp.concatenate((3 * image[:1] - 3 * image[1:2] + image[2:3],
                             image, 3 * image[-1:] - 3 * image[-2:-1] + image[-3:-2]), axis=0)
    return xp.concatenate((3 * padded[:, :1] - 3 * padded[:, 1:2] + padded[:, 2:3],
                           padded, 3 * padded[:, -1:] - 3 * padded[:, -2:-1] + padded[:, -3:-2]), axis=1)


def _cubic_weights(t):
    t2, t3 = t * t, t * t * t
    return (-0.5 * t + t2 - 0.5 * t3, 1 - 2.5 * t2 + 1.5 * t3,
            0.5 * t + 2 * t2 - 1.5 * t3, -0.5 * t2 + 0.5 * t3)


def _sample_cubic(padded, x, y, shape, xp):
    """Zero-based cubic convolution with zero extrapolation (not a B-spline)."""
    valid = (x >= 0) & (x <= shape[1] - 1) & (y >= 0) & (y <= shape[0] - 1)
    result = xp.zeros(x.shape, dtype=xp.float64)
    # Evaluate only the rays inside the image. The upper endpoint uses the
    # preceding interval, ensuring all four samples fit the padded array.
    xv, yv = x[valid].astype(xp.float64), y[valid].astype(xp.float64)
    ix = xp.minimum(xp.floor(xv).astype(xp.int64), shape[1] - 2)
    iy = xp.minimum(xp.floor(yv).astype(xp.int64), shape[0] - 2)
    wx, wy = _cubic_weights(xv - ix), _cubic_weights(yv - iy)
    values = xp.zeros(xv.shape, dtype=xp.float64)
    for row in range(4):
        for col in range(4):
            values += wy[row] * wx[col] * padded[iy + row, ix + col]
    result[valid] = values
    return result


def cubic_interpolate(image, x, y, *, backend="numpy"):
    """MATLAB interp2(V, x, y, 'cubic', 0), with ONE-based query coordinates."""
    xp = array_module(backend)
    image = xp.asarray(image, dtype=xp.float64)
    x, y = xp.broadcast_arrays(xp.asarray(x, dtype=xp.float64), xp.asarray(y, dtype=xp.float64))
    if image.ndim != 2 or min(image.shape) < 3:
        raise ValueError("cubic interpolation requires a 2D image of at least 3 x 3")
    if not bool(xp.isfinite(x).all() & xp.isfinite(y).all()):
        raise ValueError("interpolation coordinates must be finite")
    return _host(_sample_cubic(_pad_cubic(image, xp), x - 1, y - 1, image.shape, xp), xp)


def forward_project(image, detector: Detector, geometry=Geometry(), *, backend="numpy",
                    chunk_rays=64, view_indices=None, output=None, progress: Progress | None = None):
    """Integrate cubic samples along rays, in detector chunks to bound memory.

    Optional view_indices selects original view angles for validation/previews;
    it does not redistribute the angles. An output path streams float32 views
    directly to a MATLAB-compatible raw file. Return shape is (selected views,
    detectors), using a memmap when output is supplied.
    """
    _detector(detector)
    if not isinstance(chunk_rays, int) or isinstance(chunk_rays, bool) or chunk_rays < 1:
        raise ValueError("chunk_rays must be a positive integer")
    image = np.asarray(image, dtype=np.float64)
    if image.shape != (geometry.image_size,) * 2 or not np.isfinite(image).all():
        raise ValueError("image must be finite and match geometry.image_size")
    indices = np.arange(geometry.views) if view_indices is None else np.asarray(view_indices)
    if indices.ndim != 1 or not len(indices) or not np.issubdtype(indices.dtype, np.integer):
        raise ValueError("view_indices must be a nonempty sequence of integer indices")
    if np.any(indices < 0) or np.any(indices >= geometry.views):
        raise ValueError("view_indices are outside the configured angular sampling")
    xp = array_module(backend)
    shape = (len(indices), geometry.detectors)
    if output is None:
        projections = np.empty(shape, dtype=np.float32)
    else:
        projections = np.memmap(_output_path(output), dtype="<f4", mode="w+", shape=shape)
    # MATLAB forward geometry is single; fread produces a double image.
    single = np.float32
    sid, sdd, du = map(single, (geometry.sid, geometry.sdd, geometry.du))
    dx = single(geometry.dx)
    steps = int(np.ceil(single(sid + sdd) / single(dx / single(2))))
    # Mixed MATLAB arithmetic rounds ray coordinates to single precision.
    t = xp.asarray((np.arange(steps, dtype=np.float64) + 0.5) / steps, dtype=xp.float32)
    offsets = xp.asarray(np.arange(geometry.detectors) - geometry.detectors / 2 + 0.5, dtype=xp.float32)
    offsets *= single(du / sdd) if detector == "curved" else du
    padded = _pad_cubic(xp.asarray(image), xp)
    angles = geometry.angles
    for row, index in enumerate(indices):
        theta = angles[index]
        ct, st = single(np.cos(theta)), single(np.sin(theta))
        xs, ys = single(sid * ct), single(sid * st)
        for start in range(0, geometry.detectors, chunk_rays):
            stop = min(start + chunk_rays, geometry.detectors)
            offset = offsets[start:stop]
            if detector == "curved":
                angle = single(theta) + offset
                xd = xs - sdd * xp.cos(angle)
                yd = ys - sdd * xp.sin(angle)
            else:
                xd = xs - sdd * ct + offset * single(np.cos(theta + np.pi / 2))
                yd = ys - sdd * st + offset * single(np.sin(theta + np.pi / 2))
            dxr, dyr = xd - xs, yd - ys
            lengths = xp.sqrt(dxr * dxr + dyr * dyr)
            xt = xs + dxr[:, None] * t[None, :]
            yt = ys + dyr[:, None] * t[None, :]
            ix = xt / dx + single(geometry.image_size / 2 + 0.5)
            iy = yt / dx + single(geometry.image_size / 2 + 0.5)
            # MATLAB fread/interp2 keeps double sample values. The subsequent
            # multiplication by single ds rounds the double sum to single.
            samples = _sample_cubic(padded, ix - single(1), iy - single(1), image.shape, xp)
            values = xp.sum(samples, axis=1, dtype=xp.float64).astype(xp.float32) * (lengths / single(steps))
            projections[row, start:stop] = _host(values, xp)
        if progress is not None:
            progress(row + 1, len(indices))
    if isinstance(projections, np.memmap):
        projections.flush()
    return projections


def ramp_filter(projections, detector: Detector, geometry=Geometry(), *, backend="numpy"):
    """Reference detector weights and symmetric ramp convolution; CPU result."""
    _detector(detector)
    projections = np.asarray(projections, dtype=np.float32)
    if projections.shape != (geometry.views, geometry.detectors) or not np.isfinite(projections).all():
        raise ValueError("projections must be finite with shape (geometry.views, geometry.detectors)")
    xp = array_module(backend)
    count = geometry.detectors
    # MATLAB builds these axes in double, then moves them to single arrays.
    axis64 = (np.arange(count) - count / 2 + 0.5) * geometry.du
    spacing = float(np.float32(geometry.du) / np.float32(geometry.sdd)) if detector == "curved" else geometry.du
    coordinate = (np.arange(count, dtype=np.float32) - np.float32(count / 2 - 0.5)) * np.float32(spacing) if detector == "curved" else axis64
    coordinate = xp.asarray(coordinate, dtype=xp.float32)
    kernel = np.zeros(2 * count - 1, dtype=np.float32)
    k = np.arange(1, count, 2, dtype=np.float64)
    if detector == "curved":
        sine = np.sin(k.astype(np.float32) * np.float32(spacing))
        values = np.float32(-1 / np.pi ** 2) / (sine * sine)
    else:
        values = -1 / (k * np.pi * spacing) ** 2
    kernel[count - 1 - k.astype(int)] = values
    kernel[count - 1 + k.astype(int)] = values
    kernel[count - 1] = 1 / (2 * spacing) ** 2
    kernel = xp.asarray(kernel)
    if detector == "curved":
        weights = xp.cos(coordinate)
    else:
        sdd = np.float32(geometry.sdd)
        weights = sdd * sdd / xp.sqrt(sdd * sdd + coordinate * coordinate)
    data = xp.asarray(projections) * weights
    filtered = xp.empty_like(data)
    for view in range(geometry.views):
        full = xp.convolve(data[view], kernel, mode="full")
        filtered[view] = full[count - 1:2 * count - 1] * np.float32(spacing)
    return _host(filtered, xp)


def backproject(filtered, detector: Detector, geometry=Geometry(), *, backend="numpy",
                progress: Progress | None = None):
    """Port of the reference backprojection, including its boundary fallback."""
    _detector(detector)
    filtered = np.asarray(filtered, dtype=np.float32)
    if filtered.shape != (geometry.views, geometry.detectors) or not np.isfinite(filtered).all():
        raise ValueError("filtered must be finite with shape (geometry.views, geometry.detectors)")
    xp = array_module(backend)
    data = xp.asarray(filtered)
    axis = (np.arange(geometry.image_size) - geometry.image_size / 2 + 0.5) * geometry.dx
    x, y = xp.meshgrid(xp.asarray(axis, dtype=xp.float32), xp.asarray(axis, dtype=xp.float32))
    angles = xp.asarray(geometry.angles, dtype=xp.float32)
    spacing = float(np.float32(geometry.du) / np.float32(geometry.sdd)) if detector == "curved" else geometry.du
    if detector == "curved":
        coords = (np.arange(geometry.detectors, dtype=np.float32) - np.float32(geometry.detectors / 2 - 0.5)) * np.float32(spacing)
    else:
        coords = (np.arange(geometry.detectors) - geometry.detectors / 2 + 0.5) * spacing
    coords = xp.asarray(coords, dtype=xp.float32)
    img = xp.zeros((geometry.image_size,) * 2, dtype=xp.float32)
    sid, sdd = np.float32(geometry.sid), np.float32(geometry.sdd)
    for view in range(geometry.views):
        ct, st = xp.cos(angles[view]), xp.sin(angles[view])
        if detector == "curved":
            distance = xp.sqrt((x - sid * ct) ** 2 + (y - sid * st) ** 2)
            distance2 = distance * distance
            position = xp.arcsin((x * st - y * ct) / distance)
        else:
            distance = sid - x * ct - y * st
            distance2 = distance * distance
            position = sdd * (-x * st + y * ct) / distance
        bins = xp.floor((position - coords[0]) / np.float32(spacing)).astype(xp.int64)
        valid = (bins >= 0) & (bins < geometry.detectors - 1)
        bins = xp.where(valid, bins, 0)
        # Deliberately NOT divided by spacing; preserve the MATLAB formula.
        weight = xp.where(valid, position - coords[bins], np.float32(0))
        values = (1 - weight) * data[view, bins] + weight * data[view, bins + 1]
        img += sid / distance2 * values * np.float32(2 * np.pi / geometry.views)
        if progress is not None:
            progress(view + 1, geometry.views)
    return _host(img / np.float32(2), xp)


def reconstruct(projections, detector: Detector, geometry=Geometry(), *, backend="numpy",
                progress: Progress | None = None):
    """Weight, ramp-filter, then backproject a complete projection matrix."""
    filtered = ramp_filter(projections, detector, geometry, backend=backend)
    return backproject(filtered, detector, geometry, backend=backend, progress=progress)
