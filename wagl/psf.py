#! /usr/bin/env python3

"""
Adjacency filter derivation utilities
-------------------------------------

The marine atmospheric correction requires correction of spatial
pixel mixing of radiance from nearby(adjacent pixels) caused by
the atmospheric scattering. This program is written to derive
the adjacency filter to account for the pixel mixing due to nearby
pixels and correct for its affects. Initial program was written
by Dr. Li and Dr. Jubb in Fortran. The derivation of the kernel
is based on the ATBD written by Dr. David Jubb and Fuqin Li 2018
for Geosciences Australia's Analysis Ready Data program.

"""

from os.path import join as pjoin
from typing import Optional, Union
from pathlib import Path
import shutil
import subprocess
import tempfile
import enum

import numpy as np
from wagl.constants import ALBEDO_FMT, Albedos, BandType, POINT_ALBEDO_FMT, POINT_FMT
from wagl.modtran_profiles import MIDLAT_SUMMER_ALBEDO, TROPICAL_ALBEDO


def read_tp7(tp7_file: Path, nbands: Union[int, None]) -> dict:
    """
    program to read *.tp7 to extract point spread function.
    Point Spread Function is *.tp7 file are stored at the end
    of the large text file beginning above the line with str '-9999'.
    The PSF written per row (till nbands) in a text files with
    largest band dataset at the row above the line with '-9999' index.

    :param tp7_file: a 'Path' to a .tp7 file from MODTRAN output
    :param nbands: total number of bands used in computing Point Spread
                Function in a MODTRAN program. It is assumed that *.tp7
                has its PSF data (per row to a band) written with high
                numbered bands from the bottom of the file.
    :return
        A dict containing the band number as a key and psf value as a value
    """
    with open(tp7_file, "r") as fid:
        lines = fid.readlines()[::-1]
        if nbands is None:
            for idx, line in enumerate(lines):
                if "Bands=" in line:
                    nbands = int(line.split("=")[1])
        print(nbands)
        for idx, line in enumerate(lines):
            if "-9999." in line:
                return {
                    "Band_{}".format(nbands - _band): np.asarray(
                        [float(val) for val in lines[_band + idx + 1].split()]
                    )
                    for _band in range(nbands)
                }


def compute_fwhm(psf_data: np.ndarray, prange: np.ndarray, hstep: float) -> float:
    """
    calculate Full Width at Half Maximum (FWHM) values for psf data.

    :param psf_data: 1 dimensional numpy array containing the psf values for a band
    :param prange: 1 dimentional numpy array containing integer increment of hstep
    :param hstep: a MODTRAN step size used in computing the psf data
                   a default value is 10.0, unless new step size is configured in
                   future MODTRAN run to compute psf, hstep is a constant value of 10.0.
    :return
        A 'float' value for FWHM of psf data for a spectral band
    """
    accum = 0.0
    bigp = 0.0
    pival = 4.0 * np.arctan(1.0)

    for idx, val in enumerate(psf_data):
        accum = accum + val * prange[idx]
        if val > bigp:
            bigp = val

    return np.sqrt(2.0 * pival * hstep * accum / bigp)


def _max_filter_size(xres: float, yres: float, nlarge: int) -> int:
    """
    compute the maximum filter filter size

    :param xres: pixel size in x-direction
    :param yres: pixel size in y-direction
    :param nlarge: maximum allowed filter size

    :return
        A maximum filter size
    """
    return 2 * np.int(nlarge * 25.0 / max(xres, yres)) + 1


def compute_filter_matrix(
    psf_data: np.ndarray, xres: int, yres: int, nlarge: int, hstep: float
) -> np.ndarray:
    """
    compute the adjacency 2D filter matrix based on the 1D FWHM of the PSF as the radius
    of the matrix. If the pixel size is [xres, yres] for the x and y direction, then the
    number of pixels in the matrix away from the center is computed (refer to ATBD for water
    atmospheric correction by David and Fuqin.

    :param psf_data: Point Spread Function data for a spectral band
    :param xres: pixel size in x-direction
    :param yres: pixel size in y-direction
    :param nlarge: maximum filter size set by the user
    :param hstep: a MODTRAN step size used in computing the psf data

    :return:
        A np.ndarray containing adjacency kernel
    """

    # compute prange
    prange = np.arange(len(psf_data)) * hstep

    def _interp(_range):
        """ compute psf as _range by linear interpolation of psf_data"""
        if _range < 0.0 or _range > max(prange):
            return -1.0

        idx = np.int(_range / hstep)
        if idx == len(psf_data) - 1:
            return psf_data[idx]
        return (
            psf_data[idx] * (prange[idx + 1] - _range)
            + psf_data[idx + 1] * (_range - prange[idx])
        ) / hstep

    # compute max filter size
    max_filter_size = _max_filter_size(xres, yres, nlarge)

    # compute fwhm for spectral psf data
    fwhm = compute_fwhm(psf_data, prange, hstep)

    # get the number of pixel from the center of a pixel in x and y direction
    num_pixel_x = np.int(fwhm / xres)
    num_pixel_y = np.int(fwhm / yres)

    # check if number of pixels in a filter is greater than maximum pixel allowed in a filter
    # and reset to a maximum allowed
    if 2 * num_pixel_x + 1 > max_filter_size:
        num_pixel_x = (max_filter_size - 1) / 2
    if 2 * num_pixel_y + 1 > max_filter_size:
        num_pixel_y = (max_filter_size - 1) / 2

    # check if psf needs interpolation or not
    interp_psf = False
    if np.sqrt(xres * yres) <= 1.5 * hstep:
        interp_psf = True

    adj_filter = np.zeros((2 * num_pixel_y + 1, 2 * num_pixel_x + 1), dtype=float)

    # compute the kernel matrix
    for y in range(num_pixel_y * 2 + 1):
        ycoord = float(y - num_pixel_y) * yres
        for x in range(num_pixel_x * 2 + 1):
            xcoord = float(x - num_pixel_x) * xres
            val = np.sqrt(xcoord ** 2 + ycoord ** 2)
            if val > fwhm:
                continue
            if interp_psf:
                adj_filter[y, x] = _interp(val)
                continue

            val1 = np.sqrt((xcoord + xres / 2.0) ** 2 + (ycoord + yres / 2.0) ** 2)
            val2 = np.sqrt((xcoord + xres / 2.0) ** 2 + (ycoord - yres / 2.0) ** 2)
            val3 = np.sqrt((xcoord - xres / 2.0) ** 2 + (ycoord + yres / 2.0) ** 2)
            val4 = np.sqrt((xcoord - xres / 2.0) ** 2 + (ycoord - yres / 2.0) ** 2)
            adj_filter[y, x] = (
                _interp(val)
                + (_interp(val1) + _interp(val2) + _interp(val3) + _interp(val4)) / 4.0
            ) / 2.0
    return adj_filter


def get_adjacency_kernel(
    psf_data_dict: dict,
    xres: float,
    yres: float,
    max_filter_size: int,
    hstep: Optional[float] = 10.0,
) -> dict:
    """
    Reads a PSF from a *.tp7 MODTRAN 5.4 output and computes adjacency filter.

    :param psf_data_dict: A dict with Band Number as a key and PSF data as a value
    :param xres: pixel size in x-direction
    :param yres: pixel size in y-direction
    :param max_filter_size: maximum filter size set by the user
    :param hstep: a MODTRAN step size used in computing the psf data
               a default value is 10.0, unless new step size is configured in
               future MODTRAN run to compute psf, hstep is a constant value of 10.0.
    return:
        A dict with acquisition band as a key and kernel (np.ndarray) as a value
    """

    return {
        band: compute_filter_matrix(_psf_data, xres, yres, max_filter_size, hstep)
        for band, _psf_data in psf_data_dict.items()
    }


def prepare_modtran54(
    acquisitions: list, coordinate: int, albedo: enum, basedir: Path, modtran_exe: Path
) -> None:
    """
    Prepares the working directory for a MODTRAN 5.4 execution.

    :param acquisitions: An 'instance' of list of 'acquisition'
    :param coordinate: An 'int' indicating the position of MODTRAN run
    :param albedo: A 'instance' of enum 'Albedos'
    :param basedir: A 'Path', base directory to setup MODTRAN execution
    :param modtran_exe: A 'Path' to MODTRAN 5.4 executable
    """

    data_dir = modtran_exe.parent.joinpath("DATA")
    if not data_dir.exists():
        raise OSError("Cannot find MODTRAN 5.4")

    point_dir = basedir.joinpath(POINT_FMT.format(p=coordinate))
    acq = [acq for acq in acquisitions if acq.band_type == BandType.REFLECTIVE][0]
    modtran_work = point_dir.joinpath(ALBEDO_FMT.format(a=albedo.value))

    if not modtran_work.exists():
        modtran_work.mkdir(parents=True)

    with open(str(modtran_work.joinpath("mod5root.in")), "w") as src:
        src.write(POINT_ALBEDO_FMT.format(p=coordinate, a=albedo.value) + "\n")

    # create a symbolic link to a MODTRAN data directory
    symlink_dir = modtran_work.joinpath("DATA")
    if symlink_dir.exists():
        symlink_dir.unlink()
    symlink_dir.symlink_to(data_dir)

    shutil.copy(
        str(acq.spectral_filter_filepath),
        str(modtran_work.joinpath(acq.spectral_filter_name)),
    )


def format_tp5(json_data: dict) -> str:
    """
    Creates a string formatted tp5 file for albedo (0)
    using the input_data from json dat
    :param json_data: json data containing the MODTRAN input parameters
    :return:
        A 'str' formatted tp5 data
    """
    input_data = json_data["MODTRAN"][0]["MODTRANINPUT"]
    kwargs = {
        "albedo": float(Albedos.ALBEDO_0.value),
        "water": input_data["ATMOSPHERE"]["H2OSTR"],
        "ozone": input_data["ATMOSPHERE"]["O3STR"],
        "filter_function": input_data["SPECTRAL"]["FILTNM"],
        "visibility": input_data["AEROSOLS"]["VIS"],
        "elevation": input_data["GEOMETRY"]["H2ALT"],
        "sat_height": input_data["GEOMETRY"]["H1ALT"],
        "sat_view": input_data["GEOMETRY"]["OBSZEN"],
        "sat_azimuth": input_data["GEOMETRY"]["TRUEAZ"],
        "doy": input_data["GEOMETRY"]["IDAY"],
        "lat": input_data["GEOMETRY"]["PARM1"],
        "lon": input_data["GEOMETRY"]["PARM2"],
        "time": input_data["GEOMETRY"]["GMTIME"],
    }
    atm_model = input_data["ATMOSPHERE"]["MODEL"]
    if atm_model == "ATM_TROPICAL":
        data = TROPICAL_ALBEDO.format(**kwargs)
    elif atm_model == "ATM_MIDLAT_SUMMER":
        data = MIDLAT_SUMMER_ALBEDO.format(**kwargs)
    else:
        raise ValueError(f"{atm_model} is not recognized atmosphere model")
    return data


def run_modtran54(
    acqs: list,
    json_data: dict,
    nvertices: int,
    modtran54_exe: str,
    num_bands: Optional[int] = None,
) -> None:
    """
    Run MODTRAN 5.4 and extract Point Spread Function data from *.tp7 file
    :param acqs: An 'instance' of list of 'acquisition'
    :param json_data: A 'dict' containing all the MODTRAN 6.0 inputs for nvertices
    :param nvertices: An 'int' indicating total number of vertices
    :param modtran54_exe: A 'str' MODTRAN 5.4 executable path
    :param num_bands: An 'int' indicating total number of spectral bands in PSD data
    :return:
        A 'dict' containing the band number as a key and psf value as a value
    """
    tp5_fmt = pjoin(POINT_FMT, ALBEDO_FMT, "".join([POINT_ALBEDO_FMT, ".tp5"]))
    center_point = np.floor(nvertices / 2)

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir = Path(tmp_dir)
        tp5_data = format_tp5(json_data[(center_point, Albedos.ALBEDO_0)])
        prepare_modtran54(
            acqs, center_point, Albedos.ALBEDO_0, tmp_dir, Path(modtran54_exe)
        )
        tp5_fname = tmp_dir.joinpath(
            tp5_fmt.format(p=center_point, a=Albedos.ALBEDO_0.value)
        )
        with open(tp5_fname, "w") as src:
            src.writelines(tp5_data)
        work_path = tmp_dir.joinpath(
            POINT_FMT.format(p=center_point),
            ALBEDO_FMT.format(a=Albedos.ALBEDO_0.value),
        )
        subprocess.check_call([modtran54_exe], cwd=str(work_path))
        tp7_file = list(work_path.glob("*.tp7"))[0]

        return read_tp7(tp7_file, num_bands)
