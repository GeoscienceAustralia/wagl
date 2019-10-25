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
from typing import Optional, Tuple, Iterable, Dict
from pathlib import Path
import shutil
import subprocess
import tempfile
import enum
from posixpath import join as ppjoin
import logging
import numpy as np
import h5py
from wagl.hdf5 import write_dataframe, write_scalar
from wagl.constants import (
    ALBEDO_FMT,
    Albedos,
    BandType,
    DatasetName,
    GroupName,
    POINT_ALBEDO_FMT,
    POINT_FMT,
)
from wagl.hdf5.compression import H5CompressionFilter
from wagl.hdf5 import attach_attributes
from wagl.modtran_profiles import MIDLAT_SUMMER_ALBEDO, TROPICAL_ALBEDO

MAX_FILTER_SIZE = 101
_LARGE_FILTER_SIZE = 40
_HSTEP = 10.0
_TP5_FMT = pjoin(POINT_FMT, ALBEDO_FMT, "".join([POINT_ALBEDO_FMT, ".tp5"]))

_LOG = logging.getLogger(__name__)


def read_tp7(tp7_file: Path) -> dict:
    """Returns a Point Spread Function (PSF) from a '*.tp7' file.

    :param tp7_file: A full path '*.tp7' file from MODTRAN 5.4 output.
    """
    psf_dict = dict()
    with open(tp7_file, "r") as fid:
        lines = fid.readlines()
        for idx, line in enumerate(lines):
            if line.startswith("BAND-"):
                psf_dict[line.strip()] = np.asarray(
                    [float(val) for val in lines[idx + 1].split()]
                )

    return psf_dict


def compute_fwhm(psf_data: np.ndarray, prange: np.ndarray, hstep: float) -> float:
    """Returns a Full Width at Half Maximum (FWHM) computed from PSF data.

    :param psf_data: A numpy array containing the psf values for a specific band.
    :param prange: A numpy array containing integer increment of hstep.
    :param hstep: A MODTRAN 5.4 step size used in computing the psf data (default=10.0).
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
    """Returns a maximum filter filter size.

    :param xres: pixel size in x-direction.
    :param yres: pixel size in y-direction.
    :param nlarge: maximum allowed filter size.
    """

    return 2 * np.int(nlarge * 25.0 / max(xres, yres)) + 1


def compute_filter_matrix(
    psf_data: np.ndarray, xres: float, yres: float, nlarge: int, hstep: float
) -> np.ndarray:
    """Returns a adjacency kernel derived from Point Spread Function.

    compute the adjacency 2D filter matrix based on the 1D FWHM of the PSF as the radius
    of the matrix. If the pixel size is [xres, yres] for the x and y direction, then the
    number of pixels in the matrix away from the center is computed (refer to ATBD for water
    atmospheric correction by David and Fuqin.

    :param psf_data: Point Spread Function data for a spectral band.
    :param xres: Pixel size in x-direction.
    :param yres: Pixel size in y-direction.
    :param nlarge: A maximum filter size set by the user.
    :param hstep: A MODTRAN step size used in computing the psf data.
    """

    # compute prange
    prange = np.arange(len(psf_data)) * hstep

    def _interp(_range):
        """Returns a linear interpolation of psf_data at _range"""
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
    max_filter_size = min(_max_filter_size(xres, yres, nlarge), MAX_FILTER_SIZE)

    # compute fwhm for spectral psf data
    fwhm = compute_fwhm(psf_data, prange, hstep)

    # get the number of pixel from the center of a pixel in x and y direction
    num_pixel_x = np.int(fwhm / xres)
    num_pixel_y = np.int(fwhm / yres)

    # check if number of pixels in a filter is greater than maximum pixel allowed in a filter
    # and reset to a maximum allowed
    if 2 * num_pixel_x + 1 > max_filter_size:
        num_pixel_x = np.int((max_filter_size - 1) / 2)
    if 2 * num_pixel_y + 1 > max_filter_size:
        num_pixel_y = np.int((max_filter_size - 1) / 2)

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


def prepare_modtran54(
    acquisition: object, coordinate: int, albedo: enum, basedir: Path, modtran_exe: Path
) -> None:
    """Prepares the working directory for a MODTRAN 5.4 execution.

    :param acquisition: 'acquisition' object.
    :param coordinate: A point to execute MODTRAN run at.
    :param albedo: A 'instance' of enum 'Albedos'.
    :param basedir: A full path to a base directory to setup MODTRAN execution.
    :param modtran_exe: A full path to a MODTRAN 5.4 executable.
    """

    data_dir = modtran_exe.parent.joinpath("DATA")
    if not data_dir.exists():
        raise OSError("Cannot find MODTRAN 5.4")

    point_dir = basedir.joinpath(POINT_FMT.format(p=coordinate))
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
        str(acquisition.spectral_filter_filepath),
        str(modtran_work.joinpath(acquisition.spectral_filter_name)),
    )


def format_tp5(json_data: dict) -> Tuple[dict, dict]:
    """A utility to format MODTRAN 5.4 '*.tp5' input.

    :param json_data: json data containing the MODTRAN 5.4 input parameters.
    :return:
        A 'str' formatted tp5 data at Albedo (0) and dict with input parameters.
    """

    input_data = json_data["MODTRAN"][0]["MODTRANINPUT"]
    kwargs = {
        "albedo": float(Albedos.ALBEDO_0.value),
        "water": input_data["ATMOSPHERE"]["H2OSTR"],
        "ozone": input_data["ATMOSPHERE"]["O3STR"],
        "filter_function": input_data["SPECTRAL"]["FILTNM"],
        "aerosol_type": json_data["aerosol_type"],
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
    return data, kwargs


def create_dataset(
    group: h5py.Group,
    band_name: str,
    shape: Iterable[int],
    attrs: Dict,
    dtype=np.float32,
    chunks: Tuple = (),
    compression: H5CompressionFilter = H5CompressionFilter.LZF,
    filter_opts: Optional[Dict] = None,
):
    """ creates dataset and attaches attributes for h5 object. """

    if filter_opts is None:
        filter_opts = {}
    else:
        filter_opts = filter_opts.copy()

    if "chunks" not in filter_opts:
        filter_opts["chunks"] = True if not chunks else chunks

    kwargs = compression.config(**filter_opts).dataset_compression_kwargs()
    ds = group.create_dataset(band_name, shape=shape, dtype=dtype, **kwargs)
    attach_attributes(ds, attrs)

    return ds


def compute_adjacency_filter(
    container: object,
    granule: str,
    json_data: dict,
    nvertices: int,
    modtran54_exe: str,
    out_group: object,
    aerosol_type: str,
) -> None:
    """Computes adjacency filter from a MODTRAN 5.4 output.

    This method generates a '*.tp5' input to a MODTRAN 5.4 program
    at a central point in nvertices. MODTRAN 5.4 is executed to compute
    a Point Spread Function. The MODTRAN 5.4 output '*.tp7' is read and
    used in computing adjacency kernel for all the visible and NIR bands
    available in acquisitions inside a container. The intermediate results,
    per band PSF data are written to a h5py object along with associated
    '*.tp5' input data. The adjacency filter (not normalized) are also
    written to h5py object for bands in PSF data.

    :param container: An instance with acquisitions class
    :param granule:  A name of a granule.
    :param json_data: MODTRAN 6.0 inputs from all points in vertices.
    :param nvertices: The total number of points in vertices (nvertices).
    :param modtran54_exe: A  full path to a MODTRAN 5.4 executable.
    :param out_group: A `File` object from which to write the dataset to.
    :param aerosol_type: An 'instance' of AerosolModel to configure MODTRAN *.tp5 input.
    """

    # get list of all the acquistions within a container
    acqs = sum(
        [
            container.get_acquisitions(
                granule=granule, group=grp_name, only_supported_bands=False
            )
            for grp_name in container.supported_groups
        ],
        [],
    )

    supported_band_acqs = [acq for acq in acqs if acq.supported_band]
    reflective_acqs = [
        acq for acq in supported_band_acqs if acq.band_type == BandType.REFLECTIVE
    ]

    center_point = int(np.floor(nvertices / 2))

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir = Path(tmp_dir)

        # subset the json data to select only center point data
        json_data = json_data[(center_point, Albedos.ALBEDO_0)]
        json_data["aerosol_type"] = aerosol_type.value

        # generate tp5 data needed for MODTRAN 5.4 execution
        tp5_data, input_data = format_tp5(json_data)

        # write tp5 data to h5 file object
        dname = ppjoin(
            POINT_FMT.format(p=center_point),
            ALBEDO_FMT.format(a=Albedos.ALBEDO_0.value),
            DatasetName.TP5.value,
        )

        if out_group is None:
            out_group = h5py.File("atmospheric-inputs.h5", "w")

        if GroupName.ATMOSPHERIC_INPUTS_GRP.value not in out_group:
            out_group.create_group(GroupName.ATMOSPHERIC_INPUTS_GRP.value)

        group_name = out_group[GroupName.ATMOSPHERIC_INPUTS_GRP.value]
        iso_time = supported_band_acqs[0].acquisition_datetime.isoformat()
        group_name.attrs["acquisition-datetime"]: iso_time
        write_scalar(np.string_(tp5_data), dname, group_name, input_data)

        # prepare directory and data symlinks for MODTRAN 5.4 execution
        prepare_modtran54(
            reflective_acqs[0],
            center_point,
            Albedos.ALBEDO_0,
            tmp_dir,
            Path(modtran54_exe),
        )
        tp5_fname = tmp_dir.joinpath(
            _TP5_FMT.format(p=center_point, a=Albedos.ALBEDO_0.value)
        )
        with open(tp5_fname, "w") as src:
            src.writelines(tp5_data)

        work_path = tmp_dir.joinpath(
            POINT_FMT.format(p=center_point),
            ALBEDO_FMT.format(a=Albedos.ALBEDO_0.value),
        )

        # execute MODTRAN 5.4
        subprocess.check_call([modtran54_exe], cwd=str(work_path))

        # read tp7 output file from MODTRAN 5.4
        tp7_file = list(work_path.glob("*.tp7"))[0]
        psf_data = read_tp7(tp7_file)

        # determine the output group/file
        if out_group is None:
            out_group = h5py.File(
                "atmospheric-results.h5", driver="core", backing_store=False
            )

        if GroupName.ATMOSPHERIC_RESULTS_GRP.value not in out_group:
            out_group.create_group(GroupName.ATMOSPHERIC_RESULTS_GRP.value)

        group_name = out_group[GroupName.ATMOSPHERIC_RESULTS_GRP.value]

        for band, _psf_data in psf_data.items():

            # only process the support band acquisitions
            if band not in [acq.band_name for acq in supported_band_acqs]:
                continue

            # write psf data into the h5 object
            attrs = {
                "band_name": band,
                "description": f"{band} Point Spread Function data from MODTRAN .tp7 output",
            }

            dname_psf = ppjoin(DatasetName.PSF.value, band)
            psf_ds = create_dataset(group_name, dname_psf, _psf_data.shape, attrs)
            psf_ds[:] = _psf_data.astype('float32')

            acq = [acq for acq in supported_band_acqs if acq.band_name == band][0]
            xres, yres = acq.resolution
            filter_matrix = compute_filter_matrix(
                _psf_data, xres, yres, _LARGE_FILTER_SIZE, _HSTEP
            )
            # check if filter_matrix is symmetric
            if not np.allclose(
                filter_matrix,
                filter_matrix[tuple([slice(None, None, -1)] * filter_matrix.ndim)],
            ):
                _LOG.warning(f"Adjacency kernel is not symmetric for {band}")

            attrs["description"] = f"{band} Adjacency filter derived from Point Spread Function"
            dname_filter = ppjoin(DatasetName.ADJACENCY_FILTER.value, acq.band_name)

            filter_ds = create_dataset(group_name, dname_filter, filter_matrix.shape, attrs)
            filter_ds[:] = filter_matrix.astype('float32')
