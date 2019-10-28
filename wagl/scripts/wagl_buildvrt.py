#!/usr/bin/env python

"""
A script build vrt files from the image datasets contained within
the HDF5 file produced by the wagl workflow.
We could make this more generic, and simply trawl the HDF5 file and group
imagery together based on their resolution and spatial extent. However, time
is short, so we'll simply use the DatasetName's and GroupName's defined within
wagl.constants to create several vrt files.

ADDITIONAL NOTES:
the gdalbuildvrt commands -tr and -te whilst being set, were overridden by
GDAL when it couldn't detect an extent (wierd though, as we're providing one).
The -a_srs was unaffected. However, as we'll be loading the vrt to write the
transform, we might as well write the crs as the same time.
Input filenames are assumed to be {granule_id}.wagl.h5, and the output
filenames will replace {wagl} with the base group name the datasets come from.
"""

import argparse
import subprocess
from pathlib import Path
import h5py
import rasterio

from wagl.constants import GroupName
from wagl.geobox import GriddedGeoBox

GROUPS = [
    GroupName.ELEVATION_GROUP,
    GroupName.EXITING_GROUP,
    GroupName.INCIDENT_GROUP,
    GroupName.INTERP_GROUP,
    GroupName.LON_LAT_GROUP,
    GroupName.REL_SLP_GROUP,
    GroupName.SAT_SOL_GROUP,
    # GroupName.SHADOW_GROUP,  # ignore; bool dtype is not supported by GDAL
    GroupName.SLP_ASP_GROUP,
    GroupName.STANDARD_GROUP,
]

# these two groups have an additional group layer before the actual datasets
OTHER = [
    GroupName.INTERP_GROUP,
    GroupName.STANDARD_GROUP,
]

PATH_FMT = 'HDF5:"{}":/{}'
FMT = '{}-{}'


def _buildvrt(out_fname, dataset_paths, band_names, nodata, crs, transform,
              verbose):
    """
    Build a vrt from wagl dataset paths and set the crs and transform.
    """
    # while gdalbuildvrt support "-999 -999 -999" for nodata for 3 bands,
    # it does not support "None None None" for a 3 band image
    # therfore if one of the bands has None for nodata, use a single instance
    if None in nodata:
        nodata = [None]

    cmd = [
        'gdalbuildvrt',
        '-separate',
        '-srcnodata',
        '{}'.format(' '.join(['{}'.format(i) for i in nodata])),
        '{}'.format(str(out_fname)),
    ]
    cmd.extend(dataset_paths)

    if verbose:
        print('executing command:')
        print(cmd)

    subprocess.check_call(cmd)

    with rasterio.open(out_fname, 'r+') as src:
        src.crs = crs
        src.write_transform(transform.to_gdal())

        # create band names
        for i, name in enumerate(band_names):
            src.set_band_description(i+1, name)


def _construct_out_filename(fname, group_name):
    """
    Construct a specifically formatted output filename.
    The vrt will be placed adjacent to the HDF5 file, as
    such write access is required.
    """
    basedir = fname.absolute().parent
    basename = fname.with_suffix('.vrt').name.replace(
        'wagl',
        group_name
    )

    out_fname = basedir.joinpath(Path(basename))

    return out_fname


def run(fname, verbose):
    fname = Path(fname)
    with h5py.File(str(fname), 'r') as fid:
        for grn_path, granule in fid.items():
            # resolution groups
            for grp_name in [g for g in granule if 'RES-GROUP-' in g]:
                res_group = granule[grp_name]
                rg_name = 'RG{}'.format(grp_name.split('-')[-1])

                # image groups
                for img_grp in GROUPS:
                    if img_grp == GroupName.INTERP_GROUP:
                        # interpolated atmospheric coefficients
                        dataset_paths = []
                        band_names = []
                        nodata = []
                        for gname, group in res_group[img_grp.value].items():
                            for dname, dataset in group.items():
                                vrt_path = PATH_FMT.format(
                                    fname.name,
                                    dataset.name
                                )
                                bname = FMT.format(gname, dname)
                                dataset_paths.append(vrt_path)
                                band_names.append(bname)
                                nodata.append(dataset.attrs.get('no_data_value'))

                        # output filename
                        out_fname = _construct_out_filename(
                            fname, 
                            FMT.format(rg_name, GroupName.INTERP_GROUP.value)
                        )

                        # all datasets from the interp group should have the
                        # same extent and crs
                        geobox = GriddedGeoBox.from_dataset(dataset)

                        # build/create vrt
                        _buildvrt(
                            out_fname,
                            dataset_paths,
                            band_names,
                            nodata,
                            geobox.crs.ExportToWkt(),
                            geobox.transform,
                            verbose
                        )
                    elif img_grp == GroupName.STANDARD_GROUP:
                        # reflectance, thermal ...
                        for _, grp1 in res_group[img_grp.value].items():
                            # lambertian, nbar, nbart
                            for gname, group in grp1.items():
                                dataset_paths = []
                                band_names = []
                                nodata = []
                                for dname, dataset in group.items():
                                    vrt_path = PATH_FMT.format(
                                        fname.name,
                                        dataset.name
                                    )
                                    bname = FMT.format(gname, dname)
                                    dataset_paths.append(vrt_path)
                                    band_names.append(bname)
                                    nodata.append(dataset.attrs.get('no_data_value'))

                                # output product group
                                out_fname = _construct_out_filename(
                                    fname,
                                    FMT.format(rg_name, gname)
                                )

                                # all datasets in this res-group and product
                                # group will have the same extent and crs
                                geobox = GriddedGeoBox.from_dataset(dataset)

                                # build/create vrt
                                _buildvrt(
                                    out_fname,
                                    dataset_paths,
                                    band_names,
                                    nodata,
                                    geobox.crs.ExportToWkt(),
                                    geobox.transform,
                                    verbose
                                )
                    else:
                        # everything else should have same heirarchy levels
                        dataset_paths = []
                        band_names = []
                        nodata = []
                        for name, obj in res_group[img_grp.value].items():
                            if not obj.attrs.get('CLASS') == 'IMAGE':
                                continue
                            vrt_path = PATH_FMT.format(
                                fname.name,
                                obj.name
                            )
                            band_names.append(name)
                            dataset_paths.append(vrt_path)
                            nodata.append(obj.attrs.get('no_data_value'))

                            # within a group, each IMAGE dataset will have the
                            # extent and crs
                            geobox = GriddedGeoBox.from_dataset(obj)

                        # output data group
                        out_fname = _construct_out_filename(
                            fname,
                            FMT.format(rg_name, img_grp.value)
                        )

                        # build/create vrt
                        _buildvrt(
                            out_fname,
                            dataset_paths,
                            band_names,
                            nodata,
                            geobox.crs.ExportToWkt(),
                            geobox.transform,
                            verbose
                        )


def _parser():
    """ Argument parser. """
    description = "Build a VRT from a HDF5 file output by wagl."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--filename", required=True,
                        help="The file from which to list the contents.")
    parser.add_argument("--verbose", action='store_true',
                        help="If set, then print each vrt command.")

    return parser


def main():
    """ Main execution. """
    parser = _parser()
    args = parser.parse_args()
    run(args.filename, args.verbose)
