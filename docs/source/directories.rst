ULA3 Directory Structure
========================

The processes implemented in ULA3 the following ancilliary data:

* Modtran
	* aerosol data,
	* ozone data, and
	* water vapour data.

* NBAR and Terrain Correction
	* elevation data
		* Both NBAR and Terrain correction use the DEM found in ``DIR_DEM``,
		* Terrain Correction also uses a special 3 second grid (prepared by David Jupp at CSIRO) for slope calculations (this has not been 			put in place yet).
	* ... what else?

Other ancilliary data is also required, but the author is not sure where this is used (hopefully someone else will also get into this crazy idea of documenting things and update it)

* Unknown uses
	* Earth-Sun distance,
	* Solar irradiance,
	* Land-Sea rasters (probably required for pixel quality).

The location of the data is partially defined in processor.conf, and partially (implicitly) within the code. There is 

Ancilliary data
---------------

The location of ancillary data is specified in the configuration file for the system, which lives in the file processor.conf within the package :py:mod:`ULA3.image_processor`. This is loaded and managed by :py:class:`ULA3.image_processor.ProcessManager`. The structure of the ancillary data was decided by members of the ULA team at some stage, in particular Leo Lymburner and Wenjun Wu. Note that this structure differs between systems. The following reflects the structure on VAYU at the time of writing.

Details on the file formats and structures within the following directories does not appear to be documented in general and one needs to infer this from the code that accesses the data. This has been noted where possible.

* ``ANCI_ROOT (= %(MAIN_ROOT)s/eoancillary)``: The root of the ancillary data. by default, all other ancillary data lives in directories under this. It is referenced heavily in other paths in the configuration file and directly in :py:mod:`ULA3.modtran`, where it is used to build a path to a dafault MODTRAN input file (*.tp5). At the time of writing, this was around line 608.

|
* ``EPHEM_DATA_ROOT (= %(ANCI_ROOT)s/ephemeris)``: Satellite ephemeris data (`Two Line Element Sets (TLEs) <http://en.wikipedia.org/wiki/Two-line_element_set>`_) for Landsat satellites. It is referred to in :py:func:`ULA3.image_processor.utils.calc_satellite_grids` which utilises :py:class:`ULA3.geodesic.Satellite` which was written by `Alex Ip <mailto:alex.ip@ga.gov.au>`_.

|
* ``DIR_Aerosol (= %(ANCI_ROOT)s/aerosol/AATSR/2.0)``: Aerosol data. This is referenced by :py:func:`image_processor.nbar.get_ancillary_data.get_aerosol_ancillary_data.process` which calls :py:func:`ULA3.ancillary.aerosol.get_aerosol_data`.

|
* ``DIR_BRDF (= %(ANCI_ROOT)s/BRDF/CSIRO_mosaic/MCD43A1.005)``: MODIS BRDF data. This is used when running MODTRAN to calculate various atmospheric paramters that feed into both :ref:`nbar-algorithm-label` and :ref:`tc-algorithm-label`.

|
* ``DIR_DEM (= %(ANCI_ROOT)s/elevation/world_1deg)``

|
* ``DIR_DEM_TC (= %(ANCI_ROOT)s/elevation/tc_aus_3sec)``: Contains a DSM of Australia which is used in the Terrain Correction algorithm (see :py:mod:`image_processor.nbar.radiative_transfer.run_tc`, which was written by `Simon Knapp <simon.knapp@ga.gov.au>`_). Note that while the name implies a 3 second DSM, it is apparently a 1 second DSM (i.e. 30m). In the process of running the Terrain correction algorithm, a region is clipped and resampled from this data using the GDAL program `gdal_warp <http://www.gdal.org/gdalwarp.html>`_. This implies that the DSM must be readable by GDAL and have an extent large enough to contain any region for which Terrain Correction is to be performed.

|
* ``DIR_EarthSun_LUT (= %(ANCI_ROOT)s/lookup_tables/earthsun_distance)``

|
* ``DIR_Ozone_LUT (= %(ANCI_ROOT)s/lookup_tables/ozone)``

|
* ``DIR_SatFilter (= %(ANCI_ROOT)s/lookup_tables/satellite_filter)``

|
* ``DIR_SolarIrradianceLUT (= %(ANCI_ROOT)s/lookup_tables/solar_irradiance)``

|
* ``DIR_WaterVapour (= %(ANCI_ROOT)s/water_vapour)``

|
* ``DIR_LandSea (= /short/v10/tmp/Land_Sea_Rasters)``
