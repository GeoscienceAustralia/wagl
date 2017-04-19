gaip 5.0 release notes
======================

* Version 5.0 has extra satellite/sensor functionality with the inclusion of Sentinel-2a.
* HDF5 is used as the backend file format which has resulted in a simplified workflow, reduction in the number of files output during a process run, consolodated the file I/O, and various other advantages.
* Surface Brightness Temperature (SBT) is now available as an accompanying product to NBAR and PQ.
* Extra layers of hierarchy have been added to the return of what originally was a single list of *Acquisitions*. An *AcquisitionsContainer* is now returned instead, as it was necessary in order to handle scenes comprising of multiple *Granules*, and multiple *Resolutions*.
* Additional bilinear interpolation functions. The intent is to deprecate the FORTRAN version.
* Command line utilities:
    * gaip_convert  *An unpacking/converting utility that converts HDF5 Tables to CSV, HDF5 images to GeoTiff, and metadata to yaml files.*
    * test_calculate_angles *Compares and evaluates each dataset contained within a **satellite-solar.h5** file against the same datasets contained within another file.*
    * test_dsm *Compare and evaluate two **dsm-extract.h5** files.*
    * test_exiting_angles *Compare and evaluate two **exiting-angles.h5** files.*
    * test_incident_angles *Compare and evaluate two **incident-angles.h5** files.*
    * test_relative_slope *Compare and evaluate two **relative-slope.h5** files.*
    * test_shadow_masks *Compare and evaluate two **shadow-masks** files.*
    * test_slope_aspect *Compare and evaluate two **slope-aspect.h5** files.*


Python compatibility
--------------------
Support for Python 2.7* has been dropped, and gaip now requires Python 3.5 or later.


gaip 4.0 release notes
======================

Version 4.0 has seen a major undertaking in regards to streamlining image I/O through 3rd party libraies, and for the better part, a single point of source.
The directory structure has been overhauled in order to simplify the code, not only for the users, but also for the developers.
A new workflow using Luigi has been introduced, and will simplify the job workflow.


Python compatibility
--------------------
gaip requires Python 2.7.6 or later.  Support for Python 3 is not yet complete.


What's new in gaip 4.0
----------------------

Most image based I/O is handled through rasterio, which itself is a wrapper for the GDAL library.  A new algorithm to calculate the solar and satellite angles was introduced, but unlike the previous Fortran code, this was built as a Python module via F2Py. This practice will continue and eventually replace the existing Fortran 77 code.  Other Fortran modules were replaced with existing Python libraries such as the gaussian filter.
Unsed C and Fortran 77 files were removed, with the remaining files collected into a single `src` directory.
An almost comlete unittesting framework has been inroduced, with features such as allowing users to test individual modules using pre-existing outputs. Thereby allowing the tests to pick up outputs generated at previous steps.
Almost all intermediate calculations are saved to disk.  This has allowed for a more rapid QA/QC to take place without infringing on quality or program run times.

* SceneDataset has been replaced by Acquisitions, which brings forth a simpler and more maintainable code base, whilst
* GriddedGeoBox brings a spatial association to a NumPy array.  Supplying features such as x & y array shape, CRS, co-ordinate transforms between the image domain and real world domain and vice vers, but also to different CRS's.  Bounding box extents for the array in not only the current frame of reference, but also to WGS84 Lon/Lat co-ordinates.
* Luigi has replaced the existing workflow, and presents a much friendlier understanding as to how the reflectance, pixel quality and fractional cover algorithms are run.
* The NBAR algorithm has been updated to yield a more accurate method of calculating surface reflectance.  Outputs from the reflectance algorithm includes Lambertian, NBAR, and terrain corrected NBAR.
* Individual workflows are now defined by nbar.py fc.py and pq.py.
* Makefiles have been simplified, with most of the F2Py Fortran modules now built through setup.py which is a more up-to-date method of installing Python based projects.
* Fortran files and modules have been renamed which allows an easier and more recogniseable way understanding the files.
* Large workflows have been broken up into smaller individual workflows, thus allowing for easier maintainence, but has also had the added benefit of memory reduction.
* Most string formatting uses the .format() method rather than C-style % insertion which presents a more readable raw form.
* Code base has been re-vamped to be PEP8 compliant.
* Ancillary data extraction has been simplified.
