Acquisitions
============

The core class for the wagl package is the *Acquisition* object. An *Acquistion* contains much of the required information about the spectral band,
but also the satellite, sensor, datetime, gain, bias, spectral reponse etc.

An *Acquisition* instance is created for every band contained within a given group, such as a particular resolution.
The *Acquisition* contains geospatial context, and has direct read access to the spectral band that the *Acquisition* refers to.

From a top down perspective, a scene can contain zero or more granules, and each granule contains at least 1 *Acquisition*. The *AcquisitionsContainer*
is the class object containing all of this together, and can be thought of as similar to a HDF5 file which itself is structured similar to a UNIX directory.
This HDF5 stylisation is part of the design choice in moving to using HDF5 as the main backend for the algorithmic workflow, in that similar concepts and model designs can be
shared throughout the codebase.

A scene can be thought of as the root (/) level, and if the scene is tiled, i.e. made up of multiple granules/tiles, then these *Granules* would be the
next level.  Each *Granule* can contain multiple *Acquisitions*, each *Acquisition* instance is a direct reference to a given spectral band available
for that resolution group.

The following picture describes the rough layout of an *AcquisitionsContainer*.

.. image:: diagrams/acquisitions_container.png

The *AcquisitionsContainer* will detail whether or not it tiled, in which case how many *Granules* make up the entire geospatial footprint of the scene,
and how many *Groups* are contained within each *Granule*.

* *AcquisitionsContainer.tiled* Returns whether or not the scene is tiled and comprises of multiple *Granules*.
* *AcquisitionContainer.granules* Returns the *Granule* names associated with the scene. Returns *[None]* if there are no *Granules*.
* *AcquisitionContainer.groups* Returns the *Group* names associated with the scene. **Note!**:

       **A *Group* name is required when defining a new satellite/sensor acquisitions read routine.**

* *AcquisitionsContainer.get_acquisitions* Given a *Granule* and *Group*, return a list of *Acquisition* objects.
* *AcquisitionsContainer.get_granule* Returns a *dictionary* containing lists of *Acquisition* objects for each *Group* associated with the given *Granule*.
* *AcquisitionsContainer.get_root* Return a *str* to the root level for a given *Granule* and *Group*. Example:

       .. code-block:: python

          >>> scene = acquisitions('S2A_USER_PRD_MSIL2A_PDMC_20160120T071902_R016_V20160120T003331_20160120T003331.SAFE')
          >>> print(scene.get_root('my/work/directory', granule='S2A_USER_MSI_L2A_TL_SGS__20160120T053143_A003016_T55KBQ_N02.01', group='R10m')

An *Acquisition* instance also defines the tiling/chunking logic of any given image processing routine, as well as how the image result is stored on disk.
This is to provide a more dynamic and flexible processing capability for a variety of sensors that become supported by wagl.
The property `tile_size` describes the underlying storage tiles for the acquisition, and the method `tiles` returns a tile generator that can be looped over
and return a tile/chunk of data for processing.
