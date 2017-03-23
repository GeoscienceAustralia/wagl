Acquisitions
============

The core class for the gaip package is the *Acquisition* object. An *Acquistion* contains much of the required information about the spectral band,
but also the satellite, sensor, datetime, gain, bias, spectral reponse etc.

An *Acquisition* instance is created for every band contained within a given group, such as a particular resolution.
The *Acquisition* contains geospatial context, and has direct read access to the spectral band that the *Acquisition* refers to.

From a top down perspective, a scene can contain zero or more granules, and each granule contains at least 1 *Acquisition*. The *AcquisitionsContainer*
is the class object containing all of this together, and can be thought of as similar to a HDF5 file which itself is structured similar to a UNIX directory.

A scene can be thought of as the root (/) level, and if the scene is tiled, i.e. made up of multiple granules/tiles, then these *Granules* would be the
next level.  Each *Granule* can contain multiple *Acquisitions*, each *Acquisition* instance is a direct reference to a given spectral band available
for that resolution group.

The following picture describes the rough layout of an *AcquisitionsContainer*.

.. image:: acq.png

The *AcquisitionsContainer* will detail whether or not it tiled, in which case how many *Granules* make up the entire geospatial footprint of the scene,
and how many *Groups* are contained within each *Granule*.
