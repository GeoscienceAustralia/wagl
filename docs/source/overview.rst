Overview
========

wagl is a Python package for standardising imagery captured in the visible, near infrared and thermal infrared spectral channels. The imagery is standarised by correcting for atmospheric conditions, as well as the satellite and solar geometry.

For the visible, near infrared and mid infrared channels, the imagery is standardised to surface reflectance for the following models:

    * lambertian
    * BRDF corrected (NBAR)
    * BRDF + Terrain Illumination correction (NBART)

For thermal channels, the imagery is currently only standardised to:

    * Surface Brightness Temperature (SBT)

However, emissivity and Land Surface Temperature will soon become available.

The aim of wagl is to allow anyone to be able to process satellite imagery to a standardised form allowing for comparative analysis across space and time, with minimal effort. The api allows for users who want finer grained control for testing individual components, as well the entire workflow.

The luigi task manager for the *wagl.multifile_workflow* allows the user to escape having to re-run the entire workflow from scratch each time. This allows a user to delete results for the task they're executing, then re-execute without having to re-calculate all the prior dependencies. This is extremely useful for testing and validation, and allows for rapid development.

The *wagl.singlefile_workflow*, allows for mainstream mass production whilst not overloading the filesystem with hundreds/thousands of files being produced **per scene**. It also serves for direct archiving of results, easier file distribution/sharing with coleagues and users, as well as a simplified structure for testing and evaluation between versions and/or archives of wagl outputs.

The api provides command line access for executing any/all tasks of the *wagl.multifile_workflow*, as well as a bunch of tools useful for inspecting/converting/testing/comparing the file outputs created by wagl. And as usual for a Python package, all wagl modules and functions can be called interactively from within Python, allowing users to compute, test and develop using the standard Python interpretor. However, the command line access is recommended and allows the user to not worry about filenames and file I/O. In fact the primary requirement for wagl was to have command line access.

HDF5 serves as the backing store for all wagl outputs, and serves a multi-purpose storage utility for wagl:

    * Utilising a single I/O api suitable for Imagery and Tables, both of which wagl frequently outputs
    * Less files
    * A variety of compression filters (fast and light compression, to slow and heavy compression)
    * Archiving capability
    * Easier dataset comparisons between versions
    * Allows the storing of comparison evaluations results alongside both inputs (test and reference), which is useful for knowing exactly what inputs were used for testing a given verion
    * Storing metadata alongside the datasets, i.e. descriptions of what the result represents (previous versions required in-depth knowledge of the system
    * Simplified workflow

Currently only MODTRAN is supported for evaluating the radiative transfer, but in time it is anticipated that 6S would also be supported.

Most imagery related functions internal to wagl are processed in a tiled/chunked fashion in order to process images larger than available memory. With the idea being that wagl could be deployed on any compute architecture (desktop, laptop, compute cluster, cloud, etc).
