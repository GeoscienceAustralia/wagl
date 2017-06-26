Atmospherics
============

Specifically dealing with evaluating the radiative transfer for a defined set of spectra. The *gaip* package utilises `MODTRAN <http://modtran.spectral.com/>`_ for evalutating the radiative transfer.

The default number of points to evaluate the radiative transfer at is 9 for *nbar*, & 25 for *sbt*. If the *standard* model is chosen, then 25 points will be evaluated for both *nbar* and *sbt*. This is to make the workflow simpler and the output directories can be shared by both models, and share the same dependencies. The resulting *lambertian*, *nbar* and *nbart* datasets maybe slightly more accurate when using 25 points than 9 points.

Currently, the number of points are defined to cover the extents of an entire *scene*, or *granule* if the scene is *tiled*. For example a Sentinel-2a scene containing 15 granules, will subsequently have 15 * n points to evaluate the radiative transfer.

For the *nbar* model, the labels used for each radiative transfer calculation are:

    * albedo-0
    * albedo-1
    * albedo-t

For the *sbt* model, the label used for radiative transfer calculation is:

    * albedo-th

The resulting coefficient factors for the *nbar* model are:

    * dir:

      * Direct irradiance at the surface

    * dif:

      * Diffuse irradiance at the surface

    * ts:

      * Direct transmittance in the solar direction

    * fv:

      * Direct fraction in the view direction; where :math:`fv = tv / (tv + tdv)`

        * tv:

          * Direct transmittance in the view direction

        * tdv:

          * Diffuse transmittance in the view direction

    * fs:

      * Direct fraction in the solar direction; where :math:`fs = ts / (ts + tds)`

        * ts:

          * Direct transmittance in the solar direction

        * tds:

          * Diffuse transmittance in the solar direction

    * b:

      * Path radiance due to atmospheric scattering

    * s:

      * Atmospheric albedo

    * a:

      * :math:`a = (dir + dif) / pi * (tv + tdv)`


The resulting coefficient factors for the *sbt* model are:

    * path_up
    * path_down
    * transmittance_up

For scenes that are *tiled* such as Sentinel-2a, the ancillary data used for input into the radiative transfer calculation, are averaged. This is to get around higher resolution datasets, as well as the low resolution ancillary.

For the *sbt* model, ancillary data is gathered at each point. This is different to the *nbar* model which (depending on the ancillary) is gathered at a single location, such as the scene centre, and used for input into each point location for the radiative transfer calculations.
