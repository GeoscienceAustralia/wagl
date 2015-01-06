Algorithms implemented in ULA3
==============================

.. _nbar-algorithm-label:

Nadir BRDF-adjusted reflectance (NBAR)
--------------------------------------

An image taken by a satellite is affected by a range of factors, including:

* the angle at which a satellite views a pixel (view angle),
* illumination (determined largely by the solar angle),
* characteristics of the satellite, its sensor and its orbit,
* relative position of the sun and earth, and
* atmospheric conditions.

The NBAR algorithm attempts to normalise the imagery collected by the satellite to take account of these artifacts and hence render it more suitable for further analysis. The normalised imagery should appear as it would if the satellite were directly overhead (i.e. at `nadir <http://en.wikipedia.org/wiki/Nadir>`_) and at a specific solar angle (in the products produced here, 45°).

The NBAR algorithm was developed by staff at GA and CSIRO and is described in :download:`An Evaluation of the Use of Atmospheric and BRDF Correction to Standardize Landsat Data <auxiliary/li_etal_2010_05422912.pdf>`.



.. _tc-algorithm-label:

Terrain Correction (TC)
-----------------------

Another major factor which affects the appearance of remotely sensed imagery is terrain. This changes the angle of incidence and introduces shadows. The TC algorithm extends the NBAR algorithm to take account of these effects.

The algorithm is described in :download:`A physics-based atmospheric and BRDF correction for Landsat data over mountainous terrain <auxiliary/1-s2.0-S0034425712002544-main.pdf>`. The document :download:`terrain process <auxiliary/terrain_process.doc>` was provided by `Fuqin Li <mailto:fuqin.li@ga.gov.au>`_ describing the (Fortran) executables and scripts used in running the Terrain Correction algorithm (I think this is an updated version of a similar document which described the original NBAR processing). At the time of writing, the code could be found in various directories under :file:`\\win-satsan\TERRAIN_op`.



.. _pq-algorithm-label:

Pixel Quality (PQ)
------------------

Pixel Quality is a combined set of algorithms used to assess the quality of a pixel. Factors that can impede a pixels usefullness in analysis include environmental conditions such as cloud and cloud shadow, as well as sensor limitations such as saturation. Quality aspects included within this assessment are:
Saturation (Over and Under), Spectral Contiguity, Land/Sea discrimination, Cloud, Cloud Shadow.

The PQ documentation is described in :download:`Pixel Quality Algorithm Theoretical Basis Document <auxiliary/Pixel_Quality_ATBD.pdf>`.

Two cloud algorithms are incorporated into PQ. Automated Cloud Cover Assessment  (ACCA) and Function of mask (Fmask).
ACCA is described in :download:`Characterization of the Landsat-7 ETM+ Automated Cloud-Cover Assessment (ACCA) Algorithm <auxiliary/ACCA_Special_Issue_Final_Characterization_of_the_Landsat-7_ETM_Automated_Cloud-Cover_Assessment_Algorithm.pdf>`. 
Fmask is described at `Object-based cloud and cloud shadow detection in Landsat imagery <http://www.sciencedirect.com/science/article/pii/S0034425711003853>`_.
The cloud shadow algorithm is described in :download:`Cloud Shadow Algorithm Theoretical Basis Document <auxiliary/Cloud_Shadow_ATBD.pdf>`.


.. _fc-algorithm-label:

Fractional Cover (FC)
---------------------

The Fractional Cover algorithm utilises spectral bands 2,3,4,5,7 for either of the TM or ETM+ sensors. The spectral band inputs are sourced from the Australian Reflectance Grid 25m (ARG25); processed by Geoscience Australia. The bare soil, green vegetation and non-green vegetation endmenbers are calculated using models linked to an intensive field sampling program whereby more than 600 sites covering a wide variety of vegetation, soil and climate types were sampled to measure overstorey and ground cover following the procedure outlined in Muir et al (2011).
                    
Individual band information:

* PV: Green Cover Fraction. Fraction of green cover including green groundcover and green leaf material over all strata, within the Landsat pixel. A DN of 10000 is 100%.

* NPV: Non-green cover fraction. Fraction of non green cover including litter, dead leaf and branches over all strata, within the Landsat pixel. A DN of 10
000 is 100%.

* BS: Bare ground fraction. Fraction of bare ground including rock, bare and disturbed soil, within the Landsat pixel. A DN of 10
000 is 100%.

* UE: Unmixing Error. The residual error, defined as the Euclidean Norm of the Residual Vector. High values express less confidence in the fractional components.

Landcover fractions representing the proportions of green, non-green and bare cover retrieved by inverting multiple linear regression estimates and using synthetic endmembers in a constrained non-negative least squares unmixing model.

The following sources cover the algorithm theoretical development:

* Scarth, P., Roeder, A., and Schmidt, M., 2010, 'Tracking grazing pressure and climate interaction - the role of Landsat fractional cover in time series analysis' ,in Proceedings of the 15th Australasian Remote Sensing and Photogrammetry Conference, Alice Springs, Australia, 13-17 September 2010.

* Danaher, T., Scarth, P., Armston, J., Collet, L., Kitchen, J. & Gillingham, S. (2010). Remote sensing of tree-grass systems: The Eastern Australian Woodlands. In: Ecosystem Function in Savannas: Measurement and Modelling at Landscape to Global Scales, Eds. M.J. Hill and N.P. Hanan. CRC Press, Boca Raton. 

* Muir, J., Schmidt, M., Tindall, D., Trevithick, R., Scarth, P. & Stewart, J.B. (2011). Field measurement of fractional ground cover: a technical handbook supporting ground cover monitoring for Australia, prepared by the Queensland Department of Environment and Resource Management for the Australian Bureau of Agricultural and Resource Economics and Sciences, Canberra, November.
