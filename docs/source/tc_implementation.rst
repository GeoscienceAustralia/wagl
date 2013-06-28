Issues and Comments From Implementing TC
========================================

- Tight coupling between all modules - I want to run 'tc' and I get problems in 'fc'.

- ``CONFIG.code_root`` is not 'portable'. If I move 'common' relative to 'image_processor', things break. I've hacked in a fix for this that requires the setting of an environment variable (though I think I was clever enough not to break things if that variable is not defined.

- with call:
  ``../process.py --sequential --work /short/v10/tmp/nbar/nbar_work/sk/tc --process_level tc --tc /g/data/v10/tmp/sk_frac_tests/tc --l1t /g/data/v10/L1/2009-01/LS5_TM_OTH_P51_GALPGS01-002_092_086_20090115``

  I get error:
  ``option 'tc_data_root' is not defined in /home/547/sok547/workspace/ULA3/common/ULA3/common/processor.conf``
  not sure what this is needed for - but I added it just to get things going.

- It is not clear what ``<prefix>_path``, ``<prefix>_root`` and ``<prefix>_data_root`` are (meant to be) used for. While brief descriptions of the command line arguments to process.py are given in `ULA3.image_processor.ProcessorConfig`, I think more detail needs to be given on how some of the various 'derived' parameters are, well, derived and what they are used for.

- check_nbar_output.py does not take account of the input paths (I copied this to tc and it kept looking for the nbar output, even though the input paths we different).

- The catching of exceptions in the 'dynamic' import step and in process.py makes debugging hard, as it hides where exceptions occur (further, there seems to be a bug in the exception handling code at the bottom of process.py which raises an exception, further obfuscating things).

- Many errors produce a return code of 0.

- Assuming that the following files are functionally the same:
  - bilinear_ortho.f,
  - brdf_simbin.f

- Have updated fortran files (assuming functionality for previously used factors is identical to previous versions):
  - coefficient.f
  - read_modtran.f

- Where should the (national) dsm data go (regions are clipped from this)?

- Where are we going to setup test inputs (we need this to setup the unit tests)?

- Gathering of anciliary data does not work on small regions.

- :py:func:`ULA3.modtran.prepare_modtran_input.write_modis_brdf_files` and the word doc describing the process do not appear to be consistent in the order of bias and gain in brdf_modis_band<i>.txt. Need to check that these are correct.

- not sure if I've done the right thing with the threshold and average reflectance (last step - calculated not read from file. See image_processor.nbar.radiative_transfer.run_tc.process).

- I think each type of ancilliary data should have its own object type. This would allow the return values from the function to be documented properly and would help in making that code useful outside of the image processor.

- Having to know the type of the items stored in :py:class:`ULA3.data_manager.DataManager` is a real pain. I think some other mechanism for accessing them should be developed - perhaps a dictionary of dictionaries, the first indexed by name, second by type (and hence would only need to be accessed in cases of name clashes... which, come to think of it, would cause undefined behaviour anyway, since the same name implies the data gets serialised to the same file). I actually think the whole mechanism should be replaced... but that is a different matter.

- Terrain Correction depends on a 3 second DSM which has been preened by David Jupp. I have put this on VAYU (see config parameter DIR_DEM_TC in processor.conf), but it will need to be put on other systems also.

- The image processor seems to consume about 6~7GB of space that I cannot seem to identify. This may be something to do with VAYU (where I have observed it) or it may be something else.
