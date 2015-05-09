====================
Logging Memory usage
====================

When writing code for the NCI HPC platform it is extremely beneficial to 
keep the memory usage per process (specifically the Resident Set Size or RSS)
below 2 GB in total. 

If the RSS of a process is kept below 2GB, up to 16 of these processes may be run
in parallel on a single NCI node without the risk of exceeding the node's maximum 
physical memory (typically 32GB). In such cases node utilisation can approach 100%.

If the RSS of a process exceeds 2GB, it not be possible to run 16 processes on one node
and some of the node's CPUs will not be used. In such cases maximum node CPU utilisation will
be significantly less than 100%. For example, a process requiring 4GB can only expect to 
achieve a maximum node CPU utilisation of 50%.

Logging memory usage
--------------------

The simplest was to monitor memory usage during code deveopment is to log the
memory utilisation of the process during development

.. code:: python

    import gaip
    
    ...
    logging.debug("about to build big array, RSS=%f MB" % (gaip.get_rss_mbytes(), ))
    array = numpy.zeros_like(big_array)  
    logging.debug("now have big array, RSS=%f MB" % (gaip.get_rss_mbytes(), ))
    ...

The code above will report the total physical memory used by the process. Try different
coding strategies to keep the RSS value as low as possible.
