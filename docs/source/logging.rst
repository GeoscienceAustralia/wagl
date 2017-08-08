Logging
=======

Configuration for logging is controlled via a configuration file, which luigi parses and instantiates. Luigi can be pointed to the logging configuration file via the luigi.cfg under the **core** section, eg:

[core]
------
logging_conf_file = /path/to/logging.cfg

See gaip's logging.cfg :doc:`logging.cfg <../../configs/logging.cfg>` for the default gaip setup, additional loggers can be defined here and used in any gaip development work one might be likely to do.
The default setup is recommended as it'll capture any task errors to disk, as well as recording luigi's and gaip's general status. Any testing and development work can make use of the *gaip-status* qualname to retrieve the *gaip-status.log* initialised by luigi.

Any Task error's that occur when using either the *gaip.multifile_workflow* or the *gaip.singlefile_workflow* are captured, and written to the *gaip-errors.log* file. The logging performed by gaip is structured to output JSON formatted logging which is not only more readable, but allows a user to directly query the logs and retrieve the matching results using tools such as `jq <https://stedolan.github.io/jq/>`_. Errors will be output to *gaip-errors.log*.

Task history
------------

The luigi history and task status can also be stored in a database, thus providing another source of status information that in conjuction with the logs could be useful to any number of users.
To output the task history to an sqlite database, modify your luigi.cfg file with the following contents:

[scheduler]
-----------
record_task_history = True

[task_history]
--------------
db_connection = sqlite:///luigi-task-history.db

And luigi will output the task history to an sqlite database named luigi-task-history.db relative to the directory that the scheduler was launched from.
