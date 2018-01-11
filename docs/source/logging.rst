Logging
=======

Configuration for logging is controlled via a configuration file, which luigi parses and instantiates. Luigi can be pointed to the logging configuration file via the *luigi.cfg* under the **core** section, eg:

.. code-block:: cfg

   [core]
   logging_conf_file = /path/to/logging.cfg

See wagl's logging.cfg `logging.cfg <https://github.com/GeoscienceAustralia/wagl/blob/develop/configs/logging.cfg>`_ for the default wagl setup, additional loggers can be defined here and used in any wagl development work one might be likely to do.
The default setup is recommended as it'll capture any task errors to disk, as well as recording luigi's and wagl's general status. Any testing and development work can make use of the *wagl-status* qualname to retrieve the *wagl-status.log* initialised by luigi.

Any Task error's that occur when using either the *wagl.multifile_workflow* or the *wagl.singlefile_workflow* are captured, and written to the *wagl-errors.log* file. The logging performed by wagl is structured to output JSON formatted logging which is not only more readable, but allows a user to directly query the logs and retrieve the matching results using tools such as `jq <https://stedolan.github.io/jq/>`_. Errors will be output to *wagl-errors.log*.

Task history
------------

The luigi history and task status can also be stored in a database, thus providing another source of status information that in conjuction with the logs could be useful to any number of users.
To output the task history to an sqlite database, modify your *luigi.cfg* file with the following contents:

.. code-block:: cfg

   [scheduler]
   record_task_history = True
   
   [task_history]
   db_connection = sqlite:///luigi-task-history.db

And luigi will output the task history to an sqlite database named *luigi-task-history.db* relative to the directory that the scheduler was launched from.
