"""
Logging configuration for wagl logs

Defines structured logging for:
    * Task message      -- qualname task
    * Status messages   -- qualname status
    * Luigi interface   -- qualname luigi-interface
"""

import logging

import structlog
from structlog.processors import JSONRenderer


COMMON_PROCESSORS = [
    structlog.stdlib.add_log_level,
    structlog.processors.TimeStamper(fmt="ISO"),
    structlog.processors.StackInfoRenderer(),
    structlog.processors.format_exc_info,
    JSONRenderer(sort_keys=True),
]


def get_wrapped_logger(logger_name: str = "root", **kwargs):
    """ Returns a struct log equivalent for the named logger """
    return structlog.wrap_logger(
        logging.getLogger(logger_name), COMMON_PROCESSORS, **kwargs
    )


class FormatJSONL(logging.Formatter):
    """ Prevents printing of the stack trace to enable JSON lines output """

    def formatException(self, ei):
        """ Disables printing separate stack traces """
        return


TASK_LOGGER = get_wrapped_logger("task")
STATUS_LOGGER = get_wrapped_logger("status")

INTERFACE_LOGGER = logging.getLogger("luigi-interface")
