import logging

from structlog import wrap_logger
from structlog.processors import JSONRenderer

ERROR_LOGGER = wrap_logger(logging.getLogger('errors'),
                           processors=[JSONRenderer(indent=1, sort_keys=True)])

INTERFACE_LOGGER = logging.getLogger('luigi-interface')

STATUS_LOGGER = wrap_logger(logging.getLogger('status'),
                            processors=[JSONRenderer(indent=1, sort_keys=True)])

