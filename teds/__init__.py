# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

import logging as _logging


# Set up the default logger for this project
class TedsFormatter(_logging.Formatter):

    def __init__(self):
        super().__init__()
        formats = {
            _logging.DEBUG:
            '[%(asctime)s] [%(levelname)s] [%(name)s:%(lineno)d] %(message)s',
            _logging.INFO:
            '[%(asctime)s] %(message)s',
            _logging.WARNING:
            '[%(asctime)s] [%(levelname)s] %(message)s',
            _logging.ERROR:
            '[%(asctime)s] [%(levelname)s] %(message)s',
            _logging.CRITICAL:
            '[%(asctime)s] [%(levelname)s] %(message)s',
        }
        self.formatters = {}
        for f in formats:
            self.formatters[f] = _logging.Formatter(formats[f],
                                                    datefmt='%H:%M:%S')

    def format(self, record):
        log_fmt = record.levelno
        formatter = self.formatters[log_fmt]
        return formatter.format(record)


_handler = _logging.StreamHandler()
_handler.setFormatter(TedsFormatter())

log = _logging.getLogger(__name__)
log.addHandler(_handler)
log.setLevel(_logging.DEBUG)
