# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

import logging as _logging
import sys

# Set up the default logger for this project
class TedsFormatter(_logging.Formatter):

    def __init__(self):
        super().__init__()
        formats = {
            _logging.DEBUG:
            '%(asctime)s : %(name)s : %(module)s : %(lineno)d : %(levelname)s : %(message)s',
            _logging.INFO:
            '%(asctime)s : %(name)s : %(module)s: %(levelname)s : %(message)s',
            _logging.WARNING:
            '%(asctime)s : %(name)s : %(module)s : %(lineno)d : %(levelname)s : %(message)s',
            _logging.ERROR:
            '%(asctime)s : %(name)s : %(module)s : %(lineno)d : %(levelname)s : %(message)s',
            _logging.CRITICAL:
            '%(asctime)s : %(name)s : %(module)s : %(lineno)d : %(levelname)s : %(message)s',
        }
        self.formatters = {}
        for f in formats:
            self.formatters[f] = _logging.Formatter(formats[f],
                                                    datefmt='%H:%M:%S')

    def format(self, record):
        log_fmt = record.levelno
        formatter = self.formatters[log_fmt]
        return formatter.format(record)


_handler = _logging.StreamHandler(sys.stdout)
_handler.setFormatter(TedsFormatter())

log = _logging.getLogger(__name__)
log.addHandler(_handler)
# Do we want DEBUG as the default???????
log.setLevel(_logging.DEBUG)
