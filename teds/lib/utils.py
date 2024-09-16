# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""Commonly used utility functions.

"""

from teds import log


def merge_configs(config_full: dict, config_updates: dict) -> None:
    """Update a dictionary with contents of another dictionary.

    A dictionary, representing all valid configuration parameters of a
    given TEDS module, is merged with another dictionary such that
    only values that are present in both dictionaries are
    overwritten. The idea is to use a full set of configuration
    parameters with default values as a starting point and overwrite
    only those that the user has explicitly given.

    The configuration is recursively modified in place.

    Args:
      config_full:
        Full set of configuration parameters with their default values. It is
        recursively modified in place.
      config_updates:
        Set of configuration parameters that should be overwritten in
        config_full. This is read from the user-given input file.

    """
    for key in config_updates.keys():
        if key not in config_full:
            log.warning(f'unrecognized input parameter: {key}')
            continue
        if isinstance(config_updates[key], dict):
            merge_configs(config_full[key], config_updates[key])
        else:
            config_full[key] = config_updates[key]
