# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.
"""IO related operations."""

from importlib.resources import files
import yaml

from teds import log


def merge_dicts(config_full: dict, config_updates: dict) -> None:
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
        config_full. This is read from a user written input file.

    """
    for key in config_updates.keys():
        if key not in config_full:
            log.warning(f'unrecognized input parameter: {key}')
            continue
        if isinstance(config_updates[key], dict):
            merge_dicts(config_full[key], config_updates[key])
        else:
            config_full[key] = config_updates[key]


def merge_config_with_default(config: dict | None, teds_module: str) -> dict:
    """Merge user configuration with default configuration parameters.

    If user configuration dictionary is None then print a list of all
    valid configuration parameters and exit. That means a TEDS module
    was called without a configuration dictionary.

    Parameters
    ----------
    config
        Configuration dictionary written from a user written input file.
    teds_module
        TEDS module (example: 'teds.gm') that calls this function. A
        YAML file with the default settings is assumed to be in the
        directory of the module in a file called default_config.yaml.

    Returns
    -------
        Configuration settings merged with the default values.

    """
    default_config_path = str(files(teds_module) / 'default_config.yaml')
    if config is None:
        print(open(default_config_path).read())
        log.info(f'Stopping because the module {teds_module} was called '
                 'without an input file')
        exit(0)
    assert isinstance(config, dict)
    config_full: dict = yaml.safe_load(open(default_config_path))
    merge_dicts(config_full, config)
    return config_full
