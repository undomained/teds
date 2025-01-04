# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

from importlib.resources import files
from yaml import safe_load


consts_file = files("teds.lib").joinpath("constants_outputvariables.yaml")
with open(consts_file, "r") as file:
    variable_dict = safe_load(file)


def writevariablefromname(grp, _name, dims, data):
    """Write variables to netcdf groups.

    Parameters
    ----------
    grp : Netcdf group
        Points to netcdf group
    _name : String
        Name of the variable.
    dims : String
        Dimensions of the variable.
    data : Any
        Data of the variable.

    """

    attr = variable_dict.get(_name)
    var = grp.createVariable(attr["name"], data.dtype, dims)
    for _ky, val in attr.items():
        if _ky != "name":
            var.setncattr(_ky, val)
    var[:] = data
    return var
