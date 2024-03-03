import numpy as np
import netCDF4 as nc
from importlib.resources import files
from yaml import safe_load
from .attributes import Attribute

consts_file = files("teds.lib").joinpath("constants_outputvariables.yaml")
with open(consts_file, "r") as file:
    variable_dict = safe_load(file)
consts_ckd_file = files("teds.lib").joinpath("ckd_variables_no2_2_test.yaml")
with open(consts_ckd_file, "r") as file:
    variable_ckd_dict = safe_load(file)

class Variable:
    """
        The Variable class deals with netcdf variables

        Members:
        - self._logger
        - self._name
        - self._value
        - self._level
        - self._dimensions
        - self._dtype
        - self._attribute_list
        Methodes:
        - __init__(self,name, value)
        - get_name(self)
        - get_value(self)
        - get_level(self)
        - get_dimensions(self)
        - get_dtype(self)
        - get_attribute_list(self)
        - set_logger(self, logger)
        - set_name(self, name)
        - set_value(self, value)
        - set_level(self, level)
        - set_dimensions(self, dimension)
        - set_dtype(self, dtype)
        - set_attribute_list(self, attribute_list)
        - write(self, parent)
        - read(self, parent, name)
    """

    def __init__(self, logger, name, value=None, level='main', dimensions=None, dtype=None, attribute_list = None):
        """
            initialise Variable class
            Arguments: 
                      - logger:  logger pointer
                      - name:  name of the variable
                      - value: value of the variable
                      - level: indicates whether variable is attached to main or to a group
                      - dimensions: dimensions of the variable
                      - dtype: type of the variable
                      - attribute_list: attribute list of the variable
            Two members:
            - self._logger
            - self._name
            - self._value
            - self._level
            - self._dimensions
            - self._dtype
            - self._attribute_list
        """
        self.set_logger(logger)
        self.set_name(name)
        self.set_value(value)
        self.set_level(level)
        self.set_dimensions(dimensions)
        self.set_dtype(dtype)
        self.set_attribute_list(attribute_list)


        attributes = variable_dict.get(name)
        if attributes is not None:
            if 'name' in attributes:
                # SRON way: if name is set in the variables yaml file use that to name NetCDF variable
                self._name = attributes['name']
            # Now make Atrributes of the remaining attributes
            attributes.pop('name',None)
            for attribute, value in attributes.items():
                atr = Attribute(attribute, value)
                self._attribute_list.append(atr)
        attributes_ckd = variable_ckd_dict.get(name)
        if attributes_ckd is not None:
            if 'name' in attributes_ckd:
                # SRON way: if name is set in the variables yaml file use that to name NetCDF variable
                self._name = attributes_ckd['name']
            # Now make Atrributes of the remaining attributes
            attributes_ckd.pop('name',None)
            for attribute, value in attributes_ckd.items():
                atr = Attribute(attribute, value)
                self._attribute_list.append(atr)

        self._logger.debug(f"INIT var {self.__dict__}")


    def __str__(self):
        """
           Human readable printstatement.
        """
        n_indents = 2
        if self._level == 'group':
            n_indents *= 4
        pre = " "*n_indents
        var_string = f"{pre}- Variable with name: {self._name} with shape of data {self._value.shape}, dtype {self._dtype}, dimensions {self._dimensions}\n"
        if len(self._attribute_list)>0:
            var_string += f"{pre}### With attributes:\n"
            for attribute in self._attribute_list:
                var_string += str(attribute)
        else:
            var_string += f"{pre}### Without attributes\n"
        return var_string

    def __repr__(self):
        """
            A valid Python expression that can be used to recreate the object.
        """
        return f"Variable('{self._name}', {self._value}, dtype={self._dtype}, dimensions={self._dimensions}, attribute_list={self._attribute_list})"


    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def get_name(self):
        """
            Retrun the name of the variable
        """
        return self._name

    def get_value(self):
        """
            Retrun the value of the variable
        """
        return self._value

    def get_level(self):
        """
            Retrun the level of the variable
        """
        return self._level

    def get_dimensions(self):
        """
            Retrun the dimensions of the variable
        """
        return self._dimensions

    def get_dtype(self):
        """
            Retrun the dtype of the variable
        """
        return self._dtype

    def get_attribute_list(self):
        """
            Retrun the attribute_list of the variable
        """
        return self._attribute_list

    def write(self, parent):
        """
            Write given variable belonging to given parent to ntcdf file
            Also write the attributes belonging to this variable
        """
        var = parent.createVariable(self._name, self._dtype, self._dimensions)
        for attr in self._attribute_list:
            attr.write(var)
        var[:] = self._value
        return 

    def read(self, parent):
        """
            Read given variable information belonging to given parent from ntcdf file
        """
        self._value = parent.variables[self._name][:]
        self._dtype = self._value.dtype
        variable = parent.variables[self._name] #????
        self._dimensions = variable.dimensions
        nattr = variable.ncattrs()
        self._attribute_list = []
        for attrname in nattr:
            attr = Attribute(attrname,level='variable')
            attr.read(variable)
            self._attribute_list.append(attr)

        return

    def set_logger(self,logger):
        """
            Set the logger of the varibale to the given logger
        """
        self._logger = logger
        return

    def set_name(self,name):
        """
            Set the name of the varibale to the given name
        """
        self._name = name
        return

    def set_level(self,level):
        """
            Set the level of the variable to the given level
        """
        self._level = level
        return

    def set_value(self, value):
        self._value = value
        return
    
    def set_dtype(self, dtype):
        self._dtype = dtype
        return
    
    def set_dimensions(self, dimensions):
        self._dimensions = dimensions
        return
    
    def set_attribute_list(self, attribute_list):
        self._attribute_list = attribute_list
        if attribute_list is None:
            self._attribute_list = []
        return
    
    def find(self, name):
        """
            Get the attribute object corresponding to name name.
        """
        search_list = self._attribute_list
        for item in search_list:
            if item.get_name() == name:
                # item found
                return item
        #If we get here the item with name name is not found.
        self._logger.warning(f"Oops attribute with name {name} is not found in variable {self._name}")
        return None

    def remove(self, name):
        """
            remove attribute object corresponding to name name from netcdf file
        """
        search_list = self._attribute_list

        found_index = 9999
        for index, item in enumerate(search_list):
            if item.get_name() == name:
                found_index = index
                break
        if found_index != 9999:
            del search_list[found_index]
            # Todo: Does this also work on original list?

        return

    def add(self,name, value):
        """
            The only thing that can be added to a variable is an attribute
        """
        attr = Attribute(name, value)
        self._attribute_list.append(attr)

        return


