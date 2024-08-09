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
        - self._fillvalue
        - self._attribute_list
        Methodes:
        - __init__(self, logger, name, value, level, dimensions, dtype, fillvalue, attribute_list)
        - __str__(self)
        - __repr__(self)
        - __eq__(self)
        - __ne__(self)
        - __del__(self)
        - get_name(self)
        - get_value(self)
        - get_level(self)
        - get_dimensions(self)
        - get_dtype(self)
        - get_fillvalue(self)
        - get_attribute_list(self)
        - set_logger(self, logger)
        - set_name(self, name)
        - set_value(self, value)
        - set_level(self, level)
        - set_dimensions(self, dimension)
        - set_dtype(self, dtype)
        - set_fillvalue(self, fillvalue)
        - set_attribute_list(self, attribute_list)
        - write(self, parent)
        - read(self, parent, name)
        - add(self,name, value)
    """

    def __init__(self, logger, name, value=None, level=2, dimensions=None, dtype=None, fillvalue=None, attribute_list = None):
        """
            initialise Variable class
            Arguments: 
                      - logger:  logger pointer
                      - name:  name of the variable
                      - value: value of the variable
                      - level: indicates whether variable is attached to main or to a group
                      - dimensions: dimensions of the variable
                      - dtype: type of the variable
                      - fillvalue: fill value for this variable
                      - attribute_list: attribute list of the variable
            Members:
            - self._logger
            - self._name
            - self._value
            - self._level
            - self._dimensions
            - self._dtype
            - self._fillvalue
            - self._attribute_list
        """
        self.set_logger(logger)
        self.set_name(name)
        self.set_value(value)
        self.set_level(level)
        self.set_dimensions(dimensions)
        self.set_dtype(dtype)
        self.set_fillvalue(fillvalue)
        self.set_attribute_list(attribute_list)


        attributes = variable_dict.get(name)
        if attributes is not None:
            if 'name' in attributes:
                # SRON way: if name is set in the variables yaml file use that to name NetCDF variable
                self._name = attributes['name']
            if 'FillValue' in attributes:
                # FillValue is a special attribute
                if self._fillvalue is None:
                    self.set_fillvalue(attributes['FillValue'])
                else:
                    warning_msg = f"A fill value was provided for initialization, but there is "\
                                  f"also a FillValue provided as attribute in the file {consts_file}! "\
                                  f"The init fill value takes precedent!"
                    self.logger.warning(warning_msg)

            # Now make Atrributes of the remaining attributes
            attributes.pop('name',None)
            # should we pop FillValue from the attributes? Probably.
            attributes.pop('FillValue',None)
            for attribute, value in attributes.items():
                atr = Attribute(attribute, value)
                self._attribute_list.append(atr)
        attributes_ckd = variable_ckd_dict.get(name)
        if attributes_ckd is not None:
            if 'name' in attributes_ckd:
                # SRON way: if name is set in the variables yaml file use that to name NetCDF variable
                self._name = attributes_ckd['name']
            if 'FillValue' in attributes_ckd:
                # FillValue is a special attribute
                if self._fillvalue is None:
                    self.set_fillvalue(attributes_ckd['FillValue'])
                else:
                    warning_msg = f"A fill value was provided for initialization, but there is "\
                                  f"also a FillValue provided as attribute in the file {consts_ckd_file}! "\
                                  f"The init fill value takes precedent!"
                    self.logger.warning(warning_msg)

            # Now make Atrributes of the remaining attributes
            attributes_ckd.pop('name',None)
            # should we pop FillValue from the attributes? Probably.
            attributes_ckd.pop('FillValue',None)
            for attribute, value in attributes_ckd.items():
                atr = Attribute(attribute, value)
                self._attribute_list.append(atr)

        self._logger.debug(f"INIT var {self.__dict__}")


    def __str__(self):
        """
           Human readable printstatement.
        """
        pre = " "*self._level
        var_string = f"{pre}- Variable with name: {self._name} with shape of data {self._value.shape}, dtype {self._dtype}, dimensions {self._dimensions}, fill value: {self._fillvalue}\n"
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
        return f"Variable('{self._name}', {self._value}, dtype={self._dtype}, dimensions={self._dimensions}, fillvalue={self._fillvalue}, attribute_list={self._attribute_list})"


    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __del__(self):

        # Delete the corresponding attibutes
        attr_list = self._attribute_list
        for item in attr_list:
           del item
        self._attribute_list = []

        return

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

    def get_fillvalue(self):
        """
            Retrun the fillvalue of the variable
        """
        return self._fillvalue

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
        if self._fillvalue is not None:
            # if a fill value has been set add it as an argument.
            var = parent.createVariable(self._name, self._dtype, self._dimensions, fill_value=self._fillvalue)
        elif self._dimensions is not None:
            var = parent.createVariable(self._name, self._dtype, self._dimensions)
        else:
            var = parent.createVariable(self._name, self._dtype)
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
        level = self._level+4
        for attrname in nattr:
            if attrname == '_FillValue':
                # Note this is a special attribute, do not add to attribute list but set the
                # self._fillvalue member.
                fillvalue = variable.getncattr(attrname)
                self.set_fillvalue(fillvalue)
            else:
                # Normal attribute, add to attributelist
                attr = Attribute(attrname,level=level)
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
    
    def set_fillvalue(self, fillvalue):
        self._fillvalue = fillvalue
        return
    
    def set_dimensions(self, dimensions):
        self._dimensions = dimensions
        return
    
    def set_attribute_list(self, attribute_list):
        self._attribute_list = attribute_list
        if attribute_list is None:
            self._attribute_list = []
        return
    
    def add(self,name, value):
        """
            The only thing that can be added to a variable is an attribute
        """
        attr = Attribute(name, value)
        self._attribute_list.append(attr)

        return


