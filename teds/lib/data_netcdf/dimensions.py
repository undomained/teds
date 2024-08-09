import numpy as np
import netCDF4 as nc

class Dimension:
    """
        The Dimension class deals with netcdf dimensions

        Two members:
        - self._name
        - self._value
        Methodes:
        - __init__(self,name, value, level)
        - __str__(self)
        - __repr__(self)
        - __eq__(self)
        - __ne__(self)
        - get_name(self)
        - get_value(self)
        - get_level(self)
        - set_name(self,name)
        - set_value(self,value)
        - set_level(self,level)
        - write(self, parent)
        - read(self, parent)
    """

    def __init__(self, name, value=None, level=2):
        """
            initialise Dimension class
            Arguments: 
                      - name:  name of the dimension
                      - value: value of the dimension
                      - level: indicates if dimension is atttached to main or to a group
            Three members:
            - self._name
            - self._value
            - self._level
        """
        self.set_name(name)
        self.set_value(value)
        self.set_level(level)

    def __str__(self):
        """
           Human readable printstatement.
        """
        pre = " "*self._level
        return f"{pre}-Dimension {self._name} with value {self._value}\n"

    def __repr__(self):
        """
            A valid Python expression that can be used to recreate the object.
        """
        return f"Dimension('{self._name}', {self._value})"

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def get_name(self):
        """
            Retrun the name of the dimension
        """
        return self._name

    def get_value(self):
        """
            Retrun the value of the dimension
        """
        return self._value

    def get_level(self):
        """
            Retrun the level of the dimension
        """
        return self._level

    def set_name(self,name):
        """
            Set the name of the dimension to the given name
        """
        self._name = name
        return

    def set_value(self,value):
        """
            Set the value of the dimension to the given value
        """
        self._value = value
        return

    def set_level(self,level):
        """
            Set the level of the dimension to the given level
        """
        self._level = level
        return

    def write(self, parent):
        """
            Write given dimension belonging to given parent to ntcdf file
        """
        parent.createDimension(self._name, self._value)
        return

    def read(self, parent):
        """
            Read given dimension belonging to given parent from ntcdf file
        """
        self._value = len(parent.dimensions[self._name])

        return 
