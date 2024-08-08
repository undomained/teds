import numpy as np
import netCDF4 as nc

class Attribute:
    """
        The Attribute class deals with netcdf attributes

        Two members:
        - self._name
        - self._value
        - self._parent
        Methodes:
        - __init__(self,name, value, level)
        - get_name(self)
        - get_value(self)
        - get_level(self)
        - set_name(self, name)
        - set_value(self, value)
        - set_level(self, level)
        - write(self, parent)
        - read(self, parent)
    """

    def __init__(self, name, value=None, level=2):
        """
            initialise Atrribute class
            Arguments: 
                      - name:  name of the attribute
                      - value: value of the attribute
                      - level: indicates if attribute is attached to main, to a group or to a variable
            Three members:
            - self._name
            - self._value
            - self._parent
        """
        self.set_name(name)
        self.set_value(value)
        self.set_level(level)

    def __str__(self):
        """
           Human readable printstatement.
        """
        pre = " "*self._level
        return f"{pre}-Attribute {self._name} with value {self._value}\n"

    def __repr__(self):
        """
            A valid Python expression that can be used to recreate the object.
        """
        return f"Attribute('{self._name}', {self._value})"

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def get_name(self):
        """
            Retrun the name of the attribute
        """
        return self._name

    def get_value(self):
        """
            Retrun the value of the attribute
        """
        return self._value

    def get_level(self):
        """
            Retrun the level of the attribute
        """
        return self._level

    def set_name(self,name):
        """
            Set the name of the attribute to the given name
        """
        self._name = name
        return

    def set_level(self,level):
        """
            Set the parent of the attribute to the given parent
        """
        self._level = level
        return

    def set_value(self,value):
        """
            Set the value of the attribute to the given value
        """
        self._value = value
        return


    def write(self, parent):
        """
            Write given attribute belonging to given parent to ntcdf file
        """
        parent.setncattr(self._name, self._value)
        return

    def read(self, parent):
        """
            Read given attribute belonging to given parent from ntcdf file
        """
        self._value = parent.getncattr(self._name)

        return 


