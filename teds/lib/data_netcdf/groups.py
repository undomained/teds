import numpy as np
import netCDF4 as nc
from .dimensions import Dimension
from .variables import Variable
from .attributes import Attribute

class Group:
    """
        The Group class deals with netcdf groups

        Members:
        - self._logger
        - self._name
        - self._dimension_list
        - self._attribute_list
        - self._variable_list
        Methodes:
        - __init__(self,name, value)
        - get_name(self)
        - get_dimension_list(self)
        - get_attribute_list(self)
        - get_variable_list(self)
        - write(self, parent)
        - read(self, parent, name)
    """

#    def __init__(self, name, parent='main', variable_list= [], dimension_list = [], attribute_list = []):
    def __init__(self, logger, name, level='main', variable_list= None, dimension_list = None, attribute_list = None):
        """
            initialise Group class
            Arguments: 
                      - logger:  logger pointer
                      - name:  name of the group
                      - level:  level of the group
                      - diminsion_list: list of dimensions of the group
                      - attribute_list: attribute list of the group
                      - variable_list: variable list of the group
            Members:
            - self._logger
            - self._name
            - self._level
            - self._attribute_list
            - self._variable_list
            - self._dimension_list
        """
        self.set_logger(logger)
        self.set_name(name)
        self.set_level(level)
        self.set_attribute_list(attribute_list)
        self.set_dimension_list(dimension_list)
        self.set_variable_list(variable_list)

    def __str__(self):
        """
           Human readable printstatement.
        """
        n_indents = 2
        if self._level == 'group':
            n_indents *= 4
        pre = " "*n_indents

        group_string = f"{pre}##Group with name: {self._name}\n"
        group_string += f"{pre}### With dimensions:\n"
        for dimension in self._dimension_list:
            group_string += str(dimension)
        group_string += f"{pre}### With attributes:\n"
        for attribute in self._attribute_list:
            group_string += str(attribute)
        group_string += f"{pre}### With variables:\n"
        for variable in self._variable_list:
            group_string += str(variable)

        return group_string

    def __repr__(self):
        """
            A valid Python expression that can be used to recreate the object.
        """
        return f"Group('{self._name}', dimension_list={self._dimension_list}, attribute_list={self._attribute_list}, variable_list={self._variable_list})"

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def get_name(self):
        """
            Retrun the name of the group
        """
        return self._name

    def get_attribute_list(self):
        """
            Retrun the attribute_list of the group
        """
        return self._attribute_list

    def get_variable_list(self):
        """
            Retrun the variable_list of the group
        """
        return self._variable_list

    def get_dimension_list(self):
        """
            Retrun the dimension_list of the group
        """
        return self._dimension_list

    def set_logger(self,logger):
        """
            Set the logger of the group to the given logger
        """
        self._logger = logger
        return

    def set_name(self,name):
        """
            Set the name of the group to the given name
        """
        self._name = name
        return

    def set_level(self,level):
        """
            Set the level of the group to the given level
        """
        self._level = level
        return

    def set_dimension_list(self,dimension_list):
        """
            Set the dimension_list of the group to the given dimension_list
        """
        self._dimension_list = dimension_list
        if dimension_list is None:
            self._dimension_list = []
        return

    def set_attribute_list(self,attribute_list):
        """
            Set the attribute_list of the group to the given attribute_list
        """
        self._attribute_list = attribute_list
        if attribute_list is None:
            self._attribute_list = []
        return

    def set_variable_list(self,variable_list):
        """
            Set the variable_list of the group to the given variable_list
        """
        self._variable_list = variable_list
        if variable_list is None:
            self._variable_list = []
        return

    def write(self, parent):
        """
            Write given group belonging to given parent to ntcdf file
            Also write the attributes, variables and dimensions belonging to this group
        """
        group = parent.createGroup(self._name)
        self._logger.debug(f"wrting group: {self._name} to file")
        for attr in self._attribute_list:
            self._logger.debug(f"wrting attributes: for this group")
            attr.write(group)
        for dim in self._dimension_list:
            self._logger.debug(f"wrting dimensions: for this group")
            dim.write(group)
        for var in self._variable_list:
            self._logger.debug(f"wrting variables: for this group")
            var.write(group)
        return 

    def read(self, parent):
        """
            Read given group belonging to given parent from ntcdf file
        """
        group = parent.groups[self._name]
        # Fill dimension list
        dimensions = parent.groups[self._name].dimensions
        self._dimension_list = []
        for dimname in dimensions:
            dim = Dimension(dimname, level='group')
            dim.read(group)
            self._dimension_list.append(dim)


        # Fill variable list
        variables = parent.groups[self._name].variables
        self._variable_list = []
        for varname in variables:
            var = Variable(self._logger, varname, level='group')
            var.read(group)
            self._variable_list.append(var)

        # Fill attribute list
        nattr = parent.groups[self._name].ncattrs() #???
        self._attribute_list = []
        for attrname in nattr:
            attr = Attribute(attrname, level='group')
            attr.read(group)
            self._attribute_list.append(attr)

        return

    def find(self, name, var=None, kind='variable'):
        """
            Get the kind of object (default is variable) corresponding to name name.
            Note: if kind is 'attribute' also variable name is needed.
        """
        if kind == 'variable':
         search_list = self._variable_list
        elif kind == 'dimension':
         search_list = self._dimension_list
        elif kind == 'attribute':
            if var is not None:
                # we are looking for a variable attribute!
               for variable in self._variable_list:
                   if variable.get_name() == var:
                       return variable.find(name)
            else:
                search_list = self._attribute_list
        else:
            possible_kinds = ['variable', 'attribute', 'dimension']
            self._logger.warning(f"The kind of information that you want to retrieve ({kind}) does not exist in netcdf. Possibilities: {possible_kinds}")

        for item in search_list:
            if item.get_name() == name:
                # item found
                return item
        #If we get here the item with name name is not found.
        self._logger.warning(f"Oops {kind} with name {name} is not found in group {self._name}")
        return None

    def remove(self, name, var=None, kind='variable'):
        """
            remove kind of item with name name from netcdf data
            Note: if kind is 'attribute' also variable name is needed.
        """

        if kind == 'variable':
         search_list = self._variable_list
#        elif kind == 'dimension':
#         search_list = self._dimension_list
        elif kind == 'attribute':
            if var is not None:
                # we are looking for a variable attribute!
               for variable in self._variable_list:
                   if variable.get_name() == var:
                       return variable.remove(name)
            else:
                search_list = self._attribute_list
        else:
#            possible_kinds = ['variable', 'attribute', 'dimension']
            possible_kinds = ['variable', 'attribute']
            self._logger.warning(f"The kind of information that you want to remove ({kind}) does not exist in netcdf. Possibilities: {possible_kinds}")
            return

        found_index = 9999
        for index, item in enumerate(search_list):
            if item.get_name() == name:
                found_index = index
                break
        if found_index != 9999:
            del search_list[found_index]
            # Todo: Does this also work on original list?

        return

    def add(self, name=None, value=None, item=None, var=None, dimensions=None, attribute_list=[], kind='variable'):
        """
           Things that can be added to a group are: variables, dimensions or attributes
           Can be added as a name,var combination OR as an object
        """
        if name is None and item is None:
            self._logger.warning("Both the name and the item arguments are not set so nothing to add")
            return
        elif name is not None and item is not None:
            self._logger.warning("Both the name and the item arguments are set. Not sure what to use so nothing to add")
            return
        elif item is not None:
            item_type = type(item).__name__
            item_name = item.get_name()
            item.set_level('group')
            if item_type == 'Variable':
                self._variable_list.append(item)
            elif item_type == 'Dimensions':
                self._dimension_list.append(item)
            elif item_type == 'Attribute':
                self._attribute_list.append(item)

        elif name is not None:

            if kind == 'variable':
                dtype = value.dtype
                var = Variable(self._logger, name, value, dtype=dtype, dimensions=dimensions, attribute_list=attribute_list, level='group')
                self._variable_list.append(var)
            elif kind == 'dimension':
                dim = Dimension(name, value, level='group')
                self._dimension_list.append(dim)
            elif kind == 'attribute':
                if var is not None:
                   for variable in self._variable_list:
                       if variable.get_name() == var:
                           variable.add(name, value)
                else:
                    attr = Attribute(name, value, level='group')
                    self._attribute_list.append(attr)
            else:
                possible_kinds = ['variable', 'attribute', 'dimension']
                self._logger.warning(f"The kind of information that you want to add ({kind}) does not exist in netcdf. Possibilities: {possible_kinds}")

        return

