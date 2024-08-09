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
        - self._group_list
        - self._dimension_list
        - self._attribute_list
        - self._variable_list
        Methodes:
        - __init__(self,name, value, level, variable_list, dimension_list, attribute_list, group_list)
        - __str__(self)
        - __repr__(self)
        - __eq__(self)
        - __ne__(self)
        - __del__(self)
        - set_name(self, name)
        - set_level(self, level)
        - set_dimension_list(self, dimension_list)
        - set_attribute_list(self, attribute_list)
        - set_variable_list(self, variable_list)
        - set_group_list(self, group_list)
        - set_list(self, kind, new_list, var=None)
        - get_name(self)
        - get_level(self)
        - get_dimension_list(self)
        - get_attribute_list(self)
        - get_variable_list(self)
        - get_group_list(self)
        - get_list(self, kind, var)
        - write(self, parent)
        - read(self, parent, name)
        - find_in_group(self,name, var, kind)
        - find_group(self, name)
    """

    def __init__(self, logger, name, level=2, variable_list=None, dimension_list=None, attribute_list=None, group_list=None):
        """
            initialise Group class
            Arguments: 
                      - logger:  logger pointer
                      - name:  name of the group
                      - level:  level of the group
                      - diminsion_list: list of dimensions of the group
                      - attribute_list: attribute list of the group
                      - variable_list: variable list of the group
                      - group_list: group list of the group
            Members:
            - self._logger
            - self._name
            - self._level
            - self._attribute_list
            - self._variable_list
            - self._dimension_list
            - self._group_list
        """
        self.set_logger(logger)
        self.set_name(name)
        self.set_level(level)
        self.set_attribute_list(attribute_list)
        self.set_dimension_list(dimension_list)
        self.set_variable_list(variable_list)
        self.set_group_list(group_list)

    def __str__(self):
        """
           Human readable printstatement.
        """
        pre = " "*self._level

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
        group_string += f"{pre}## And the following group information:\n"
        for group_group in self._group_list:
            group_string += str(group_group)

        return group_string

    def __repr__(self):
        """
            A valid Python expression that can be used to recreate the object.
        """
        return f"Group('{self._name}', dimension_list={self._dimension_list}, attribute_list={self._attribute_list}, variable_list={self._variable_list}, group_list={self._group_list})"

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __del__(self):

        # delete corresponding variables
        for item in self._variable_list:
            del item
        self._variable_list = []
        # delete corresponding dimensions
        for item in self._dimension_list:
            del item
        self._dimension_list = []
        # delete corresponding attributes
        for item in self._attribute_list:
            del item
        self._attribute_list = []
        # delete corresponding groups
        for item in self._group_list:
            del item
        self._group_list = []

        return

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

    def get_group_list(self):
        """
            Retrun the group_list of the group
        """
        return self._group_list

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

    def get_level(self):
        """
            Retrun the level of the group
        """
        return self._level

    def get_list(self, kind, var=None):
        """
           Get specified kind of list
        """
        if kind == 'group':
            return self.get_group_list()
        elif kind == 'dimension':
            return self.get_dimension_list()
        elif kind == 'attribute':
            # Group attribute or variable attribute?
            if var is not None:
                # looking for a variable attribute!
                for variable in self.get_variable_list():
                   if variable.get_name() == var:
                       return variable.get_attribute_list()
            else:
                # we are looking for a variable attribute!
                return self.get_attribute_list()
        elif kind == 'variable':
            return self.get_variable_list()
        else:
            print(f"Oops, there is no list for kind: {kind}")

        return None

    def set_list(self, kind, new_list, var=None):
        """
           Set specified kind of list
        """
        if kind == 'group':
            return self.set_group_list(new_list)
        elif kind == 'dimension':
            return self.set_dimension_list(new_list)
        elif kind == 'attribute':
            return self.set_attribute_list(new_list)
            if var is not None:
                # looking for a variable attribute!
                for variable in self.get_variable_list():
                   if variable.get_name() == var:
                       return variable.set_attribute_list(new_list)
            else:
                # we are setting group attribute list!
                return self.set_attribute_list(new_list)
        elif kind == 'variable':
            return self.set_variable_list(new_list)
        else:
            print(f"Oops, there is no list for kind: {kind}")

        return None

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

    def set_group_list(self, group_list):
        """
            Set group list of the group to given group_list
        """
        self._group_list = group_list
        if group_list is None:
            self._group_list = []
        return

    def write(self, parent):
        """
            Write given group belonging to given parent to ntcdf file
            Also write the attributes, variables and dimensions belonging to this group
        """
        if self._name == 'main':
            group = parent
        else:
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
        for group_group in self._group_list:
            self._logger.debug("writing groups...")
            group_group.write(group)

        return 

    def read(self, parent):
        """
            Read given group belonging to given parent from ntcdf file
        """
        level = self._level
        if self._name == 'main':
            group = parent
        else:
            group = parent.groups[self._name]

        # Fill dimension list
        dimensions = group.dimensions
        self._dimension_list = []
        for dimname in dimensions:
            dim = Dimension(dimname, level=level)
            dim.read(group)
            self._dimension_list.append(dim)

        # Fill variable list
        variables = group.variables
        self._variable_list = []
        for varname in variables:
            var = Variable(self._logger, varname, level=level)
            var.read(group)
            self._variable_list.append(var)

        # Fill attribute list
        nattr = group.ncattrs() #???
        self._attribute_list = []
        for attrname in nattr:
            attr = Attribute(attrname, level=level)
            attr.read(group)
            self._attribute_list.append(attr)

        # Fill group list
        group_groups = group.groups
        self._group_list = []
        level = self._level+4
        # Dive deeper if needed
        for groupname in group_groups:
            group_group = Group(self._logger, groupname, level=level)
            group_group.read(group)
            self._group_list.append(group_group)

        return

    def find_in_group(self,name, var=None, kind='variable'):
        """
            Find kind of item in the present group
        """
        search_list = []
        if kind == 'variable':
            search_list = self._variable_list
        elif kind == 'dimension':
            search_list = self._dimension_list
        elif kind == 'attribute':
            if var is not None:
                # looking for a variable attribute!
               for variable in self._variable_list:
                   if variable.get_name() == var:
                       return variable.find(name)
            else:
                # looking for a group attribute!
                search_list = self._attribute_list
        else:
            possible_kinds = ['variable', 'attribute', 'dimension']
            self._logger.warning(f"The kind of information that you want to retrieve ({kind}) does not exist in netcdf. Possibilities: {possible_kinds}")
            return None

        for item in search_list:
            if item.get_name() == name:
                # item found
                return item
        #If we get here the item with name name is not found.
        self._logger.warning(f"Oops {kind} with name {name} is not found in group {self._name}!")
        return None


    def find_group(self, name):
        """
           Find group with given name in group list of this group.
           If name is of type group1/group2/group3:
           It is checked if group1 is in group list of this group
           and find_group is called for name=group2/group3
        """

        if "/" in name:
            groups = name.split("/")
            groupname = groups[0]
            group = "/".join(groups[1:] )
            for gr in self._group_list:
                if gr.get_name() == groupname:
                    return gr.find_group(group)
        else:
            groupname = name
            for gr in self._group_list:
                if gr.get_name() == groupname:
                    return gr
        return None
