import numpy as np
import netCDF4 as nc
from .variables import Variable
from .groups import Group
from .dimensions import Dimension
from .attributes import Attribute

class DataNetCDF:
    """
        The DataNetCDF class deals with netcdf data

        The members:
        - self._file_name
        - self._title
        - self._group_list
        - self._dimension_list
        - self._attribute_list
        - self._variable_list
        Methodes:
        - __init__(self,file_name, title, mode, group_list, variable_list, dimension_list, attribute_list)
        - get_file_name(self)
        - get_title(self)
        - get_attribute_list(self)
        - get_dimension_list(self)
        - get_variable_list(self)
        - get_group_list(self)
        - set_file_name(self, file_name)
        - set_title(self, title)
        - set_attribute_list(self, attribute_list)
        - set_dimension_list(self, dimension_list)
        - set_variable_list(self, variable_list)
        - set_group_list(self, group_list)
        - write(self)
        - read(self, file_name)
    """

    def __init__(self, logger, file_name, title='', mode=None, group_list=None, variable_list= None, dimension_list = None, attribute_list = None):
        """
            initialise DataNetCDF class
            Arguments: 
                      - logger:  logger pointer
                      - file_name:  name of the netcdf file
                      - title:  title of the netcdf file
                      - diminsion_list: list of dimensions of the netcdf file
                      - attribute_list: attribute list of the netcdf file
                      - variable_list: variable list of the netcdf file
                      - group_list: group list of the netcdf file
            Members:
            - self._logger
            - self._file_name
            - self._title
            - self._group_list
            - self._attribute_list
            - self._variable_list
            - self._dimension_list
        """

        self.set_logger(logger)
        self.set_file_name(file_name)
        self.set_title(title)
        self.set_attribute_list(attribute_list)
        self.set_dimension_list(dimension_list)
        self.set_variable_list(variable_list)
        self.set_group_list(group_list)

        if attribute_list is None:
            self._attribute_list = []
        if dimension_list is None:
            self._dimension_list = []
        if variable_list is None:
            self._variable_list = []
        if group_list is None:
            self._group_list = []

        if mode == 'r':
            self.read(self._file_name)

    def __str__(self):
        """
           Human readable printstatement.
        """

        data_string = f"Netcdf data with filename: {self._file_name} and title: {self._title}\n"
        data_string += "## And dimensions:\n"
        for dimension in self._dimension_list:
            data_string += str(dimension)
        data_string += "## And attributes:\n"
        for attribute in self._attribute_list:
            data_string += str(attribute)
        data_string += "## And variables:\n"
        for variable in self._variable_list:
            data_string += str(variable)
        data_string += "## And the following group information:\n"
        for group in self._group_list:
            data_string += str(group)

        return data_string

    def __repr__(self):
        """
            A valid Python expression that can be used to recreate the object.
        """
        return f"DataNetCDF('{self._file_name}', title='{self._title}', group_list={self._group_list}, dimension_list={self._dimension_list}, attribute_list={self._attribute_list}, variable_list={self._variable_list})"

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def set_logger(self, logger):
        """
            Set logger
        """
        self._logger = logger
        return

    def set_file_name(self, file_name):
        """
            Set file name
        """
        self._file_name = file_name
        return

    def set_title(self, title):
        """
            Set title
        """
        self._title = title
        return

    def set_attribute_list(self, attribute_list):
        """
            Set attribute list
        """
        self._attribute_list = attribute_list
        return

    def set_dimension_list(self, dimension_list):
        """
            Set dimension list
        """
        self._dimension_list = dimension_list
        return

    def set_variable_list(self, variable_list):
        """
            Set variable list
        """
        self._variable_list = variable_list
        return

    def set_group_list(self, group_list):
        """
            Set group list
        """
        self._group_list = group_list
        return

    def get_file_name(self):
        """
            Retrun the file_name of the netcdf file
        """
        return self._file_name

    def get_title(self):
        """
            Retrun the title of the netcdf file
        """
        return self._title

    def get_group_list(self):
        """
            Retrun the group_list of the netcdf file
        """
        return self._group_list

    def get_attribute_list(self):
        """
            Retrun the attribute_list of the netcdf file (main)
        """
        return self._attribute_list

    def get_variable_list(self):
        """
            Retrun the variable_list of the netcdf file (main)
        """
        return self._variable_list

    def get_dimension_list(self):
        """
            Retrun the dimension_list of the netcdf file (main)
        """
        return self._dimension_list

    def write(self, file_name=None):
        """
            Write given all data to given file
            Note: file_name might be different from self._file_name
            if file_name = None, the output will written to self._file_name
            Also write the groups, attributes, variables and dimensions 
        """
        output_file = self._file_name
        if file_name is not None:
            output_file = file_name
        self._logger.info(f"Writing to output file: {output_file}")
        with nc.Dataset(output_file, mode="w") as file_handle:
            file_handle.title=self._title
            for attr in self._attribute_list:
                attr.write(file_handle)
            for dim in self._dimension_list:
                self._logger.debug("writing dimensions...")
                self._logger.debug(f"...{dim}...")
                dim.write(file_handle)
            for var in self._variable_list:
                var.write(file_handle)
            for group in self._group_list:
                self._logger.debug("writing groups...")
                group.write(file_handle)
        return 

    def read(self, file_name):
        """
            Read data from ntcdf file
        """
        self._file_name = file_name
        with nc.Dataset(file_name) as file_handle:

            # Fill dimension list
            dimensions = file_handle.dimensions
            self._dimension_list = []
            for dimname in dimensions:
                dim = Dimension(dimname)
                dim.read(file_handle)
                self._dimension_list.append(dim)
    
    
            # Fill variable list
            variables = file_handle.variables
            self._variable_list = []
            for varname in variables:
                var = Variable(self._logger, varname)
                var.read(file_handle)
                self._variable_list.append(var)
    
            # Fill attribute list
            nattr = file_handle.ncattrs() #???
            self._attribute_list = []
            for attrname in nattr:
                attr = Attribute(attrname)
                attr.read(file_handle)
                self._attribute_list.append(attr)
    
            # Fill group list
            groups = file_handle.groups #???
            self._group_list = []
            for groupname in groups:
                group = Group(self._logger, groupname)
                group.read(file_handle)
                self._group_list.append(group)

        return

    def find(self, name, group=None, var=None, kind='variable'):
        """
            Get the kind of  object (default is variable) corresponding to name in in given group.
            If group is None the variables in main are searched.
            Note: if kind is 'attribute' also variable name is needed.
        """
        if group is None:
            if kind == 'variable':
             search_list = self._variable_list
            elif kind == 'dimension':
             search_list = self._dimension_list
            elif kind == 'attribute':
             search_list = self._attribute_list
            elif kind == 'group':
             search_list = self._group_list
            else:
                possible_kinds = ['variable', 'attribute', 'dimension', 'group']
                self._logger.warning(f"The kind of information that you want to retrieve ({kind}) does not exist in netcdf. Possibilities: {possible_kinds}")

            for item in search_list:
                if item.get_name() == name:
                    # item found
                    return item
            #If we get here the item with name name is not found.
            self.logger.warning(f"Oops {kind} with name {name} is not found in main directory of file {self._file_name}!")
            return None
        else:
            for gr in self._group_list:
                if gr.get_name() == group:
                    return gr.find(name, var=var, kind=kind)
            #If we get here the group has not been found
            self._logger.warning(f"Oops group with name {group} is not found in main directory of file {self._file_name}! So {kind} with name {name} can not be found.")
            return None

        return None

    def update(self, name, value, group=None, var=None, kind='variable'):
        """
            Update value for kind of information (default is variable) in a certain group.
            Note: if kind is 'attribute' also variable name is needed.
        """

        item = self.find(name, group=group, var=var, kind=kind)
        if item is not None:
            item.set_value(value)
        else:
            self._logger.warning(f"{kind} with name {name} not found in group {group}. Value can not be changed")
        return

    def get(self, name, group=None, var=None, kind='variable'):
        """
            Get value for kind of information (default is variable) in a certain group.
            Note: if kind is 'attribute' also variable name is needed.
        """
        item = self.find(name, group=group, var=var, kind=kind)

        if item is not None:
            return item.get_value()
        else:
            self._logger.warning(f"{kind} with name {name} not found in group {group}. Value can not be returned")
        return None

    def remove(self, name, group=None, var=None, kind='variable'):
        """
            remove kind of item with name name from netcdf data
            Note: if kind is 'attribute' also variable name is needed.
            Note2: if a dimension is removed it should also be removed from all
                   the variables that use it. Not sure if that is a good idea.
                   So for the moment it is not possible to remove a dimension
        """
        if group is None:
            if kind == 'variable':
             search_list = self._variable_list
        #    elif kind == 'dimension':
        #     search_list = self._dimension_list
            elif kind == 'attribute':
             search_list = self._attribute_list
            elif kind == 'group':
             search_list = self._group_list
            else:
        #        possible_kinds = ['variable', 'attribute', 'dimension', 'group']
                possible_kinds = ['variable', 'attribute', 'group']
                self._logger.warning(f"The kind of information that you want to remove ({kind}) does not exist in netcdf. Possibilities: {possible_kinds}. Nothing removed.")
                return

            found_index = 9999
            for index, item in enumerate(search_list):
                if item.get_name() == name:
                    found_index = index
                    break
            if found_index != 9999:
                del search_list[found_index]
                # Todo: Does this also work on original list?

        else:

            for gr in self._group_list:
                if gr.get_name() == group:
                    return gr.remove(name, var=var, kind=kind)
            #If we get here the group has not been found
            self._logger.warning(f"Oops group with name {group} is not found in main directory of file {self._file_name}. Not able to remove requested {kind} with name {name}")

        return

    def add(self, name=None, value=None, item=None, group=None, var=None, dimensions=None, attribute_list=None, kind='variable'):
        """
           Things that can be added to a main are: groups, variables, dimensions or attributes
           Can be added as a name,var combination OR as an object
        """

        if name is None and item is None:
            self._logger.warning("Both the name and the item arguments are not set so nothing to add")
            return
        elif name is not None and item is not None:
            self._logger.warning("Both the name and the item arguments are set. Not sure what to use so nothing to add")
            return
        elif item is not None:
            if group is not None:
                for gr in self._group_list:
                    if gr.get_name() == group:
                        gr.add(name=name, value=value, item=item, var=var, dimensions=dimensions, attribute_list=attribute_list, kind=kind)
                        return
                #If we get here the group has not been found
                item_type = type(item).__name__
                item_name = item.get_name()
                self._logger.warning(f"Oops group with name {group} is not found in main directory of file {self._file_name}. Not able to add requested item of kind {item_type} and name {item_name}")
            else:
                # Look in main
                item_type = type(item).__name__
                item_name = item.get_name()
                item.set_level('main')
                if item_type == 'Variable':
                    self._variable_list.append(item)
                elif item_type == 'Dimensions':
                    self._dimension_list.append(item)
                elif item_type == 'Attribute':
                    self._attribute_list.append(item)
                elif item_type == 'Group':
                    self._group_list.append(item)

        elif name is not None:

            if group is not None:

                for gr in self._group_list:
                    if gr.get_name() == group:
                        gr.add(name, value=value, var=var, dimensions=dimensions, attribute_list=attribute_list, kind=kind)
                        return
                #If we get here the group has not been found
                self._logger.warning(f"Oops group with name {group} is not found in main directory of file {self._file_name}. Not able to add requested {kind} with name {name}")

            else:

                if kind == 'variable':
                    dtype = value.dtype
                    var = Variable(self._logger, name, value, dtype=dtype, dimensions=dimensions, attribute_list=attribute_list, level='main')
                    self._variable_list.append(var)
                elif kind == 'dimension':
                    dim = Dimension(name, value, level='main')
                    self._dimension_list.append(dim)
                elif kind == 'attribute':
                    attr = Attribute(name, value, level='main')
                    self._attribute_list.append(attr)
                elif kind == 'group':
                    grpr = Group(self._logger, name, level='main')
                    self._group_list.append(grpr)
                else:
                    possible_kinds = ['variable', 'attribute', 'dimension', 'group']
                    self._logger.warning(f"The kind of information that you want to add ({kind}) does not exist in netcdf. Possibilities: {possible_kinds}")

        return



