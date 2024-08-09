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
        - self._main_group
        Methodes:
        - __init__(self,file_name, title, mode)
        - __str__(self)
        - __repr__(self)
        - __eq__(self)
        - __ne__(self)
        - get_file_name(self)
        - get_title(self)
        - get_main_group(self)
        - set_logger(self, logger)
        - set_file_name(self, file_name)
        - set_title(self, title)
        - write(self)
        - read(self, file_name)
        - find(self, name, group, var, kind)
        - find_group(self, name)
        - update(self, name, value, group, var, kind)
        - get(self, name, group, var, kind)
        - remove(self, name, group, var, kind)
        - add(self, item, name, value, dimensions, attribute_list, group, var, kind)
    """

    def __init__(self, logger, file_name, title='', mode=None):
        """
            initialise DataNetCDF class
            Arguments: 
                      - logger:  logger pointer
                      - file_name:  name of the netcdf file
                      - title:  title of the netcdf file
                      - mode: of mode is 'r' the given file is readin
            Members:
            - self._logger
            - self._file_name
            - self._title
            - self._main_group
        """

        self.set_logger(logger)
        self.set_file_name(file_name)
        self.set_title(title)
        # group main is always there
        main = Group(self._logger, 'main')
        self._main_group = main

        if mode == 'r':
            self.read(self._file_name)

    def __str__(self):
        """
           Human readable printstatement.
        """

        data_string = f"Netcdf data with filename: {self._file_name} and title: {self._title}\n"
        data_string += str(self._main_group)

        return data_string

    def __repr__(self):
        """
            A valid Python expression that can be used to recreate the object.
        """
        return f"DataNetCDF('{self._file_name}', title='{self._title}')"

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

    def get_main_group(self):
        """
            Retrun the main group of the netcdf file
        """
        return self._main_group

    def write(self, file_name=None):
        """
            Write given all data to given file
            Note: file_name might be different from self._file_name
            if file_name = None, the output will written to self._file_name
            Also write the groups, and attributes, variables and dimensions 
        """
        output_file = self._file_name
        if file_name is not None:
            output_file = file_name
        self._logger.info(f"Writing to output file: {output_file}")
        with nc.Dataset(output_file, mode="w") as file_handle:
            file_handle.title=self._title
            self._main_group.write(file_handle)
        return 

    def read(self, file_name):
        """
            Read data from ntcdf file
            The main 'directory' in a netcdf file is treated as 
            a (special) group
        """
        self._file_name = file_name
        with nc.Dataset(file_name) as file_handle:

            # Getting meain stuff
            main = Group(self._logger, 'main')
            main.read(file_handle)
            self._main_group = main

        return

    def find(self, name, group=None, var=None, kind='variable'):
        """
            Get the kind of  object (default is Variable) corresponding to name in in given group.
            Variables, dimensions and attributes belong to a group. Attributes can also belong
            to Variables that belong to a group.
            If group is None the kind of object is searched in group main.
            First the group object is found, 
            Then the requested object list is obtained from the group,
            Then the list is searched to find the requested object with given name.
            Note: when kind is 'attribute' it can be an group attribute or a variable attribute.
            When the latter is the case also var needs to be given (name of the variable)
        """
        found_item = None

        grp = self.find_group(name=group)

        if grp is None:
            self._logger.warning(f"Given group: {groupname} can not be found. So {kind} with name {name} can not be found")
            return found_item

        # Find given object list
        if kind == 'variable':
            search_list = grp.get_variable_list()
        elif kind == 'dimension':
            search_list = grp.get_dimension_list()
        elif kind == 'attribute':
            # Group attribute or variable attribute?
            if var is not None:
                # we are looking for a variable attribute!
                var_list = grp.get_variable_list()
                for variable in grp.get_variable_list():
                   if variable.get_name() == var:
                       search_list = variable.get_attribute_list()
            else:
                # we are looking for a variable attribute!
                search_list = grp.get_attribute_list()

        # Look in list for object with given name
        for item in search_list:
            if item.get_name() == name:
                # item found
                found_item = item
    
        if found_item is None:
            var_string = ""
            if var is not None:
                var_string = f" for given {var}"
            self._logger.warning(f"{kind} with name {name} can not be found in group {group}{var_string}")

        return found_item

    def find_group(self, name):

        found_group = None

        if name is None:
            found_group = self._main_group
        else:
            grp = self._main_group
            found_group = grp.find_group(name)

        if found_group is None:
            self._logger.warning(f"Group {name} can not be found in netcdf file")

        return found_group

    def update(self, name, value, group=None, var=None, kind='variable'):
        """
            Update value for kind of information (default is variable) in a certain group.
            Note: if kind is 'attribute' it can belong to either the group or a Variable from the group.
            If latter is the case then also var (variable name) is needed.
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
        if kind == 'dimension':
            error_msg = "Removing a dimension is tricky. Also the variables that use it would need to be deleted. At the moment this is not possible"
            self._logger.error(error_msg)
            sys.exit(error_msg)
            return

        #Find the group that the item belongs to
        grp = self.find_group(name=group)
        # Find the corresponding kind list
        # Since we do not want to delete dimensions, the only kinds eligable
        # for removal are variables and attributes
        if kind == 'variable':
            search_list = grp.get_variable_list()
        elif kind == 'attribute':
            # Group attribute or variable attribute?
            if var is not None:
                # we are looking for a variable attribute!
                var_list = grp.get_variable_list()
                for variable in grp.get_variable_list():
                   if variable.get_name() == var:
                       search_list = variable.get_attribute_list()
            else:
                # we are looking for a variable attribute!
                search_list = grp.get_attribute_list()

        # Find item on list and remove from list
        found_index = 9999
        for index, item in enumerate(search_list):
            if item.get_name() == name:
                found_index = index
                break
        if found_index == 9999:
            var_string = ""
            if var is not None:
                var_string = " for variable {var}"
            self._logger.warning("{kind} with name {name} not found in group {group}{var_string}. Not removed!")
        if found_index != 9999:
            found_item = item
            del search_list[found_index]
        del found_item

        return

    def add(self, item=None, name=None, value=None, dimensions=None, attribute_list=None, group=None, var=None, kind='variable'):

        if item is not None and name is not None:
            warning_msg = f"Both item ({item}) and name ({name}) are set. I am confused and do not know what to do. No adding is performed"
            return
        if item is None and name is None:
            warning_msg = f"Both item ({item}) and name ({name}) are NOT set. I do not know what to do. No adding is performed"
            return

        #Find the group to which item needs to be added
        grp = self.find_group(name=group)
        level = grp.get_level()

        # Find the list the item needs to be added to
        if item is not None:
            item_type = type(item).__name__
            if kind is None:
                kind = item_type.lower()
            else:
                if kind != item_type.lower():
                    kind = item_type.lower()

        search_list = grp.get_list(kind, var=var)

        # Add to list
        if item is not None:
            item.set_level(level)
            search_list.append(item)
        else:
            #First create object, then add it to the list
            if kind == 'variable':
                dtype = value.dtype
                var = Variable(self._logger, name, value, dtype=dtype, dimensions=dimensions, attribute_list=attribute_list, level=level)
                search_list.append(var)
            elif kind == 'dimension':
                dim = Dimension(name, value, level=level)
                search_list.append(dim)
            elif kind == 'attribute':
                attr = Attribute(name, value, level=level)
                search_list.append(attr)
            elif kind == 'group':
                new_group = Group(self._logger, name, level=level+4)
                search_list.append(new_group)

#        # Set new list: should not be needed?????
#        grp.set_list(kind, search_list, var=var)

        return
