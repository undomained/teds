import time

class Datasets:

    def __init__(self, logger, list_name):
        self._logger = logger
        self._list_name = list_name
        self._datasets = []

        return

    def __str__(self):
        """
            Human readable printstatement
        """
        datasets_string = "The following datasets are included on this list:\n"
        for container in self._datasets:
            datasets_string += f"{container}\n"
        return datasets_string

    def add_container(self, c_name, container):

        c_type = type(container).__name__
        data_container = self.find_container( c_name)
        if data_container is None:
            # container does not yet exist, can be added.
            print(f"In ADD_CONTAINER {c_name}: c_type: {c_type}")
            data_container = Data_Container(self._logger, c_name, container, c_type)
            self._datasets.append(data_container)
        else:
            # This container already exists.
            # overwrite 
            print(f"WARNING: there exists already a container with name {c_name}. It will be overwritten")
            data_container.set_container(container)

        return

    def find_dataset(self, ds_name, c_name, group=None, kind=None):
        """
            Find dataset with given name ds_name in given group (if set) in given data container
            When data container, or group or dataset is not found, return None, otherwise
            return dataset. Could be object.
        """

        # First look for container
        container = self.find_container( c_name)
#        data_container = None
#        for container in self._datasets:
#            if container.get_name() == c_name:
#                data_container = container
#                break
        if container is None:
            self._logger.warning(f"Given data container {c_name} can not be found.")
            return None

        data_container = container.get_container()
        if container.get_type() == 'DataNetCDF':
            ds = data_container.get(ds_name, group=group, kind=kind)
            #Note: ds is still a Variable object
            return ds
        else:
            if group is not None:
                if group not in data_container:
                    self._logger.warning(f"Given group  {group} can not be found in given data container {c_name}.")
                    return None
                if ds_name not in data_container[group]:
                    self._logger.warning(f"Given dataset {ds_name} can not be found in group  {group} in given data container {c_name}.")
                    return None
                ds = data_container[group][ds_name]
                return ds
            else:
                if ds_name not in data_container:
                    self._logger.warning(f"Given dataset {ds_name} can not be found in given data container {c_name}.")
                    return None
                ds = data_container[ds_name]
                return ds
        return None



    def get_dataset(self, ds_name, c_name, group=None, kind=None):
        """
            Get dataset from given container and given group (if set).
        """
        dataset = None
        ds = self.find_dataset(ds_name, c_name, group=group, kind=kind)
        if ds is not None:
            if (type(ds).__name__ == 'Variable') or (type(ds).__name__ == 'Dimension'):
                dataset = ds.get_value(value)
            else:
                dataset = ds
        else:
            self._logger.warning(f"Dataset {ds_name} in data container {c_name} not found.")

        return dataset

    def find_container(self, c_name):

        data_container = None
        for container in self._datasets:
            if container.get_name() == c_name:
                data_container = container
                break
        if data_container is None:
            self._logger.warning(f"Given data container {c_name} can not be found.")

        return data_container

    def add_dataset(self, ds_name, c_name, data, group=None):

        container = self.find_container( c_name)
        if container is None:
            self._logger.warning(f"Given data container {c_name} can not be found. Dataset {ds_name} can not be added to it.")
            return None

        data_container = container.get_container()
        if container.get_type() == 'DataNetCDF':
            data_container.add(ds_name, value=data, group=group, kind='variable')
        else:
            if group is not None:
                if group not in data_container:
                    self._logger.warning(f"Given group  {group} can not be found in given data container {c_name}. Dataset {ds_name} can not be added.")
                    return None
                if ds_name in data_container[group]:
                    self._logger.warning(f"Given dataset {ds_name} already exists in group  {group} in given data container {c_name}. Dataset {ds_name} can not be added.")
                    return None
                data_container[group][ds_name] = data
            else:
                if ds_name in data_container:
                    self._logger.warning(f"Given dataset {ds_name} already exists in given data container {c_name}. Dataset {ds_name} can not be added.")
                    return None
                data_container[ds_name] = data
        return


    def update_dataset(self, ds_name, c_name, data, group=None, img=None):

        container = self.find_container( c_name)
        if container is None:
            self._logger.warning(f"Given data container {c_name} can not be found. Dataset {ds_name} can not be found and updated.")
            return None

        data_container = container.get_container()
        if container.get_type() == 'DataNetCDF':
            if img is not None:
                # Update one image in dataset
                # First get data
                data_to_update = data_container.get(ds_name, group=group, kind='variable')
                # Update data
                data_to_update[img,:,:] = data
                # Update dataset
                data_container.update(ds_name, value=data_to_update, group=group, kind='variable')
            else:
                # Update of whole dataset
                data_container.update(ds_name, value=data, group=group, kind='variable')
        else:
            if group is not None:
                if group not in data_container:
                    self._logger.warning(f"Given group  {group} can not be found in given data container {c_name}. Dataset {ds_name} can not be found and updated.")
                    return None
                if ds_name not in data_container[group]:
                    self._logger.warning(f"Given dataset {ds_name} can not be found in group  {group} in given data container {c_name}. Dataset {ds_name} can not be found and updated.")
                    return None
                data_container[group][ds_name] = data
            else:
                if ds_name not in data_container:
                    self._logger.warning(f"Given dataset {ds_name} can not be found in given data container {c_name}. Dataset {ds_name} can not be found and updated.")
                    return None
                data_container[ds_name] = data
        return

    def write(self):
        for container in self._datasets:
            data_container = container.get_container()
            c_name = container.get_name()
            self._logger.info(f"Writing for container {c_name}") 
            if container.get_type() == 'DataNetCDF':
                start_time_writing = time.perf_counter()
                self._logger.debug(f"{data_container}")
                self._logger.debug("#############")
                self._logger.debug("#############")

                data_container.write()
                end_time_writing = time.perf_counter()
                print(f"Writing data container {container.get_name()} to file took {(end_time_writing - start_time_writing):.6f}s")


class Data_Container:

    def __init__(self, logger, name, container, c_type):
        self._logger = logger
        self._name = name
        self._container = container
        self._type = c_type
        return

    def __str__(self):
        """
            Human readable printstatement
        """
        container_str = f"container {self._name} and content: {self._container}"
        return container_str

    def get_name(self):
        return self._name

    def get_type(self):
        return self._type

    def get_container(self):
        return self._container

    def set_container(self, container):
        self._container = container
        return


