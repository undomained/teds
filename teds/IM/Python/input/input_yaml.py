import yaml
from teds.IM.Python.input.input_base import Input

class Input_yaml(Input):

    def __init__(self, file_name):
        self._file_name = file_name
        self._data = None


    def read(self):
        stream = self.open_file()
        data = yaml.safe_load(stream)
        self.close_file(stream)
        self._data = data
        return data

    def print(self):
        """
            Loop over entries in yaml and create and return a string
        """
        data = self._data
        yaml_string = f"Contents of yaml file {self._file_name}:"
        for key, value in data.items():
            yaml_string += f"{key}: {value}\n"
        return yaml_string

