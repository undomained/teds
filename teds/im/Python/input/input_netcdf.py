from teds.im.Python.input.input_base import Input
import teds.lib.data_netcdf.data_netcdf as dn

class Input_netcdf(Input):

    def __init__(self, file_name):
        self._file_name = file_name


    def read(self):
        """
            Read data from netcdf file and return netcdf_data object
        """
        netcdf_data = dn.DataNetCDF(self._file_name, mode='r')
        return netcdf_data

    def print(self):
        netcdf_data = self.read()
        data_str = ''
        data_str += str(netcdf_data)
        return data_str
