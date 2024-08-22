import numpy as np

from teds.im.Python.algos.algo_base import Algorithm
from teds import log

class Wavemap(Algorithm):

    def __init__(self, algo_name="Wavemap"):
        
        self._algo_name = algo_name
        self._data = None


#    def check_input(self, image, ckd, wavelength):
#    def check_input(self, **input_data):
    def check_input(self, input_data):
        print(f"Check INPUT from {self._algo_name} class")

#    def execute(self, image, ckd, wavelength):
#    def execute(self, **input_data):
    def execute(self, input_data):

        #print(f"INPUT_DATA: {input_data}")  # for the love of god
        self._data = input_data.get_dataset('wavelength', c_name='measurement', kind='variable')

        return




