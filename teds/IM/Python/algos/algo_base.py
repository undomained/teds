class Algorithm:
    """
        Algortihm base class
        The members:
           - self._logger
           - self._algo_name: Name of this algorithm
           - self._data: data to be handled by the algorithm
        The Methodes:
        - __init__(self, logger, algo_name)
        - __str__(self)
        - check_input(self, input_data)
        - execute(self, input_data)
        - get_data(self)
        check_input and execute methods are to be overwritten in the Algorithm sub classes
    """

    def __init__(self, logger, algo_name):
        self._logger = logger
        self._algo_name = algo_name
        self._data = None
        self._stdev = None
        self._logger.debug("INIT ALGO_BASE")

    def __str__(self):
        """
            Human readable printstatement
        """
        algo_string = f"Algorithm {self._algo_name} is being applied\n"
        return algo_string

    def check_input(self, input_data):
        """
            Check the input data
        """
        self._logger.info("check input from base class Algorithm")

    def execute(self, input_data):
        """
            Execute the algorithm
        """
        self._data=input_data['image']
        self._logger.info("execute code from base class Algorithm")

    def get_data(self):
        """
            Retreive the data.
            Usually called after the algorithm has been run to obtained the updated data
        """
        return self._data
    def get_stdev(self):
        """
            Retreive the standard deviationdata.
            Usually called after the algorithm has been run to obtained the updated standard deviation data
        """
        return self._stdev

