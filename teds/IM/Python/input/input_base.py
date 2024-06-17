
class Input:

    def __init__(self, logger, file_name):
        self._logger = logger
        self._file_name = file_name

    def __str__(self):
        """
            Human readable printstatement
        """
        input_string = f"Input file {self._file_name} is used\n"
        input_string += self.print()

        return input_string
       
    def open_file(self, mode='r') :
        """
            Open a file and return a filehandle
        """
        fh = open(self._file_name, mode)
        return fh

    def close_file(self,fh):
        """
            Close file of given filehandle
        """
        fh.close()
        return

    def read(self):
        """
            Read data from file using given file handle
        """
        fh = self.open_file()
        data = fh.readlines()
        self.close_file(fh)
        return data

    def print(self):
        """
            Loop over each line and create and return a string
        """
        data = self.read()
        print_string = ''
        for line in data:
            print_string += f"{line}\n"
        return print_string

    # Do we want this in Input class?
    def write(self, data):
        """
            Write data to file
        """
        fh = self.open_file(mode='w')
        fh.writelines(data)
        self.close_file(fh)
        return


