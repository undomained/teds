import sys, os
import logging
import numpy as np
import math
import lib.data_netcdf.data_netcdf as dn

def get_logger():
    """
       Gets or creates a logger
    """

    log_level = logging.INFO
    log_format = '%(asctime)s : %(name)s : %(module)s : %(lineno)d : %(levelname)s : %(message)s'
    log_formatter = logging.Formatter(log_format)
    date_format = '%d/%m/%Y %H:%M:%S'

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    c_handler = logging.StreamHandler(sys.stdout)
    c_handler.setFormatter(log_formatter)
    logger.addHandler(c_handler)

    return logger

def build_simple(logger, bin_file, nrow, ncols, row_binning, col_factor=1):
    """
        Creating simple binning tables for given binning factors.
    """
    # Creating data_netcdf object to hold all the binning data information
    bin_data = dn.DataNetCDF(logger, bin_file)

    # Add dimension scales
    bin_data.add(name='column', value=ncols, kind='dimension') 
    bin_data.add(name='row', value=nrows, kind='dimension') 

    nbins_unbinned = nrows*ncols

    # Create the binning tables for the different row binnings
    for row_factor in row_binning:
        table_name = f"Table_0{row_factor}"
        #create group
        bin_data.add(name=table_name, kind='group') 

        # Add variable lineskip_arr, which depends on nrows (number of unbinned rows)
        lineskip_arr = np.ones((nrows,), dtype=np.int8)
        bin_data.add(name='lineskip_arr', dimensions=('row',), value=lineskip_arr, group=table_name, kind='variable')

        n_binned_cols = math.floor(ncols/col_factor)
        n_binned_rows = math.floor(nrows/row_factor)

        row_remainder = nrows % row_factor
        col_remainder = ncols % col_factor
        # Assume for now that all rows/cols remaining will be binned into 1 row/col
        row_factor_remaining = row_remainder
        col_factor_remaining = col_remainder

        n_binned_rows_total = n_binned_rows
        n_binned_cols_total = n_binned_cols
        if row_remainder != 0:
            n_binned_rows_total += 1
        if col_remainder != 0:
            n_binned_cols_total += 1

        # create number of binned pixels dimension
        bins = n_binned_rows_total * n_binned_cols_total
        logger.info(f"Creating binning Table: {table_name} with number of  binned pixels: {bins}, n_binned_rows_total: {n_binned_rows_total} and n_bined_cols_total: {n_binned_cols_total}")
        bin_data.add(name='bins', value=bins, group=table_name, kind='dimension') 

        # Note: SRON IM code kinda assumes no binning in spectral dimension since there is only one count_table table that
        # contains row_factor, i.e. binning on row direction
        count_table = np.ones((bins,), dtype=np.uint16)
        count_table[0:n_binned_rows] = row_factor
        count_table[n_binned_rows:] = row_factor_remaining
        bin_data.add(name='count_table', dimensions=('bins',), value=count_table, group=table_name, kind='variable')

        # First create the first binned row pixel numbering.
        # Taking into account that also binning is possible in the col dimension
        bin_row = np.ones((ncols,))
        unbinned_col = 0
        for col in range(n_binned_cols):
            start = unbinned_col
            while unbinned_col < start+col_factor:
                bin_row[unbinned_col] =  col
                unbinned_col += 1
        # Note using the fact that we do some simple binning and that if there were columns remaining they are all binned together.
        if col_remainder != 0:
            start = unbinned_col
            while unbinned_col < start+col_factor_remaining:
                bin_row[unbinned_col] =  n_binned_cols
                unbinned_col += 1

        # Now create the binning table
        # Do the same as done for the columns, but now in row direction
        binning_table = np.ones((nrows,ncols), dtype=np.uint32)

        binned_rows = range(n_binned_rows)
        unbinned_row = 0
        for row in range(n_binned_rows):
            start = unbinned_row
            while unbinned_row < start+row_factor:
                binning_table[unbinned_row,:] = bin_row+(row*n_binned_cols_total)
                unbinned_row += 1
        # Note using the fact that we do some simple binning and that if there were rows remaining they are all binned together.
        if row_remainder != 0:
            start = unbinned_row
            while unbinned_row < start+row_factor_remaining:
                binning_table[unbinned_row,:] = bin_row+(n_binned_rows*n_binned_cols_total)
                unbinned_row += 1

        # Add the binning table
        bin_data.add(name='binning_table', dimensions=('row','column'), value=binning_table, group=table_name, kind='variable')

    # Done creating binning tables. Write data to file
    bin_data.write()
    logger.info(f"Done creating and writing binning tables")

    return

if __name__ == '__main__' and __package__ is None:

    sys.path.append( os.path.dirname( os.path.dirname( os.path.abspath(__file__) ) ) )

    build_logger = get_logger()

    binning_table_file_name = "../data/no2/ckd/binning_table_no2.nc"
    # For now:
    nrows = 2870 
    ncols = 1681
#    nrows = 512
#    ncols = 640
    row_binnings = [1,2,3,4,5,6]
    col_factor = 1
    # TODO test case when col_factor != 1
    build_simple(build_logger, binning_table_file_name, nrows, ncols, row_binnings, col_factor)

