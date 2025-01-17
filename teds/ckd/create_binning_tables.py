import sys, os
from teds import log
import numpy as np
import math
import teds.lib.data_netcdf.data_netcdf as dn
import teds.lib.lib_utils as Utils

def determine_binning(factor, n_detector, requested=None):
    """
        Given binning factor and number of detector bins
        determine number of bins and count table.

        In case of simple binning (requested not set):
        The number of binned bins is n_bins. 
        n_bins is floor(n_detector/factor)
        If there are detector bins remaining, then those are added
        together in the last (extra) bin. So n_bins += 1
        Last bin has different binning factor.

        Slightly less simple binning in case requested is given.
        Two cases:
        1.: requested <= n_bins
            Remaining detector bins are distributed over the first and last binned bin
            First and last bin have different binning factor
        2.: requested > n_bins
            Determine number of bins to be filled (n_to_be_filled = requested - n_bins)
            The remaining detector bins are distributed over the number of bins to be filled.
        
    """

    n_bins = math.floor(n_detector/factor)
    n_remainder = n_detector % factor
    print(f"n_bins: {n_bins} and n_remainder: {n_remainder}")

    if requested is not None:
        #Requested number of bins to end up with

        n_binned = requested
        count_bins = np.ones((requested,), dtype=np.uint16)
        count_bins *= factor

        # two possible scenarios: requested <= n_bins or > n_bins
        if requested <= n_bins:
            # distribute remaining detector rows/cols over two bins
            # so total of detector rows/cols for these two bins:
            total_remaining = n_remainder + 2 * factor
            # is total remaining odd or even?
            remainder = math.floor(total_remaining/2.)
            count_bins[0] = remainder
            if total_remaining % 2 == 0:
                # same number of rows in first and last binned row
                count_bins[-1] = remainder
            else:
                # one more row in last binned row
                count_bins[-1] = remainder+1

        elif requested > n_bins:
            # less binned rows/cols than requested.
            # distribute remainder over number of rows/cols that need to be filled
            # Those rows/cols to be filled will have different binning factor
            n_to_be_filled = requested - n_bins
            # how many rows with different BF at beginning and end?
            n_binned_different = math.floor(n_to_be_filled/2.)

            # in each n_to_be_filled there will be remainder rows.
            remainder = math.floor(n_remainder/n_to_be_filled)
            for i in range(n_binned_different):
                count_bins[i] = remainder
            if n_to_be_filled % 2 == 0:
                # equal number of binned rows at beginning of end with remaining factor
                for i in range(n_binned_different):
                    count_bins[-i] = remainder
                # if there are some detector bins remaining, add them to the last bin
                count_bins[-1] +=  n_remainder % n_to_be_filled
            else:
                # at the end one more binned row with remining factor
                for i in range(n_binned_different+1):
                    count_bins[-i] = remainder
                # if there are some detector bins remaining, add them to the last bin
                count_bins[-1] +=  n_remainder % n_to_be_filled

    else:
        # Nothing requested
        # put all remaining rows/cols (if there are any) in one binned row/col at the end
        factor_remaining = n_remainder

        if n_remainder > 0 :
            count_bins = np.ones((n_bins+1,), dtype=np.uint16)
            count_bins *= factor
            count_bins[-1] = factor_remaining
            n_binned = n_bins + 1
        else:
            count_bins = np.ones((n_bins,), dtype=np.uint16)
            count_bins *= factor
            n_binned = n_bins


    return n_binned, count_bins

def create_bin_list(n_detector, count_bins):
    """
       Based on count table information create 1D array
       which indicats to which binned bin a detector bin belongs
    """

    n_binned_index = np.zeros((n_detector,))
    index = 0
    start = 0
    stop = 0
    for nc in count_bins:
        stop = start + nc
        n_binned_index[start:stop] = index
        index += 1
        start += nc

    return n_binned_index



#def build_simple(bin_file, nrows, ncols, row_binning, col_factor=1, requested=None):
def build_simple(config):
    """
        Creating simple binning tables for given binning factors.
    """
    # Creating data_netcdf object to hold all the binning data information
    bin_data = dn.DataNetCDF(config['binning_table_file_name'])
    nrows = config['nrows']
    ncols = config['ncols']
    row_binnings = config['row_binnings']
    col_factor = config['col_factor']

    # Add dimension scales
    bin_data.add(name='column', value=ncols, kind='dimension') 
    bin_data.add(name='row', value=nrows, kind='dimension') 

    nbins_unbinned = nrows*ncols

    # Create the binning tables for the different row binnings
    for bf in row_binnings:

        row_factor = row_binnings[bf]['factor']
        requested_binned_rows = row_binnings[bf]['requested']
        log.info(f"Creating binning table for BF: {bf} with row factor: {row_factor} and requested binned rows: {requested_binned_rows}")

        # Now naming of binning table is based on row binning only.
        # when we also want to bin in col direction name should probably become somthing like:
        # f"Table_{int(row_factor)}_{int(col_factor)}"
        table_name = f"Table_{int(row_factor)}"
        #create group
        bin_data.add(name=table_name, kind='group') 

        # Obtain the number of binned rows and corresponding count_rows
        n_binned_rows, count_rows = determine_binning(row_factor, nrows, requested=requested_binned_rows)
        # Obtain the number of binned cols and corresponding count_cols
        n_binned_cols, count_cols = determine_binning(col_factor, ncols)

        # Determine the number of binned pixels
        bins = n_binned_rows * n_binned_cols
        log.info(f"Creating binning Table: {table_name} with number of  binned pixels: {bins}, n_binned_rows: {n_binned_rows} and n_bined_cols: {n_binned_cols}")
        bin_data.add(name='bins', value=bins, group=table_name, kind='dimension') 

        # Obtain count table for the binned pixels
        count_table = np.ones((bins,), dtype=np.uint16)
        for row in range(n_binned_rows):
            for col in range(n_binned_cols):
                count_table[row*(n_binned_cols) + col] = count_rows[row]*count_cols[col]

        # Get the binned index for the cols
        n_binned_col_index = create_bin_list(ncols, count_cols)
        # Get the binned index for the rows
        n_binned_row_index = create_bin_list(nrows, count_rows)

        # Create the binning table using the binned cols index and the binned row index
        binning_table = np.zeros((nrows,ncols), dtype=np.uint32)

        for row in range(nrows):
            bin_row = n_binned_row_index[row]
            factor = bin_row*n_binned_cols
            binning_table[row,:] = n_binned_col_index + factor

        # Add count table and binning table to the binning file
        bin_data.add(name='count_table', dimensions=('bins',), value=count_table, group=table_name, kind='variable')
        # Add the binning table
        bin_data.add(name='binning_table', dimensions=('row','column'), value=binning_table, group=table_name, kind='variable')

    # Done creating binning tables. Write data to file
    bin_data.write()
    log.info(f"Done creating and writing binning tables")

    return

if __name__ == '__main__' and __package__ is None:

    sys.path.append( os.path.dirname( os.path.dirname( os.path.abspath(__file__) ) ) )

    # Get configuration info
    cfgFile = sys.argv[1]
    config = Utils.get_config(cfgFile)

    build_simple(config)

