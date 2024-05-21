#include "header.h"
#include "functions.h"
#include "netcdf_object.h"
#include "binningtable.h"

// Binning table file. All binning tables are included in one file.
// That file governs the pre-binning, which is applied for all non-trivial
// binning tables. This means that all binning tables used at the same time
// use the same prebinning, except for possibly the trivial binning table
// that is just 'no binning at all, even no prebinning'. A binning table
// with just pre-binning should be defined with a table that does no
// furtherm binning.

// Constructor.
Binningtable_file::Binningtable_file( // {{{
    Logger *creator
) : Logger(creator)
{
} // }}}
Binningtable_file::~Binningtable_file() {} // Destructor.

int Binningtable_file::read( // {{{
    string &binningtable_filename, // File with all the used binning tables.
    size_t a_dim_detector_spat,
    size_t a_dim_detector_spec
)
{

    // Open the NetCDF file. Let us hope we can keep it open. It will be
    // closed when the L1A manager dies.
    nc = make_unique<NetCDF_object>(this);
    handle(nc->open(binningtable_filename,NcFile::read));
    // Read dimensions in main group.
    netcdf_check(nc,dim_binningtable_spat = nc->ncid->getDim("row").getSize());
    netcdf_check(nc,dim_binningtable_spec = nc->ncid->getDim("column").getSize());
    // Verify that the pre-binning is with integer factor.
    dim_detector_spat = a_dim_detector_spat;
    dim_detector_spec = a_dim_detector_spec;
    prebinning_spat = dim_detector_spat / dim_binningtable_spat; // This must be an exact division.
    check_error(prebinning_spat * dim_binningtable_spat != dim_detector_spat,"Error: Prebinning factor in spatial (row) dimension is not integral: %zu to %zu.",dim_detector_spat,dim_binningtable_spat);
    prebinning_spec = dim_detector_spec / dim_binningtable_spec; // This must be an exact division.
    check_error(prebinning_spec * dim_binningtable_spec != dim_detector_spec,"Error: Prebinning factor in spectral (column) dimension is not integral: %zu to %zu.",dim_detector_spec,dim_binningtable_spec);
    // Leave the file open, it can be used by the binning tables themselves.

    return 0;

} // }}}

// Constructor.
Binningtable::Binningtable( // {{{
    Logger *creator
) : Logger(creator)
{
} // }}}
Binningtable::~Binningtable() {} // Destructor.

// Binning table read routine.
int Binningtable::read( // {{{
    Binningtable_file *bin_file, // Binning table file.
    uint8_t id // Mask identifier.
)
{

    npix_unbinned = bin_file->dim_detector_spat * bin_file->dim_detector_spec;

    NetCDF_object *nc = bin_file->nc.get(); // Short-access.
    NcGroup grp;
    string groupname = format("Table_%2.2i",id);
    if (groupname == "Table_%2.2i") {
        groupname = format("Table_{:02d}",id);
    }
    netcdf_check(nc,grp = nc->ncid->getGroup(groupname.c_str()));
    vector<uint32_t> binning_table_read(bin_file->dim_binningtable_spat*bin_file->dim_binningtable_spec);
    netcdf_check(nc,grp.getVar("binning_table").getVar(binning_table_read.data())); // From file, that is after pre-binning.
    pixelpointer.resize(npix_unbinned); // This is the pixel pointer of the combined binning.
    npix = 0; // Will be updated to highest non-fillvalue element in the binning table plus one.
    for (size_t ibin_spat=0 ; ibin_spat<bin_file->dim_binningtable_spat ; ibin_spat++) {
        for (size_t ibin_spec=0 ; ibin_spec<bin_file->dim_binningtable_spec ; ibin_spec++) {
            uint32_t index_cur = binning_table_read[ibin_spat*bin_file->dim_binningtable_spec+ibin_spec];
            if (index_cur != NC_FILL_UINT) {
                if (index_cur+1 > npix) npix = index_cur+1; // Eventually update number of binned pixels.
            }
            for (size_t iprebin_spat=0 ; iprebin_spat<bin_file->prebinning_spat ; iprebin_spat++) {
                size_t ispat = ibin_spat*bin_file->prebinning_spat + iprebin_spat;
                for (size_t iprebin_spec=0 ; iprebin_spec<bin_file->prebinning_spec ; iprebin_spec++) {
                    size_t ispec = ibin_spec*bin_file->prebinning_spec + iprebin_spec;
                    pixelpointer[ispat*bin_file->dim_detector_spec+ispec] = index_cur; // Write this index, fillvalue or not.
                }
            }
        }
    }

    return 0;

} // }}}

