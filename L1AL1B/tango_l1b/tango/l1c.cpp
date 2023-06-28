// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "logger.h"
#include "lininv.h"
#include "vector.h"
#include "bspline.h"
#include "netcdf_object.h"
#include "planet.h"
#include "l1c.h"

// Settings functions.
Settings_l1c::Settings_l1c( // {{{
    Logger *creator
) : Settings_proc(creator)
{
    tag = "l1c";
} // }}}
int Settings_l1c::init_step( // {{{
    stringstream &stream, // A string stream to use (just initialize one).
    string &key, // Name of the setting.
    string &value, // Where the value will be stored.
    bool &recognized // Return flag whether the setting is successfully recognized.
)
{

    // Recognize specific settings.
    recognize_setting(demfile); // Detailed elevation map. Leave empty for placeholder DEM function.
    recognize_setting(semi_major_axis); // Earth property (move to database).
    recognize_setting(semi_minor_axis); // Earth property (move to database).
    recognize_setting(latitude_tol); // Tolerance for iterative calculation of latitude during geolocation (radians).
    recognize_setting(l1bfile); // Input L1B file (no L1A files are used here).
    recognize_setting(l1cfile); // Output L1C file.
    recognize_setting_vector(lat); // Latitudes of target positions.
    recognize_setting_vector(lon); // Longitudes of target positions.

    return 0;

} // }}}

// Static mathematical formula.
double getSide( // {{{
    Vector &line1,
    Vector &line2,
    Vector &tar,
    Vector &zen
)
{
    Vector leg1 = tar-line1;
    Vector leg2 = line2-line1;
    return zen.dot(leg1.cross(leg2));
} // }}}

// Constructor for the L1C structure.
L1C::L1C( // {{{
    Logger *creator,
    CKD *ckd_arg
) : Processor(creator,ckd_arg)
{
    setName("l1c");
    set = make_unique<Settings_l1c>(this);
    Processor::set = set.get();
} // }}}

int L1C::process_init( // {{{
)
{

    // Check legal settings.
    check_error(set->semi_major_axis == NC_FILL_DOUBLE,"Error: Missing semi-major axis of the planet. Setting 'semi_major_axis'");
    check_error(set->semi_minor_axis == NC_FILL_DOUBLE,"Error: Missing semi-minor axis of the planet. Setting 'semi_minor_axis'");
    check_error(set->latitude_tol == NC_FILL_DOUBLE,"Error: Missing tolerance for iterative calculation of latitude during geolocation (radians). Setting 'latitude_tol'");
    check_error(set->l1bfile.compare("") == 0,"Error: Missing L1B input file. Setting 'l1bfile'");
    check_error(set->l1cfile.compare("") == 0,"Error: Missing L1C output file. Setting 'l1cfile'");
    size_t dim_l1c_obs = set->lat.size();
    check_error(set->lon.size() != dim_l1c_obs,"Error: Number of target latitudes unequal to the number of target longitudes. Settings 'lat' and 'lon'");

    // This process has no loop over images or something like that. So
    // everything is in images. The batch is absent, because there are no
    // L1A files.

    // This process init function is the whole algorithm.

    // Create target grid, using the Earth and the selected points.
    // TODO: Make a more convenient grid generator or selector.

    // Read L1B along-track dimension.
    NetCDF_object nc_l1b(this);
    handle(nc_l1b.open(set->l1bfile,NcFile::read));

    // Orbit metadata to copy.
    string startDirection;
    string endDirection;
    long long orbit_number;
    string time_coverage_start;
    string time_coverage_stop;

    { // Search global attributes for orbit data.
        multimap<string,NcGroupAtt> mp;
        netcdf_check(&nc_l1b,mp = nc_l1b.ncid->getAtts());
        multimap<string,NcGroupAtt>::iterator search;
        search = mp.find("startDirection");
        if (search == mp.end()) startDirection = "Unknown";
        else {
            netcdf_check(&nc_l1b,search->second.getValues(startDirection));
        }
        search = mp.find("endDirection");
        if (search == mp.end()) endDirection = "Unknown";
        else {
            netcdf_check(&nc_l1b,search->second.getValues(endDirection));
        }
        search = mp.find("orbit_number");
        if (search == mp.end()) orbit_number = NC_FILL_INT64;
        else {
            netcdf_check(&nc_l1b,search->second.getValues(&orbit_number));
        }
        search = mp.find("time_coverage_start");
        if (search == mp.end()) time_coverage_start = "Unknown";
        else {
            netcdf_check(&nc_l1b,search->second.getValues(time_coverage_start));
        }
        search = mp.find("time_coverage_stop");
        if (search == mp.end()) time_coverage_stop = "Unknown";
        else {
            netcdf_check(&nc_l1b,search->second.getValues(time_coverage_stop));
        }
    }

    size_t dim_along;
    netcdf_check(&nc_l1b,dim_along = nc_l1b.ncid->getDim("bins_along_track").getSize());
    check_error(dim_along < 2,"Error: Cannot collocate if along-track dimension is less than two");
    size_t dim_vp;;
    netcdf_check(&nc_l1b,dim_vp = nc_l1b.ncid->getDim("number_of_views").getSize());
    // Copy the viewport skipping array, so that viewports that are
    // skipped for another reason can be added without having the idea
    // that we hack the settings.
    vector<int> vp_skip(dim_vp,0);
    for (size_t ivp=0 ; ivp<dim_vp && ivp<set->vp_skip.size() ; ivp++) {
        vp_skip[ivp] = set->vp_skip[ivp];
    }

    // Number of across-track positions is variable per viewport.
    size_t dim_fov; // All fields of view, of all viewports combined.
    netcdf_check(&nc_l1b,dim_fov = nc_l1b.ncid->getDim("spatial_samples_per_image").getSize());

    // Number of wavelengths.
    size_t dim_pol_wave;
    netcdf_check(&nc_l1b,dim_pol_wave = nc_l1b.ncid->getDim("polarization_bands_per_view").getSize());
    size_t dim_int_wave;
    netcdf_check(&nc_l1b,dim_int_wave = nc_l1b.ncid->getDim("intensity_bands_per_view").getSize());

    // Viewport map, to be turned into the old (CKD) way of representing it.
    // The philosophy is still that we do not need the CKD in this L1C processor.
    NcGroup grp;
    NcGroup grpw;
    netcdf_check(&nc_l1b,grp = nc_l1b.ncid->getGroup("SENSOR_VIEW_BANDS"));
    vector<uint8_t> viewportmap(dim_fov);
    netcdf_check(&nc_l1b,grp.getVar("viewport_index").getVar(viewportmap.data()));
    vector<size_t> dims_across(dim_vp,0);
    vector<size_t> ifov_start(dim_vp);
    for (size_t ifov=0 ; ifov<dim_fov ; ifov++) dims_across[viewportmap[ifov]]++;
    for (size_t ivp=0 ; ivp<dim_vp ; ivp++) {
        if (vp_skip[ivp] > 0) continue;
        if (dims_across[ivp] < 2) {
            // If there are zero FOVs, the viewport was skipped by the CKD.
            // If there is one FOV, the viewport cannot be used for L1C,
            // though its L1B product is still valid.
            if (dims_across[ivp] == 1) {
                writelog(log_warning,"Warning: Cannot collocate if across-track dimension of one, which is the case for viewport %zu. Skipping that viewport.",ivp);
            }
            vp_skip[ivp] = 1;
        }
    }
    ifov_start[0] = 0;
    
    netcdf_check(&nc_l1b,grpw = nc_l1b.ncid->getGroup("OBSERVATION_DATA"));
    
    for (size_t ivp=1 ; ivp<dim_vp ; ivp++) ifov_start[ivp] = ifov_start[ivp-1] + dims_across[ivp-1];
    // Read the wavelengths.
    vector<double> pol_wavelength(dim_vp*dim_pol_wave);
    netcdf_check(&nc_l1b,grp.getVar("polarization_wavelengths").getVar(pol_wavelength.data()));
    vector<double> int_wavelength(dim_vp*dim_int_wave);
    netcdf_check(&nc_l1b,grpw.getVar("intensity_wavelengths").getVar(int_wavelength.data()));
    // TODO: Read other variables in SENSOR_VIEW_BANDS to copy and/or combine to the L1C product.

    // Open the observation data.
    netcdf_check(&nc_l1b,grp = nc_l1b.ncid->getGroup("OBSERVATION_DATA"));
    NcVar var_l1b_i;
    netcdf_check(&nc_l1b,var_l1b_i = grp.getVar("I"));
    NcVar var_l1b_i_noise;
    netcdf_check(&nc_l1b,var_l1b_i_noise = grp.getVar("I_noise"));
    NcVar var_l1b_i_polsample;
    netcdf_check(&nc_l1b,var_l1b_i_polsample = grp.getVar("I_polsample"));
    NcVar var_l1b_i_polsample_noise;
    netcdf_check(&nc_l1b,var_l1b_i_polsample_noise = grp.getVar("I_polsample_noise"));
    NcVar var_l1b_q;
    netcdf_check(&nc_l1b,var_l1b_q = grp.getVar("q"));
    NcVar var_l1b_q_noise;
    netcdf_check(&nc_l1b,var_l1b_q_noise = grp.getVar("q_noise"));
    NcVar var_l1b_u;
    netcdf_check(&nc_l1b,var_l1b_u = grp.getVar("u"));
    NcVar var_l1b_u_noise;
    netcdf_check(&nc_l1b,var_l1b_u_noise = grp.getVar("u_noise"));
    // We need the noise on the DoLP (or AoLP) so we can reconstruct the
    // noise covariance between Q and U. This coveriance term in not in
    // the L1B, because no one undestands it.
    NcVar var_l1b_dolp_noise;
    netcdf_check(&nc_l1b,var_l1b_dolp_noise = grp.getVar("DoLP_noise"));

    // Geolocation.
    netcdf_check(&nc_l1b,grp = nc_l1b.ncid->getGroup("GEOLOCATION_DATA"));
    NcVar var_l1b_lat;
    netcdf_check(&nc_l1b,var_l1b_lat = grp.getVar("latitude"));
    NcVar var_l1b_lon;
    netcdf_check(&nc_l1b,var_l1b_lon = grp.getVar("longitude"));
    NcVar var_l1b_alt;
    netcdf_check(&nc_l1b,var_l1b_alt = grp.getVar("altitude"));
    NcVar var_l1b_sza;
    netcdf_check(&nc_l1b,var_l1b_sza = grp.getVar("solar_zenith"));
    NcVar var_l1b_saa;
    netcdf_check(&nc_l1b,var_l1b_saa = grp.getVar("solar_azimuth"));
    NcVar var_l1b_vza;
    netcdf_check(&nc_l1b,var_l1b_vza = grp.getVar("sensor_zenith"));
    NcVar var_l1b_vaa;
    netcdf_check(&nc_l1b,var_l1b_vaa = grp.getVar("sensor_azimuth"));

    // Initialize L1C output. {{{

    NetCDF_object nc_l1c(this);
    handle(nc_l1c.open(set->l1cfile,NcFile::replace));
    // Write orbit metadata.
    netcdf_check(&nc_l1c,nc_l1c.ncid->putAtt("startDirection",startDirection));
    netcdf_check(&nc_l1c,nc_l1c.ncid->putAtt("endDirection",endDirection));
    netcdf_check(&nc_l1c,nc_l1c.ncid->putAtt("orbit_number",ncInt64,orbit_number));
    netcdf_check(&nc_l1c,nc_l1c.ncid->putAtt("time_coverage_start",time_coverage_start));
    netcdf_check(&nc_l1c,nc_l1c.ncid->putAtt("time_coverage_stop",time_coverage_stop));

    // TODO (again). How to derive the target grid. The IODS proposes
    // a swath-type L1C output file. We like the nadir grid, but maybe
    // there will be a common PACE grid.
    // Now, we make a one by N grid.
    size_t dim_l1c_along = dim_l1c_obs;
    size_t dim_l1c_across = 1;
    NcDim dimid_along;
    netcdf_check(&nc_l1c,dimid_along = nc_l1c.ncid->addDim("bins_along_track",dim_l1c_along));
    NcDim dimid_across;
    netcdf_check(&nc_l1c,dimid_across = nc_l1c.ncid->addDim("bins_across_track",dim_l1c_across));
    NcDim dimid_vp;
    netcdf_check(&nc_l1c,dimid_vp = nc_l1c.ncid->addDim("number_of_views",dim_vp));
    NcDim dimid_pol_wave;
    netcdf_check(&nc_l1c,dimid_pol_wave = nc_l1c.ncid->addDim("polarization_bands_per_view",dim_pol_wave));
    NcDim dimid_int_wave;
    netcdf_check(&nc_l1c,dimid_int_wave = nc_l1c.ncid->addDim("intensity_bands_per_view",dim_int_wave));
    
    netcdf_check(&nc_l1b,grpw = nc_l1b.ncid->getGroup("OBSERVATION_DATA"));

    // Group SENSOR_VIEW_BANDS.
    netcdf_check(&nc_l1c,grp = nc_l1c.ncid->addGroup("SENSOR_VIEW_BANDS"));
    netcdf_check(&nc_l1c,grp.addVar("polarization_wavelengths",ncDouble,{dimid_vp,dimid_pol_wave}).putVar(pol_wavelength.data()));
    netcdf_check(&nc_l1c,grpw.addVar("intensity_wavelengths",ncDouble,{dimid_vp,dimid_int_wave}).putVar(int_wavelength.data()));
    // TODO: More stuff to be added.

    // Group OBSERVATION_DATA.
    netcdf_check(&nc_l1c,grp = nc_l1c.ncid->addGroup("OBSERVATION_DATA"));
    vector<NcDim> dims_pol = {dimid_along,dimid_across,dimid_vp,dimid_pol_wave};
    vector<NcDim> dims_int = {dimid_along,dimid_across,dimid_vp,dimid_int_wave};
    NcVar var_l1c_i;
    netcdf_check(&nc_l1c,var_l1c_i = grp.addVar("I",ncDouble,dims_int));
    NcVar var_l1c_i_noise;
    netcdf_check(&nc_l1c,var_l1c_i_noise = grp.addVar("I_noise",ncDouble,dims_int));
    NcVar var_l1c_i_polsample;
    netcdf_check(&nc_l1c,var_l1c_i_polsample = grp.addVar("I_polsample",ncDouble,dims_pol));
    NcVar var_l1c_i_polsample_noise;
    netcdf_check(&nc_l1c,var_l1c_i_polsample_noise = grp.addVar("I_polsample_noise",ncDouble,dims_pol));
    NcVar var_l1c_q;
    netcdf_check(&nc_l1c,var_l1c_q = grp.addVar("q",ncDouble,dims_pol));
    NcVar var_l1c_q_noise;
    netcdf_check(&nc_l1c,var_l1c_q_noise = grp.addVar("q_noise",ncDouble,dims_pol));
    NcVar var_l1c_u;
    netcdf_check(&nc_l1c,var_l1c_u = grp.addVar("u",ncDouble,dims_pol));
    NcVar var_l1c_u_noise;
    netcdf_check(&nc_l1c,var_l1c_u_noise = grp.addVar("u_noise",ncDouble,dims_pol));
    NcVar var_l1c_dolp;
    netcdf_check(&nc_l1c,var_l1c_dolp = grp.addVar("DoLP",ncDouble,dims_pol));
    NcVar var_l1c_dolp_noise;
    netcdf_check(&nc_l1c,var_l1c_dolp_noise = grp.addVar("DoLP_noise",ncDouble,dims_pol));
    NcVar var_l1c_aolp;
    netcdf_check(&nc_l1c,var_l1c_aolp = grp.addVar("AoLP",ncDouble,dims_pol));
    NcVar var_l1c_aolp_noise;
    netcdf_check(&nc_l1c,var_l1c_aolp_noise = grp.addVar("AoLP_noise",ncDouble,dims_pol));

    // Group GEOLOCATION_DATA.
    netcdf_check(&nc_l1c,grp = nc_l1c.ncid->addGroup("GEOLOCATION_DATA"));
    vector<NcDim> dims_geo = {dimid_along,dimid_across};
    netcdf_check(&nc_l1c,grp.addVar("latitude",ncDouble,dims_geo).putVar(set->lat.data())); // Still in degrees.
    netcdf_check(&nc_l1c,grp.addVar("longitude",ncDouble,dims_geo).putVar(set->lon.data())); // Still in degrees.
    NcVar var_l1c_alt;
    netcdf_check(&nc_l1c,var_l1c_alt = grp.addVar("altitude",ncDouble,dims_geo));
    vector<NcDim> dims_ang = {dimid_along,dimid_across,dimid_vp}; // I think that the solar angles also depend on viewport, because of time differences.
    NcVar var_l1c_sza;
    netcdf_check(&nc_l1c,var_l1c_sza = grp.addVar("solar_zenith",ncDouble,dims_ang));
    NcVar var_l1c_saa;
    netcdf_check(&nc_l1c,var_l1c_saa = grp.addVar("solar_azimuth",ncDouble,dims_ang));
    NcVar var_l1c_vza;
    netcdf_check(&nc_l1c,var_l1c_vza = grp.addVar("sensor_zenith",ncDouble,dims_ang));
    NcVar var_l1c_vaa;
    netcdf_check(&nc_l1c,var_l1c_vaa = grp.addVar("sensor_azimuth",ncDouble,dims_ang));

    // }}}

    // Construct the Earth. We will have latitudes and longitudes, and we need
    // cartesian coordinates.
    Planet earth(this,set->semi_major_axis,set->semi_minor_axis,set->latitude_tol);
    if (set->demfile.compare("") != 0) {
        handle(earth.include_dem(set->demfile))
    }

    // Construct array of target positions. For target positions, altitudes are
    // not relevant, because the final projection is along the zenith.
    vector<Vector> target(dim_l1c_obs);
    vector<Vector> tarzen(dim_l1c_obs); // The zenith of the target.
    vector<double> alt(dim_l1c_obs);
    for (size_t iobs=0 ; iobs<dim_l1c_obs ; iobs++) {
        double lat = set->lat[iobs]*DEGREES;
        double lon = set->lon[iobs]*DEGREES;
        handle(earth.dem(lat,lon,alt[iobs]));
        earth.xyz(lat,lon,alt[iobs],&target[iobs]);
        tarzen[iobs][0] = cos(lat)*cos(lon);
        tarzen[iobs][1] = cos(lat)*sin(lon);
        tarzen[iobs][2] = sin(lat);
    }
    netcdf_check(&nc_l1c,var_l1c_alt.putVar(alt.data()));

    // From now on, we will have an along-track dimensions, which is fixed
    // for all viewports, and an across-track dimension, which may differ
    // for different viewport. Optimal would be to keep all the along
    // track things outside the for loop, but that is an optimization step
    // that may be done in the future and only if it saves significant time.
    // Keeping all this inside the loop keeps the related code together.

    // Loop over viewports (source grids).
    for (size_t ivp=0 ; ivp<dim_vp ; ivp++) {

        // Do not execute skipped viewports.
        if (vp_skip[ivp] > 0) continue;

        writelog(log_trace,"Constructing grid for viewport %zu.",ivp);

        size_t &dim_across = dims_across[ivp];
        // Interpret some dimension sizes.
        // Go to 'int' type, because we may virtually move to indices like -1 and that
        // should be recognized properly.
        int ncell_along = (int)dim_along-1;
        int ncell_across = (int)dim_across-1;
        int ncell = ncell_along*ncell_across;

        // Construct 2D B-splines over an index grid. The index grid is the same for
        // all the viewports. TODO: No. On the contrary, no.
        vector<double> gain_along((dim_along+2)*(dim_along+2));
        vector<double> gain_across((dim_across+2)*(dim_across+2));
        for (size_t al_ac=0 ; al_ac<2 ; al_ac++) {
            size_t dim = al_ac==0?dim_along:dim_across;
            double *gain = al_ac==0?gain_along.data():gain_across.data();
            // Start with the inverse of the corresponding Toeplitz matrix.
            // Initialize the Chebyshev polynomials, without recycling between along-track
            // and across track. The x-value of these is always two.
            // The Chebyshev polynomials become huge, so we save factors
            // between adjacent orders. chebfacts[i] = cheb[i+1] / cheb[i],
            // usually close to four.
            vector<double> chebfacts(dim+2); // Negative defition to get rid of the sign flip action.
            chebfacts[0] = -4.0;
            for (size_t icheb=1 ; icheb<dim+2 ; icheb++) {
                chebfacts[icheb] = -4.0 - 1.0/chebfacts[icheb-1];
            }
            vector<double> els_cur((dim+3)/2);
            els_cur[0] = -6.0/chebfacts[dim+1];
            gain[0] = els_cur[0];
            gain[(dim+1)*(dim+2)+dim+1] = els_cur[0];
            bool firsthalf = true;
            for (size_t idim1=1 ; idim1<dim+2 ; idim1++) {
                size_t sz2 = firsthalf?idim1:dim+2-idim1;
                for (size_t idim2=0 ; idim2<sz2 ; idim2++) els_cur[idim2] /= chebfacts[dim+1-idim1];
                if (firsthalf) {
                    sz2++;
                    els_cur[sz2-1] = els_cur[sz2-2] * chebfacts[sz2-2];
                }
                for (size_t idim2=0 ; idim2<sz2 ; idim2++) {
                    gain[idim1*(dim+2)+idim2] = els_cur[idim2];
                    // Possibly writing the same index again, but ifs
                    // may be more expensive than additional write actions.
                    gain[idim2*(dim+2)+idim1] = els_cur[idim2];
                    gain[(dim+1-idim1)*(dim+2)+dim+1-idim2] = els_cur[idim2];
                    gain[(dim+1-idim2)*(dim+2)+dim+1-idim1] = els_cur[idim2];
                }
                if (2*sz2 > dim+1) firsthalf = false;
            }
            // Apply perturbations.
            vector<size_t> ifirst = {0,0,0,dim+1,dim+1,dim+1};
            vector<size_t> isecond = {0,1,2,dim-1,dim,dim+1};
            vector<double> pert = {1.0/3.0,-13.0/6.0,1.0,1.0,-13.0/6.0,1.0/3.0};
            for (size_t icase=0 ; icase<6 ; icase++) {
                vector<double> one(dim+2);
                vector<double> two(dim+2);
                for (size_t idim=0 ; idim<dim+2 ; idim++) {
                    one[idim] = gain[idim*(dim+2)+ifirst[icase]];
                    two[idim] = gain[isecond[icase]*(dim+2)+idim];
                }
                double cross = gain[isecond[icase]*(dim+2)+ifirst[icase]];
                for (size_t idim_first=0 ; idim_first<dim+2 ; idim_first++) {
                    for (size_t idim_second=0 ; idim_second<dim+2 ; idim_second++) {
                        gain[idim_first*(dim+2)+idim_second] -= one[idim_first]*two[idim_second]*pert[icase] / (1.0 + cross*pert[icase]);
                    }
                }
            }
        }

        // Read latitudes and longitudes for this view port.
        vector<double> l1b_lat(dim_along*dim_across);
        vector<double> l1b_lon(dim_along*dim_across);
        vector<double> l1b_alt(dim_along*dim_across);
        const vector<size_t> l1b_start = {0,ifov_start[ivp]};
        const vector<size_t> l1b_cnt = {dim_along,dim_across};
        netcdf_check(&nc_l1b,var_l1b_lat.getVar(l1b_start,l1b_cnt,l1b_lat.data()));
        netcdf_check(&nc_l1b,var_l1b_lon.getVar(l1b_start,l1b_cnt,l1b_lon.data()));
        netcdf_check(&nc_l1b,var_l1b_alt.getVar(l1b_start,l1b_cnt,l1b_alt.data()));

        // Construct the grid in cartesian coordinates.
        vector<Vector> grid(dim_along*dim_across);

        for (size_t igrid=0 ; igrid<dim_along*dim_across ; igrid++) {
            double lat = l1b_lat[igrid]*DEGREES;
            double lon = l1b_lon[igrid]*DEGREES;
            earth.xyz(lat,lon,l1b_alt[igrid],&grid[igrid]);
        }

        // Cell zeniths (for mesh investigation). They need not be precise, only such that
        // the cells do not become complex after projection.
        vector<Vector> zens(ncell);
        for (int icell_along=0 ; icell_along<ncell_along ; icell_along++) {
            for (int icell_across=0 ; icell_across<ncell_across ; icell_across++) {
                Vector &ll = grid[icell_along*dim_across+icell_across];
                Vector &lr = grid[icell_along*dim_across+icell_across+1];
                Vector &ul = grid[(icell_along+1)*dim_across+icell_across];
                Vector &ur = grid[(icell_along+1)*dim_across+icell_across+1];
                Vector avg = 0.25*(ll+lr+ul+ur);
                double lat;
                double lon;
                double alt;
                earth.lla(avg,&lat,&lon,&alt);
                zens[icell_along*ncell_across+icell_across] = Vector(cos(lat)*cos(lon),cos(lat)*sin(lon),sin(lat));
            }
        }
        // Alternatively, use corner grid and center cells for zeniths.

        // Calibration for what is left or right.
        // We assume that if along-flight is forward and the zenith is
        // upward, across-track is to the right. This must be verified.
        // We use the zeniths of the cells. We may have used the zeniths
        // of the points that are the target points here (ATBD), but the cell
        // should also suffice.
        double cal = 0.0; // Neutral.
        for (int icell_along=0 ; icell_along<ncell_along ; icell_along++) {
            for (int icell_across=0 ; icell_across<ncell_across ; icell_across++) {
                Vector *base = &grid[icell_along*dim_across+icell_across];
                Vector &z = zens[icell_along*ncell_across+icell_across];
                size_t nright = 0;
                if (getSide(base[0],base[dim_across],base[dim_across+1],z) > 0.0) nright++;
                if (getSide(base[dim_across],base[dim_across+1],base[1],z) > 0.0) nright++;
                if (getSide(base[dim_across+1],base[1],base[0],z) > 0.0) nright++;
                if (getSide(base[1],base[0],base[dim_across],z) > 0.0) nright++;
                check_error(nright == 2,"Error: Complex quadrilateral in grid.");
                if (cal == 0.0) cal = nright>2?1.0:-1.0;
                check_error((cal > 0.0) != (nright>2),"Error: Clockwise and counter-clockwise cells co-exist.")
                if (nright == 1 || nright == 3) {
                    writelog(log_warning,"Warning: Concave quadrilateral in grid. This is rare, though it is supported.");
                }
            }
        }

        writelog(log_trace,"Collocating observations into grid for viewport %zu.",ivp);

        // The observation loop calculates a factor matrix of interpolation factors
        // from the L1B along-across grid to the observation index on L1C. If it
        // fails, the failure flag is turned on.
        vector<double> factormat(dim_l1c_obs*dim_along*dim_across,0.0);
        double *factormat_cur = factormat.data();
        vector<bool> fail(dim_l1c_obs,false);

        for (size_t iobs=0 ; iobs<dim_l1c_obs ; iobs++) {

            // Target point and its zenith, used many times.
            Vector &tar = target[iobs];
            Vector &zen = tarzen[iobs];

            // Set a-priori cell where the dot product between cell zenith
            // and target zenith is maximum. We could have taken the middle
            // of the cell and the target, but we forgot the middle of the
            // cells and we still have the zeniths for verifying the mesh.
            // The a-priori need not be precise, so it does not matter.

            int icell = 0;
            double zenith_dot_product = zen.dot(zens[0]);
            for (int icell_search=1 ; icell_search<ncell ; icell_search++) {
                double newdot = zen.dot(zens[icell_search]);
                if (newdot > zenith_dot_product) {
                    zenith_dot_product = newdot;
                    icell = icell_search;
                }
            }
            // Split up a-priori cell to along and across track coordinates.
            int icell_along = icell / ncell_across; // Integer division.
            int icell_across = icell - icell_along*ncell_across; // Or icell modulo ncell_across.

            // getSide constructions.
            double construction_left;
            double construction_right;
            double construction_down;
            double construction_up;

            // Execute search algorithm.
            // Rounding errors cannot cause the search algorithm to go
            // back and forth, in which case we will compare the badnesses,
            // which should approach zero. This can occur if the target
            // point is on a grid line. In such a case, both cells are
            // correct, but the lowest zero badness will be taken.

            // Concave quadrilaterals are handled such that, if a corner
            // is also outside the line, it is a concave quadrilateral
            // and that direction is flagged as a concave direction. Only
            // if there are desires to get out in two concave directions,
            // there is a wish to get out, in the second-worst direction.
            // If there is a concave quadrilateral and we are on a rim or
            // a vertex, we just save the badness and end at the least bad
            // cell, which should have a badness equal to a rounding error.

            // We have to save the concave (overruled) movement desires.
            vector<double> badness(ncell,0.0);
            vector<bool> visited(ncell,false);
            vector<bool> concave_left(ncell,false);
            vector<bool> concave_right(ncell,false);
            vector<bool> concave_down(ncell,false);
            vector<bool> concave_up(ncell,false);
            bool continue_investigation = false; // Flag for cycle due to rounding errors.
            for (int iter=0 ; iter<ncell+1 ; iter++) { // Maximum number of loop iterations is smaller than this.
                // Protection against infinite loop.
                check_error(iter == ncell,"Program error: Too many iterations in search loop. Infinite loop in search algorithm or zero-sized grid.");
                icell = icell_along*ncell_across + icell_across;
                if (visited[icell]) {
                    // Loop detected.
                    // For the final step, these must be rounding errors and
                    // the candidates will be compared.
                    continue_investigation = true;
                    break;
                }
                visited[icell] = true;

                // Calculate movement desires.
                // Index offset of the cell in the grid.
                Vector *base = &grid[icell_along*dim_across+icell_across];
                // Left.
                construction_left = getSide(base[0],base[dim_across],tar,zen)*cal;
                bool wish_left = construction_left < 0.0;
                // Right.
                construction_right = getSide(base[dim_across+1],base[1],tar,zen)*cal;
                bool wish_right = construction_right < 0.0;
                // Down.
                construction_down = getSide(base[1],base[0],tar,zen)*cal;
                bool wish_down = construction_down < 0.0;
                // Up.
                construction_up = getSide(base[dim_across],base[dim_across+1],tar,zen)*cal;
                bool wish_up = construction_up < 0.0;

                // Check for concave angles. In such a case, it is possible
                // that a desire to move out is calculated while the
                // target point is in the cell. More specifically, if the
                // reflex-angle condition is fulfilled and this is the
                // only desire, the point is in the cell. We will not
                // remove the wish flag but just add a concave flag under
                // the condition. Flags are saved in arrays, because
                // for numeric rounding-error loops, we may finally
                // converge to a cell that is not the last one visited.
                if (wish_left) {
                    double construction_corner1 = getSide(base[0],base[dim_across],base[1],zen)*cal;
                    double construction_corner2 = getSide(base[0],base[dim_across],base[dim_across+1],zen)*cal;
                    if (construction_corner1 < 0.0 || construction_corner2 < 0.0) concave_left[icell] = true;
                }
                if (wish_right) {
                    double construction_corner1 = getSide(base[dim_across+1],base[1],base[dim_across],zen)*cal;
                    double construction_corner2 = getSide(base[dim_across+1],base[1],base[0],zen)*cal;
                    if (construction_corner1 < 0.0 || construction_corner2 < 0.0) concave_right[icell] = true;
                }
                if (wish_down) {
                    double construction_corner1 = getSide(base[1],base[0],base[dim_across+1],zen)*cal;
                    double construction_corner2 = getSide(base[1],base[0],base[dim_across],zen)*cal;
                    if (construction_corner1 < 0.0 || construction_corner2 < 0.0) concave_down[icell] = true;
                }
                if (wish_up) {
                    double construction_corner1 = getSide(base[dim_across],base[dim_across+1],base[0],zen)*cal;
                    double construction_corner2 = getSide(base[dim_across],base[dim_across+1],base[1],zen)*cal;
                    if (construction_corner1 < 0.0 || construction_corner2 < 0.0) concave_up[icell] = true;
                }

                size_t nwish = 0;
                if (wish_left) nwish++;
                if (wish_right) nwish++;
                if (wish_down) nwish++;
                if (wish_up) nwish++;
                if (nwish == 0) {
                    // Perfect conversion.
                    break;
                }
                if (nwish == 1) {
                    // Wish to move in one side. If this movement is under discussion
                    // because of a reflex angle, we know we are in the cell and we are
                    // ready. Otherwise, move in that direction.
                    if (wish_left && concave_left[icell]) break; // Ready.
                    if (wish_right && concave_right[icell]) break; // Ready.
                    if (wish_down && concave_down[icell]) break; // Ready.
                    if (wish_up && concave_up[icell]) break; // Ready.
                    // Not converged, thus badness is written for if it is a rounding error.
                    if (wish_left) badness[icell] = construction_left;
                    if (wish_right) badness[icell] = construction_right;
                    if (wish_down) badness[icell] = construction_down;
                    if (wish_up) badness[icell] = construction_up;

                    // Move in that direction.
                    if (wish_left) icell_across--;
                    if (wish_right) icell_across++;
                    if (wish_down) icell_along--;
                    if (wish_up) icell_along++;
                }
                if (nwish == 2) {
                    // If these are opposite, it is an error unless one side is concave. Then,
                    // go to the other side.
                    if (wish_up && wish_down) {
                        if (concave_down[icell] == concave_up[icell]) {
                            writelog(log_warning,"Error: Complex quadrilateral encountered during search.\n");
                            fail[iobs] = true;
                            break;
                        }
                        if (concave_down[icell]) icell_along++;
                        else icell_along--;
                        continue;
                    }
                    if (wish_left && wish_right) {
                        if (concave_left[icell] == concave_right[icell]) {
                            writelog(log_warning,"Error: Complex quadrilateral encountered during search.\n");
                            fail[iobs] = true;
                            break;
                        }
                        if (concave_left[icell]) icell_across++;
                        else icell_across--;
                        continue;
                    }
                    // Now, we assume there is a desire to move in two perpendicular directions.
                    // Objection flags.
                    bool object_across = false;
                    bool object_along = false;
                    // In theory, it is impossible to have a desire to move in two
                    // perpendicular and one is concave. In such a case, the other
                    // concave border should also trigger. Here, however, it is
                    // most important to save badness for eventual rounding-error
                    // driven infinite loop and stay resillient against numeric
                    // errors anyway.
                    // Concave wish always has lowest priority.
                    object_across = (concave_left[icell] || concave_right[icell]);
                    object_along = (concave_down[icell] || concave_up[icell]);
                    // Save badness. If there is one objection, save the other.
                    // If there are two objections, save the least bad one.
                    // If there are no objections, save the worst one.
                    double badness_across = wish_left?construction_left:construction_right; // Flag wish_left equals !wish_right here.
                    double badness_along = wish_down?construction_down:construction_up; // Flag wish_down equals !wish_up here.
                    if (object_across) {
                        if (object_along) {
                            // Save least bad one. Low (strongly negative) is bad.
                            // Also, overrule concaveness for the least bad badness. Only the
                            // worst one can be concave. This is important if this ends up
                            // being the correct cell and that the one of the two movement
                            // desires is a rounding error and the other one in concave.
                            // Then, only the true concave direction should be flagged as concave
                            // for the interpolation part.
                            if (badness_across < badness_along) {
                                // Across are concave, along are not.
                                badness[icell] = badness_along;
                                concave_down[icell] = false;
                                concave_up[icell] = false;
                            } else {
                                // Along are concave, across are not.
                                badness[icell] = badness_across;
                                concave_left[icell] = false;
                                concave_right[icell] = false;
                            }
                        } else badness[icell] = badness_along;
                    } else {
                        if (object_along) badness[icell] = badness_across;
                        else {
                            // Save worst one. Low (strongly negative) is bad.
                            if (badness_across < badness_along) badness[icell] = badness_across;
                            else badness[icell] = badness_along;
                        }
                    }
                    // Second rule.
                    if (object_across == object_along) {
                        // Prefer moving in along-track direction, except when this involves wrapping or
                        // revisiting.
                        // If across-track movement also involves wrapping or revisiting, execute the wrap or
                        // revisit in along-track direction.

                        // Move virtually.
                        int icell_along_desire = icell_along;
                        int icell_across_desire = icell_across;
                        if (wish_left) icell_across_desire--;
                        if (wish_right) icell_across_desire++;
                        if (wish_down) icell_along_desire--;
                        if (wish_up) icell_along_desire++;

                        // Object if you move out of the grid or to a visited cell.
                        object_across = icell_across_desire == -1 || icell_across_desire == ncell_across;
                        if (!object_across) object_across = visited[icell_along*ncell_across+icell_across_desire]; // If-clause to prevent any chance of out-of-bound errors in reading visited.
                        object_along = icell_along_desire == -1 || icell_along_desire == ncell_along;
                        if (!object_along) object_along = visited[icell_along_desire*ncell_across+icell_across]; // If-clause to prevent any chance of out-of-bound errors in reading visited.
                    }
                    // These are the final objection flags. If it is still a draw, take
                    // the along-track direction.
                    if (object_along && !object_across) {
                        if (wish_left) icell_across--;
                        if (wish_right) icell_across++;
                    } else {
                        if (wish_down) icell_along--;
                        if (wish_up) icell_along++;
                    }
                }
                if (nwish == 3) {
                    // Do opposite the one that is not there, even if this involves wrapping.
                    if (!wish_left) icell_across++;
                    if (!wish_right) icell_across--;
                    if (!wish_down) icell_along++;
                    if (!wish_up) icell_along--;
                    // Do not write badness, it will have a badness of zero, which, in the
                    // interpretation will be flagged as infinitely bad. It is not possible
                    // to have three rounding errors unless the grid is ill.
                }
                if (nwish == 4) {
                    writelog(log_warning,"Warning: Current cell is on the opposite side as target. This should be impossible or nearly impossible, although it should not harm the algorithm.");
                    // This is rare. You are exactly on the opposite of the globe. Do somehting,
                    // it does not matter if you wrap, but do not die if not needed.
                    // We can assume that you did not wrap to this position, so moving in
                    // any direction along track should be safe.
                    icell_along++;
                    // Do not write badness, it will have a badness of zero, which, in the
                    // interpretation will be flagged as infinitely bad. It is not possible
                    // to have four rounding errors unless the grid is ill.
                }

                // Detect moving out of the orbit.
                if (icell_across == -1 || icell_across == ncell_across || icell_along == -1 || icell_along == ncell_along) {
                    fail[iobs] = true;
                    break;
                }

            }

            // If the continue investigation flag is turned on, it is
            // a loop by rounding erors and the badness array is consulted.
            // Terribly bad cells will not have a badness written, hence
            // their array value is at zero. These are indeed interpreted as
            // very bad and will not be considered. For the rest, the best
            // cell is the one with the smallest negative badness (highest
            // numeric value of badness).
            if (continue_investigation) {
                double best_badness = 0.0;
                int icell_best = -1; // Nothing found.
                for (int icell_loop=0 ; icell_loop<ncell ; icell_loop++) {
                    if (badness[icell_loop] == 0.0) continue; // Not visited or flagged as chanceless (badness is zero).
                    if (best_badness == 0.0 || badness[icell_loop] > best_badness) {
                        icell_best = icell_loop;
                        best_badness = badness[icell_loop];
                    }
                }
                check_error(best_badness == 0.0,"Error: Continue investigation only involved cells flagged as chanceless.");
                // Set 2D and 1D indices accordingly.
                icell = icell_best;
                icell_along = icell/ncell_across;
                icell_across = icell - icell_along*ncell_across;
            }

            // Failed collocations are ready. No floating indices need to be calculated.
            if (!fail[iobs]) {

                check_error(icell < 0 || icell >= ncell,"Program error: Variable icell not in boundaries after search. It is a miracle that no worse error occured before.");

                // Find the weight factors so that the floating part of the index
                // can be retrieved.

                // Start iterative interpolation.
                double weightleft = 0.5;
                double weightright = 0.5;
                double weightdown = 0.5;
                double weightup = 0.5;
                double mv = 0.25;
                Vector &lowleft = grid[icell_along*dim_across+icell_across];
                Vector &lowright = grid[icell_along*dim_across+icell_across+1];
                Vector &upleft = grid[(icell_along+1)*dim_across+icell_across];
                Vector &upright = grid[(icell_along+1)*dim_across+icell_across+1];
                for (size_t i=0 ; i<200 ; i++) { // Iterative loop. TODO: Set loop limit.
                    Vector left = weightdown*lowleft + weightup*upleft;
                    Vector right = weightdown*lowright + weightup*upright;
                    Vector down = weightleft*lowleft + weightright*lowright;
                    Vector up = weightleft*upleft + weightright*upright;
                    bool force_left = false;
                    bool force_right = false;
                    if (concave_left[icell]) {
                        double construction_corner1 = getSide(down,up,lowright,zen)*cal;
                        double construction_corner2 = getSide(down,up,upright,zen)*cal;
                        if (construction_corner1 < 0.0 || construction_corner2 < 0.0) force_right = true;
                    }
                    if (concave_right[icell]) {
                        double construction_corner1 = getSide(down,up,lowleft,zen)*cal;
                        double construction_corner2 = getSide(down,up,upleft,zen)*cal;
                        if (construction_corner1 > 0.0 || construction_corner2 > 0.0) force_left = true;
                    }
                    // Only if nothing is forced here, do normal.
                    bool go_left;
                    if (force_left) go_left = true;
                    else if (force_right) go_left = false;
                    else {
                        double construction = getSide(down,up,tar,zen)*cal;
                        go_left = construction < 0.0;
                    }
                    // Repeat joke for down and up.
                    bool force_down = false;
                    bool force_up = false;
                    if (concave_down[icell]) {
                        double construction_corner1 = getSide(left,right,upleft,zen)*cal;
                        double construction_corner2 = getSide(left,right,upright,zen)*cal;
                        if (construction_corner1 > 0.0 || construction_corner2 > 0.0) force_up = true;
                    }
                    if (concave_up[icell]) {
                        double construction_corner1 = getSide(left,right,lowleft,zen)*cal;
                        double construction_corner2 = getSide(left,right,lowright,zen)*cal;
                        if (construction_corner1 < 0.0 || construction_corner2 < 0.0) force_down = true;
                    }
                    // Only if nothing is forced here, do normal.
                    bool go_up;
                    if (force_down) go_up = false;
                    else if (force_up) go_up = true;
                    else {
                        double construction = getSide(left,right,tar,zen)*cal;
                        go_up = construction < 0.0;
                    }
                    // Execute movement.
                    if (go_left) {
                        weightleft += mv;
                        weightright -= mv;
                    } else {
                        weightleft -= mv;
                        weightright += mv;
                    }
                    if (go_up) {
                        weightdown -= mv;
                        weightup += mv;
                    } else {
                        weightdown += mv;
                        weightup -= mv;
                    }
                    mv /= 2.0;
                }

                vector<double> splinecoefs_along(4);
                vector<double> splinecoefs_across(4);
                // Zeroth element: decreasing part of the low section.
                splinecoefs_along[0] = 1.0/6.0 * pow(weightdown,3.0);
                splinecoefs_across[0] = 1.0/6.0 * pow(weightleft,3.0);
                // First element: decreasing part of the high section.
                splinecoefs_along[1] = -1.0/2.0 * pow(weightdown,3.0) + 1.0/2.0 * pow(weightdown,2.0) + 1.0/2.0 * weightdown + 1.0/6.0;
                splinecoefs_across[1] = -1.0/2.0 * pow(weightleft,3.0) + 1.0/2.0 * pow(weightleft,2.0) + 1.0/2.0 * weightleft + 1.0/6.0;
                // Second element: increasing part of high section.
                splinecoefs_along[2] = -1.0/2.0 * pow(weightup,3.0) + 1.0/2.0 * pow(weightup,2.0) + 1.0/2.0 * weightup + 1.0/6.0;
                splinecoefs_across[2] = -1.0/2.0 * pow(weightright,3.0) + 1.0/2.0 * pow(weightright,2.0) + 1.0/2.0 * weightright + 1.0/6.0;
                // Third element: increasing part of low section.
                splinecoefs_along[3] = 1.0/6.0 * pow(weightup,3.0);
                splinecoefs_across[3] = 1.0/6.0 * pow(weightright,3.0);

                for (size_t isect_along=0 ; isect_along<4 ; isect_along++) {
                    for (size_t isect_across=0 ; isect_across<4 ; isect_across++) {
                        for (size_t ialong=0 ; ialong<dim_along ; ialong++) {
                            for (size_t iacross=0 ; iacross<dim_across ; iacross++) {
                                factormat_cur[ialong*dim_across+iacross] += splinecoefs_along[isect_along]*splinecoefs_across[isect_across] * gain_along[(icell_along+isect_along)*(dim_along+2)+ialong+1] * gain_across[(icell_across+isect_across)*(dim_across+2)+iacross+1];
                            }
                        }
                    }
                }
            }

            // Move running pointer.
            factormat_cur += dim_along*dim_across;

        } // Loop over observations.

        // Apply the factor matrix to everything that is being collocated.

        writelog(log_trace,"Interpolating intensity product.");
        // Radiance.
        {
            vector<double> l1b_radiance(dim_along*dim_across*dim_int_wave);
            vector<double> l1b_radiance_noise(dim_along*dim_across*dim_int_wave);
            // Slice out one viewport.
            vector<size_t> l1b_start = {0,ifov_start[ivp],0};
            vector<size_t> l1b_cnt = {dim_along,dim_across,dim_int_wave};
            netcdf_check(&nc_l1b,var_l1b_i.getVar(l1b_start,l1b_cnt,l1b_radiance.data()));
            netcdf_check(&nc_l1b,var_l1b_i_noise.getVar(l1b_start,l1b_cnt,l1b_radiance_noise.data()));
            factormat_cur = factormat.data();
            for (size_t iobs=0 ; iobs<dim_l1c_obs ; iobs++) {

                // Where to write.
                vector<size_t> l1c_start = {iobs,0,ivp,0};
                vector<size_t> l1c_cnt = {1,1,1,dim_int_wave};

                if (fail[iobs]) {
                    // Write fill values.
                    vector<double> fillvalues(dim_int_wave,NC_FILL_DOUBLE);
                    netcdf_check(&nc_l1c,var_l1c_i.putVar(l1c_start,l1c_cnt,fillvalues.data()));
                    netcdf_check(&nc_l1c,var_l1c_i_noise.putVar(l1c_start,l1c_cnt,fillvalues.data()));
                } else {
                    vector<double> l1c_radiance(dim_int_wave,0.0);
                    vector<double> l1c_radiance_noise(dim_int_wave,0.0);
                    double *l1b_radiance_cur = l1b_radiance.data();
                    double *l1b_radiance_noise_cur = l1b_radiance_noise.data();
                    for (size_t ialong=0 ; ialong<dim_along ; ialong++) {
                        for (size_t iacross=0 ; iacross<dim_across ; iacross++) {
                            for (size_t iwave=0 ; iwave<dim_int_wave ; iwave++) {
                                l1c_radiance[iwave] += factormat_cur[ialong*dim_across+iacross] * l1b_radiance_cur[iwave];
                                l1c_radiance_noise[iwave] += pow(factormat_cur[ialong*dim_across+iacross] * l1b_radiance_noise_cur[iwave],2.0); // L1C radiance noise is here L1C radiance variance.
                            }
                            l1b_radiance_cur += dim_int_wave;
                            l1b_radiance_noise_cur += dim_int_wave;
                        }
                    }
                    // Convert variance to noise.
                    for (size_t iwave=0 ; iwave<dim_int_wave ; iwave++) l1c_radiance_noise[iwave] = sqrt(l1c_radiance_noise[iwave]);

                    // Write the results.
                    netcdf_check(&nc_l1c,var_l1c_i.putVar(l1c_start,l1c_cnt,l1c_radiance.data()));
                    netcdf_check(&nc_l1c,var_l1c_i_noise.putVar(l1c_start,l1c_cnt,l1c_radiance_noise.data()));
                }
                factormat_cur += dim_along*dim_across;
            }
        }

        writelog(log_trace,"Interpolating polarization product.");
        // Polarization. TODO: Discuss on what to collocate, small or large Q/U?
        // TODO: Discuss averaging terhnique with or without averaging out noise.
        {
            vector<double> l1b_i_polsample(dim_along*dim_across*dim_pol_wave);
            vector<double> l1b_i_polsample_noise(dim_along*dim_across*dim_pol_wave);
            vector<double> l1b_q(dim_along*dim_across*dim_pol_wave);
            vector<double> l1b_q_noise(dim_along*dim_across*dim_pol_wave);
            vector<double> l1b_u(dim_along*dim_across*dim_pol_wave);
            vector<double> l1b_u_noise(dim_along*dim_across*dim_pol_wave);
            vector<double> l1b_dolp_noise(dim_along*dim_across*dim_pol_wave);
            vector<double> l1b_qu_noise_covar(dim_along*dim_across*dim_pol_wave);
            // Slice out one viewport.
            vector<size_t> l1b_start = {0,ifov_start[ivp],0};
            vector<size_t> l1b_cnt = {dim_along,dim_across,dim_pol_wave};
            netcdf_check(&nc_l1b,var_l1b_i_polsample.getVar(l1b_start,l1b_cnt,l1b_i_polsample.data()))
            netcdf_check(&nc_l1b,var_l1b_i_polsample_noise.getVar(l1b_start,l1b_cnt,l1b_i_polsample_noise.data()))
            netcdf_check(&nc_l1b,var_l1b_q.getVar(l1b_start,l1b_cnt,l1b_q.data()))
            netcdf_check(&nc_l1b,var_l1b_q_noise.getVar(l1b_start,l1b_cnt,l1b_q_noise.data()))
            netcdf_check(&nc_l1b,var_l1b_u.getVar(l1b_start,l1b_cnt,l1b_u.data()))
            netcdf_check(&nc_l1b,var_l1b_u_noise.getVar(l1b_start,l1b_cnt,l1b_u_noise.data()))
            netcdf_check(&nc_l1b,var_l1b_dolp_noise.getVar(l1b_start,l1b_cnt,l1b_dolp_noise.data()));
            // Acquire noise covariance from all noise terms.
            for (size_t imat=0 ; imat<dim_along*dim_across*dim_pol_wave ; imat++) {
                l1b_qu_noise_covar[imat] = (pow(l1b_dolp_noise[imat],2.0) * (pow(l1b_q[imat],2.0) + pow(l1b_u[imat],2.0)) - pow(l1b_q[imat]*l1b_q_noise[imat],2.0) - pow(l1b_u[imat]*l1b_u_noise[imat],2.0)) / (l1b_q[imat]*l1b_u[imat]);
            }
            factormat_cur = factormat.data();
            for (size_t iobs=0 ; iobs<dim_l1c_obs ; iobs++) {

                // Where to write.
                vector<size_t> l1c_start = {iobs,0,ivp,0};
                vector<size_t> l1c_cnt = {1,1,1,dim_pol_wave};

                if (fail[iobs]) {

                    // Write fill values.
                    vector<double> fillvalues(dim_pol_wave,NC_FILL_DOUBLE);
                    netcdf_check(&nc_l1c,var_l1c_i_polsample.putVar(l1c_start,l1c_cnt,fillvalues.data()));
                    netcdf_check(&nc_l1c,var_l1c_i_polsample_noise.putVar(l1c_start,l1c_cnt,fillvalues.data()));
                    netcdf_check(&nc_l1c,var_l1c_q.putVar(l1c_start,l1c_cnt,fillvalues.data()));
                    netcdf_check(&nc_l1c,var_l1c_q_noise.putVar(l1c_start,l1c_cnt,fillvalues.data()));
                    netcdf_check(&nc_l1c,var_l1c_u.putVar(l1c_start,l1c_cnt,fillvalues.data()));
                    netcdf_check(&nc_l1c,var_l1c_u_noise.putVar(l1c_start,l1c_cnt,fillvalues.data()));
                    netcdf_check(&nc_l1c,var_l1c_dolp.putVar(l1c_start,l1c_cnt,fillvalues.data()));
                    netcdf_check(&nc_l1c,var_l1c_dolp_noise.putVar(l1c_start,l1c_cnt,fillvalues.data()));
                    netcdf_check(&nc_l1c,var_l1c_aolp.putVar(l1c_start,l1c_cnt,fillvalues.data()));
                    netcdf_check(&nc_l1c,var_l1c_aolp_noise.putVar(l1c_start,l1c_cnt,fillvalues.data()));

                } else {

                    vector<double> l1c_i_polsample(dim_pol_wave,0.0);
                    vector<double> l1c_i_polsample_noise(dim_pol_wave,0.0);
                    vector<double> l1c_q(dim_pol_wave,0.0);
                    vector<double> l1c_q_noise(dim_pol_wave,0.0);
                    vector<double> l1c_u(dim_pol_wave,0.0);
                    vector<double> l1c_u_noise(dim_pol_wave,0.0);
                    vector<double> l1c_qu_noise_covar(dim_pol_wave,0.0);
                    vector<double> l1c_dolp(dim_pol_wave,0.0);
                    vector<double> l1c_dolp_noise(dim_pol_wave,0.0);
                    vector<double> l1c_aolp(dim_pol_wave,0.0);
                    vector<double> l1c_aolp_noise(dim_pol_wave,0.0);
                    double *l1b_i_polsample_cur = l1b_i_polsample.data();
                    double *l1b_i_polsample_noise_cur = l1b_i_polsample_noise.data();
                    double *l1b_q_cur = l1b_q.data();
                    double *l1b_q_noise_cur = l1b_q_noise.data();
                    double *l1b_u_cur = l1b_u.data();
                    double *l1b_u_noise_cur = l1b_u_noise.data();
                    double *l1b_qu_noise_covar_cur = l1b_qu_noise_covar.data();
                    for (size_t ialong=0 ; ialong<dim_along ; ialong++) {
                        for (size_t iacross=0 ; iacross<dim_across ; iacross++) {
                            for (size_t iwave=0 ; iwave<dim_pol_wave ; iwave++) {
                                l1c_i_polsample[iwave] += factormat_cur[ialong*dim_across+iacross] * l1b_i_polsample_cur[iwave];
                                l1c_i_polsample_noise[iwave] += pow(factormat_cur[ialong*dim_across+iacross] * l1b_i_polsample_noise_cur[iwave],2.0); // Turned into variances.
                                l1c_q[iwave] += factormat_cur[ialong*dim_across+iacross] * l1b_q_cur[iwave];
                                l1c_q_noise[iwave] += pow(factormat_cur[ialong*dim_across+iacross] * l1b_q_noise_cur[iwave],2.0); // Turned into variances.
                                l1c_u[iwave] += factormat_cur[ialong*dim_across+iacross] * l1b_u_cur[iwave];
                                l1c_u_noise[iwave] += pow(factormat_cur[ialong*dim_across+iacross] * l1b_u_noise_cur[iwave],2.0); // Turned into variances.
                                l1c_qu_noise_covar[iwave] += pow(factormat_cur[ialong*dim_across+iacross],2.0) * l1b_qu_noise_covar_cur[iwave]; // This is already a variance. Still averages with averaging out noise assuming there is no correlation between one Q and the other U.
                            }
                            l1b_i_polsample_cur += dim_pol_wave;
                            l1b_i_polsample_noise_cur += dim_pol_wave;
                            l1b_q_cur += dim_pol_wave;
                            l1b_q_noise_cur += dim_pol_wave;
                            l1b_u_cur += dim_pol_wave;
                            l1b_u_noise_cur += dim_pol_wave;
                            l1b_qu_noise_covar_cur += dim_pol_wave;
                        }
                    }
                    for (size_t iwave=0 ; iwave<dim_pol_wave ; iwave++) {
                        l1c_dolp[iwave] = sqrt(pow(l1c_q[iwave],2.0) + pow(l1c_u[iwave],2.0));
                        l1c_aolp[iwave] = 0.5*atan2(l1c_u[iwave],l1c_q[iwave]) / DEGREES;
                        l1c_dolp_noise[iwave] = sqrt(
                            pow(l1c_q[iwave] / l1c_dolp[iwave],2.0) * l1c_q_noise[iwave] + // Last term is still a variance.
                            pow(l1c_u[iwave] / l1c_dolp[iwave],2.0) * l1c_u_noise[iwave] + // Last term is still a variance.
                            l1c_q[iwave]*l1c_u[iwave] / pow(l1c_dolp[iwave],2.0) * l1c_qu_noise_covar[iwave]
                        );
                        l1c_aolp_noise[iwave] = 0.5 * sqrt(
                            pow(l1c_u[iwave] / pow(l1c_dolp[iwave],2.0),2.0) * l1c_q_noise[iwave] + // Last term is still a variance.
                            pow(l1c_q[iwave] / pow(l1c_dolp[iwave],2.0),2.0) * l1c_u_noise[iwave] - // Last term is still a variance.
                            l1c_q[iwave]*l1c_u[iwave] / pow(l1c_dolp[iwave],4.0) * l1c_qu_noise_covar[iwave]
                        ) / DEGREES;
                        // Now turn the variances on I, q and u into noise.
                        l1c_i_polsample_noise[iwave] = sqrt(l1c_i_polsample_noise[iwave]);
                        l1c_q_noise[iwave] = sqrt(l1c_q_noise[iwave]);
                        l1c_u_noise[iwave] = sqrt(l1c_u_noise[iwave]);
                    }
                    // Write the results.
                    netcdf_check(&nc_l1c,var_l1c_i_polsample.putVar(l1c_start,l1c_cnt,l1c_i_polsample.data()));
                    netcdf_check(&nc_l1c,var_l1c_i_polsample_noise.putVar(l1c_start,l1c_cnt,l1c_i_polsample_noise.data()));
                    netcdf_check(&nc_l1c,var_l1c_q.putVar(l1c_start,l1c_cnt,l1c_q.data()));
                    netcdf_check(&nc_l1c,var_l1c_q_noise.putVar(l1c_start,l1c_cnt,l1c_q_noise.data()));
                    netcdf_check(&nc_l1c,var_l1c_u.putVar(l1c_start,l1c_cnt,l1c_u.data()));
                    netcdf_check(&nc_l1c,var_l1c_u_noise.putVar(l1c_start,l1c_cnt,l1c_u_noise.data()));
                    netcdf_check(&nc_l1c,var_l1c_dolp.putVar(l1c_start,l1c_cnt,l1c_dolp.data()));
                    netcdf_check(&nc_l1c,var_l1c_dolp_noise.putVar(l1c_start,l1c_cnt,l1c_dolp_noise.data()));
                    netcdf_check(&nc_l1c,var_l1c_aolp.putVar(l1c_start,l1c_cnt,l1c_aolp.data()));
                    netcdf_check(&nc_l1c,var_l1c_aolp_noise.putVar(l1c_start,l1c_cnt,l1c_aolp_noise.data()));
                }
                factormat_cur += dim_along*dim_across;
            }
        }

        writelog(log_trace,"Interpolating solar and viewing geometry.");
        // Geometry.
        {
            vector<double> l1b_sza(dim_along*dim_across);
            vector<double> l1b_saa(dim_along*dim_across);
            vector<double> l1b_vza(dim_along*dim_across);
            vector<double> l1b_vaa(dim_along*dim_across);
            vector<Vector> l1b_sun(dim_along*dim_across);
            vector<Vector> l1b_view(dim_along*dim_across);

            // Slice out one viewport.
            vector<size_t> l1b_start = {0,ifov_start[ivp]};
            vector<size_t> l1b_cnt = {dim_along,dim_across};
            netcdf_check(&nc_l1b,var_l1b_sza.getVar(l1b_start,l1b_cnt,l1b_sza.data()));
            netcdf_check(&nc_l1b,var_l1b_saa.getVar(l1b_start,l1b_cnt,l1b_saa.data()));
            netcdf_check(&nc_l1b,var_l1b_vza.getVar(l1b_start,l1b_cnt,l1b_vza.data()));
            netcdf_check(&nc_l1b,var_l1b_vaa.getVar(l1b_start,l1b_cnt,l1b_vaa.data()));
            // Turn angles into something interpolatable.
            for (size_t imat=0 ; imat<dim_along*dim_across ; imat++) {
                // Go to the correct units.
                l1b_sza[imat] *= DEGREES;
                l1b_saa[imat] *= DEGREES;
                l1b_vza[imat] *= DEGREES;
                l1b_vaa[imat] *= DEGREES;
                // Now, turn them into vectors.
                // We have to interpret zenith as altitude and azimuth
                // as longitude. But for the zenith, we need to flip
                // sine and cosine, because it is a zero (or 180 degrees)
                // angle that is in the singularity where the azimuth
                // does not matter.
                l1b_sun[imat] = Vector (
                    sin(l1b_sza[imat]) * cos(l1b_saa[imat]) ,
                    sin(l1b_sza[imat]) * sin(l1b_saa[imat]) ,
                    cos(l1b_sza[imat])
                );
                l1b_view[imat] = Vector (
                    sin(l1b_vza[imat]) * cos(l1b_vaa[imat]) ,
                    sin(l1b_vza[imat]) * sin(l1b_vaa[imat]) ,
                    cos(l1b_vza[imat])
                );
            }
            factormat_cur = factormat.data();
            for (size_t iobs=0 ; iobs<dim_l1c_obs ; iobs++) {

                // Where to write.
                vector<size_t> l1c_start = {iobs,0,ivp};
                vector<size_t> l1c_cnt = {1,1,1};

                if (fail[iobs]) {

                    // Write fill values.
                    double fillvalue = NC_FILL_DOUBLE;
                    netcdf_check(&nc_l1c,var_l1c_sza.putVar(l1c_start,l1c_cnt,&fillvalue));
                    netcdf_check(&nc_l1c,var_l1c_saa.putVar(l1c_start,l1c_cnt,&fillvalue));
                    netcdf_check(&nc_l1c,var_l1c_vza.putVar(l1c_start,l1c_cnt,&fillvalue));
                    netcdf_check(&nc_l1c,var_l1c_vaa.putVar(l1c_start,l1c_cnt,&fillvalue));

                } else {

                    Vector l1c_sun(0.0,0.0,0.0);
                    Vector l1c_view(0.0,0.0,0.0);
                    Vector *l1b_sun_cur = l1b_sun.data();
                    Vector *l1b_view_cur = l1b_view.data();
                    for (size_t ialong=0 ; ialong<dim_along ; ialong++) {
                        for (size_t iacross=0 ; iacross<dim_across ; iacross++) {
                            l1c_sun += factormat_cur[ialong*dim_across+iacross] * *l1b_sun_cur;
                            l1c_view += factormat_cur[ialong*dim_across+iacross] * *l1b_view_cur;
                            l1b_sun_cur++;
                            l1b_view_cur++;
                        }
                    }
                    // Convert vectors to angles.
                    // Note that the averged vector is not normalized.
                    l1c_sun.normalize();
                    l1c_view.normalize();
                    // Now they are.
                    double l1c_sza = acos(l1c_sun[2]) / DEGREES;
                    double l1c_saa = atan2(l1c_sun[1],l1c_sun[0]) / DEGREES;
                    double l1c_vza = acos(l1c_view[2]) / DEGREES;
                    double l1c_vaa = atan2(l1c_view[1],l1c_view[0]) / DEGREES;

                    // Write the results.
                    netcdf_check(&nc_l1c,var_l1c_sza.putVar(l1c_start,l1c_cnt,&l1c_sza));
                    netcdf_check(&nc_l1c,var_l1c_saa.putVar(l1c_start,l1c_cnt,&l1c_saa));
                    netcdf_check(&nc_l1c,var_l1c_vza.putVar(l1c_start,l1c_cnt,&l1c_vza));
                    netcdf_check(&nc_l1c,var_l1c_vaa.putVar(l1c_start,l1c_cnt,&l1c_vaa));
                }
                factormat_cur += dim_along*dim_across;
            }
        }

    } // Loop over viewports.

    return 0;

} // }}}

