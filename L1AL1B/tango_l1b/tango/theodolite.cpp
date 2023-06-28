// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "vector.h"
#include "matrix.h"
#include "lininv.h"
#include "netcdf_object.h"
#include "theodolite.h"

// Constructor.
Theodolite::Theodolite( // {{{
    Logger *creator
) : Logger(creator)
{
} // }}}
Theodolite::~Theodolite(){}

// Do what you can do with the reference file.
int Theodolite::calculate_reference( // {{{
    size_t a_dim_vp,
    string reference_filename,
    double *viewport_alt_angles,
    double alt_angle_tol,
    size_t altinv_maxiter,
    double altinv_tol,
    double altinv_steptolerance,
    double altinv_initial_step_reducer,
    double altinv_reducerfactor_fail,
    double altinv_reducerfactor_success,
    double altinv_reducer_limit
)
{

    // Save number of viewports as quickly as possible.
    dim_vp = a_dim_vp;

    // File contents.
    size_t nmeas;
    vector<double> act;
    vector<double> alt;
    vector<double> ang;
    vector<double> stst; // Not used in reference file.

    handle(read(
        reference_filename,
        nmeas,
        act,
        alt,
        tar,
        ang,
        stst
    ));
    // We will ignore (and probably not have) stst.

    // Recognize viewports by ALT angle. Meanwhile, the viewports should
    // be ordered properly. It is possible (and true) that the viewports
    // from 0 to dim_vp-1 are not in ascending order of ALT stage angle.

    // Therefore, we will organize the measurements to the closest expected
    // ALT angle of the viewport. If the difference is higher than a
    // tolerance, the measurement is discarded and a warning is raised.

    vector<vector<size_t>> indexarray_vp(dim_vp); // Index array per viewport.
    for (size_t imeas=0 ; imeas<nmeas ; imeas++) {
        // Only one viewport can meet the requirement. That has already
        // been verified.
        bool found = false;
        for (size_t ivp=0 ; ivp<dim_vp ; ivp++) {
            double diff = abs(alt[imeas] - viewport_alt_angles[ivp]);
            if (diff < alt_angle_tol) {
                indexarray_vp[ivp].emplace_back(imeas);
                found = true;
                break; // Not needed, but why not.
            }
        }
        if (!found) {
            writelog(log_warning,"Warning: Theodolite file contains measurement with stage ALT angle %.7f, which is not recognized in any of the viewports.",alt[imeas]);
        }
    }

    // Fit linear relationship between theodolite angles for each viewport
    // for differences in ACT stage angle.
    diff_theo_angs_act.resize(dim_vp*ntar*nhv);
    // Save offsets for ALT model.
    vector<double> theo_angs_act0(dim_vp*ntar*nhv);
    // Meanwhile, get the actual ALT angles for each viewport.
    actual_alt_angles.resize(dim_vp,0.0);
    for (size_t ivp=0 ; ivp<dim_vp ; ivp++) {
        size_t nmeas_vp = indexarray_vp[ivp].size();
        size_t nstate = 2; // Constant and linear term.

        // Fill Jacobian. It is constant, because it only depends on
        // stage ACT angles.
        vector<double> jac(nstate*nmeas_vp);
        for (size_t imeas_vp=0 ; imeas_vp<nmeas_vp ; imeas_vp++) {
            size_t &imeas = indexarray_vp[ivp][imeas_vp];
            jac[imeas_vp] = 1.0; // State parameter 0: Constant term.
            jac[nmeas_vp+imeas_vp] = act[imeas]; // Linear term: Equal to abscissa.
            // Meanwhile, get the mean value of the ALT stage angle.
            actual_alt_angles[ivp] += alt[imeas] / static_cast<double>(nmeas_vp);

        }
        for (size_t itar=0 ; itar<ntar ; itar++) {
            for (size_t ihv=0 ; ihv<nhv ; ihv++) {
                // Fill measurement vector.
                vector<double> meas(nmeas_vp);
                for (size_t imeas_vp=0 ; imeas_vp<nmeas_vp ; imeas_vp++) {
                    size_t &imeas = indexarray_vp[ivp][imeas_vp];
                    meas[imeas_vp] = ang[imeas*ntar*nhv + itar*nhv + ihv];
                }
                // Linear inversion.
                vector<double> res(nstate);
                check_error(linear_invert(nstate,nmeas_vp,jac.data(),OPT_NONE,NULL,NULL,meas.data(),res.data()) != 0,"Error: Failed fitting linear relationship of theodolite %s angle for IAC %s face for viewport %zu.",hvnames[ihv].c_str(),targetnames[tar[itar]].c_str(),ivp);
                // Write slope result.
                diff_theo_angs_act[ivp*ntar*nhv + itar*nhv + ihv] = res[1];
                // Also write the offset result. Those are, namely the
                // theodolite angles with ACT=0, which we want to have
                // for the model for ALT-dependence.
                theo_angs_act0[ivp*ntar*nhv + itar*nhv + ihv] = res[0];
            }
        }
    }

    // Solve ALT-model.
    // The measurements are the IAC orientations (expressed as matrix
    // or quaternion) for each viewport with ACT=0.
    // The state is a zero-orientation for ALT=0 and a rotation axis
    // with which it is rotated about with prescribed magnitude.
    // There are five DFS in the state, three for zero-orientation
    // and two for the rotation axis.
    // To do an inversion, we desire to have restrict to legal situations
    // while allowing the five state parameters any value.
    // For the base quaternion, this is doable. We set a real part to one,
    // and let the three imaginary parts do what they want. Then, we
    // normalize the quaternion and use it. This forces the real part to
    // stay positive, which is not a restriction.
    // For the rotation axis, we have to do something similar, forcing
    // the X, Y or Z component to be either positive or negative. Angles
    // are defined with sign, so involves a restriction. To avoid
    // removing the true solution from the state space, we test which
    // of the six axes: +/- X, Y or Z is closest. That is the restriction
    // that we take. So, if -Y is the best of the six, we take rotation
    // axis X, -1, Z and then normalize it with X and Z as state variables.
    // For this investigation, we do not have a zero-position, so we
    // force the zero-position to the nadir viewport and express the
    // angles as the differences with that viewport. This overweighs one
    // point, but for a pre-fit, that should not matter too much.

    // First, convert the five measurements into quaternions.
    vector<Quaternion> alt_model_iacs(dim_vp);
    for (size_t ivp=0 ; ivp<dim_vp ; ivp++) {
        handle(theo_to_quat(&theo_angs_act0[ivp*ntar*nhv],&alt_model_iacs[ivp]));
    }

    // Find out which viewport is the most nadir.
    size_t ivp_nadir = 0;
    for (size_t ivp=1 ; ivp<dim_vp ; ivp++) {
        if (abs(actual_alt_angles[ivp]) < abs(actual_alt_angles[ivp_nadir])) ivp_nadir = ivp; // Overwrites.
    }
    // Calculate ALT angle differences.
    vector<double> alt_angle_differences(dim_vp);
    for (size_t ivp=0 ; ivp<dim_vp ; ivp++) {
        alt_angle_differences[ivp] = actual_alt_angles[ivp] - actual_alt_angles[ivp_nadir];
    }

    // Get the rotation quaternions by ALT stage change by comparing
    // IAC orientations.
    vector<Quaternion> rotations_meas(dim_vp);
    // The one for nadir will certainly be 1 (no rotation).
    for (size_t ivp=0 ; ivp<dim_vp ; ivp++) {
        rotations_meas[ivp] = alt_model_iacs[ivp] / alt_model_iacs[ivp_nadir];
    }
    // The idea is to check out which imaginary element of the measured
    // rotation is depending most signficanly on the ALT-angle in whether
    // it is ascending or descending.
    // So, we will fit the slope of the quaternion element with sin(ALT/2).
    // Let us hope that is sufficient.
    size_t nstate = 2;
    vector<double> jac(nstate*dim_vp);
    for (size_t ivp=0 ; ivp<dim_vp ; ivp++) {
        jac[ivp] = 1.0;
        jac[dim_vp+ivp] = sin(actual_alt_angles[ivp]/2.0);
    }
    double best = 0.0;
    size_t ivec_fixed = NC_FILL_UINT64;
    bool neg;
    for (size_t ivec=0 ; ivec<DIM_VEC ; ivec++) {
        vector<double> meas(dim_vp);
        for (size_t ivp=0 ; ivp<dim_vp ; ivp++) {
            meas[ivp] = alt_model_iacs[ivp][ivec];
        }
        vector<double> res(nstate);
        check_error(linear_invert(nstate,dim_vp,jac.data(),OPT_NONE,NULL,NULL,meas.data(),res.data()) != 0,"Error: Failed fitting linear relationship ofquaternion elements %zu with the sine of half the ALT angle",ivec);
        // Check if this is the best up until now.
        if (abs(res[1]) > best) {
            // Write results into 'best'.
            best = abs(res[1]);
            ivec_fixed = ivec;
            neg = res[1] < 0.0;
        }
    }
    check_error(ivec_fixed == NC_FILL_UINT64,"Error: No correlation found between rotation quaternions of ALT-model with ALT angle.");

    Altmodel_Inversion altinv(this,ivec_fixed,neg,actual_alt_angles);
    altinv.setMaxiter(altinv_maxiter);
    vector<double> apriori(altinv.getNstate(),0.0);
    altinv.setApriori(apriori.data());
    vector<double> tolerance(altinv.getNstate(),altinv_tol);
    altinv.setTol(tolerance.data());
    altinv.setSteptolerance(altinv_steptolerance);
    altinv.setInitialStepReducer(altinv_initial_step_reducer);
    altinv.setReducerfactorFail(altinv_reducerfactor_fail);
    altinv.setReducerfactorSuccess(altinv_reducerfactor_success);
    altinv.setReducerLimit(altinv_reducer_limit);
    vector<double> meas(dim_vp*DIM_QUAT);
    for (size_t ivp=0 ; ivp<dim_vp ; ivp++) {
        alt_model_iacs[ivp].get(&meas[ivp*DIM_QUAT]);
    }
    altinv.setMeas(meas.data());
    vector<double> retrieval_result(altinv.getNstate());
    int conv = altinv.inversion_calculate(retrieval_result.data());
    check_error(conv == 2,"Error: Numeric error during non-linear ALT-model retrieval.");
    check_error(conv == 1,"Error: Non-linear ALT-model retrieval did not converge.");
    // Shape final ALT-model results.
    altmodel_base = make_unique<Quaternion>();
    altmodel_rotax = make_unique<Vector>();
    // Reperat forward model with retrieval results to obtain base quaternion and rotation axis.
    altinv.translate(retrieval_result.data(),altmodel_base.get(),altmodel_rotax.get());

    return 0;

} // }}}

// Add information from spot check file.
int Theodolite::calculate_spotcheck( // {{{
    string spotcheck_filename
)
{

    // File contents.
    size_t nmeas;
    vector<double> act;
    vector<double> alt;
    vector<int> tar_spot;
    vector<double> ang;
    vector<double> stst;

    handle(read(
        spotcheck_filename,
        nmeas,
        act,
        alt,
        tar_spot,
        ang,
        stst
    ));

    for (size_t itar=0 ; itar<ntar ; itar++) {
        check_error(tar_spot[itar] != tar[itar],"Error: Inconsistent target %zu for reference (%d) and spot check (%d)",itar,tar[itar],tar_spot[itar]);
    }

    // There are two things in the spotcheck file: The offset from one or a few normal
    // theodolite measurements and a star stimulus.

    check_error(nmeas == 0,"Error: No measurements in spot check file.");

    // Get theodolite angles from model with zero offset. Compare them with the actual
    // theodolite angles and the difference becomes the offset. Averaged over measurements.
    offset = {0.0,0.0,0.0,0.0}; // H and V for two targets.
    vector<double> calculated_offset = {0.0,0.0,0.0,0.0};
    // We have to keep offset to zero, because it is used during stage_to_theo, because the
    // same routine is used when calculating swath vectors (thus with offset).
    for (size_t imeas=0 ; imeas<nmeas ; imeas++) {
        vector<double> expected_angs(ntar*nhv);
        handle(stage_to_theo(
            act[imeas],
            alt[imeas],
            expected_angs.data()
        ));
        for (size_t itar=0 ; itar<ntar ; itar++) {
            for (size_t ihv=0 ; ihv<nhv ; ihv++) {
                double diff = ang[imeas*ntar*nhv+itar*nhv+ihv] - expected_angs[itar*nhv+ihv];
                // Wrap the horizontal (longitude-like) difference.
                if (ihv == 0) {
                    while (diff > 180.0*DEGREES) diff -= 360.0*DEGREES;
                    while (diff < -180.0*DEGREES) diff += 360.0*DEGREES;
                }
                calculated_offset[itar*nhv+ihv] += diff / static_cast<double>(nmeas);
            }
        }
    }
    // Now put the offset into the member variable.
    offset = calculated_offset;

    // Star stimulus.
    // It must be there.
    check_error(stst.size() == 0,"Error: Spot check file '%s' does not contain a star stimulus.",spotcheck_filename.c_str());
    // TODO: Figure out if this vector has to be taken positive or negative. Depends on the
    // definition of a cube face sign.
    star_stimulus = make_unique<Vector> (
        // X = cos(V).
        cos(stst[1]), // The 1 stands for V.
        // Y = cos(H)*sin(v).
        cos(stst[0])*sin(stst[1]),
        // Z = sin(H)*sin(v).
        sin(stst[0])*sin(stst[1])
    );

    return 0;

} // }}}

// Get swath vector from stage angle using the theodolite measurements.
int Theodolite::calculate_swath_vector( // {{{
    double act,
    double alt,
    Vector &swathvector
)
{

    // 1. Convert the stage into expected theodolite angles.
    // 2. Convert these theodolite angles into a quaternion.
    // 3. Compare this quaternion (=matrix) with the star stimulus direction vector.

    vector<double> theo_angs(ntar*nhv);
    handle(stage_to_theo(act,alt,theo_angs.data()));
    Quaternion quat;
    handle(theo_to_quat(theo_angs.data(),&quat));
    // TODO: Verify transpose or not.
    for (size_t ax=0 ; ax<DIM_VEC ; ax++) {
        Vector vec;
        handle(quat_to_unitvec(quat,ax,vec));
        swathvector[ax] = star_stimulus->dot(vec);
    }

    return 0;

} // }}}

// Data member getters.
size_t Theodolite::getDimVp() {return dim_vp;}

// Private routines.

// Read theodolite file.
int Theodolite::read( // {{{
    string filename,
    size_t &nmeas,
    vector<double> &act,
    vector<double> &alt,
    vector<int> &tar,
    vector<double> &ang,
    vector<double> &stst
)
{

    // Acquire the information from the theodolite measurements.
    NetCDF_object nc(this);
    handle(nc.open(filename,NcFile::read));

    // Read the contents of the NetCDF file.
    netcdf_check(&nc,nmeas = nc.ncid->getDim("nmeas").getSize());
    act.resize(nmeas);
    alt.resize(nmeas);
    size_t dimcheck;
    netcdf_check(&nc,dimcheck = nc.ncid->getDim("ntar").getSize());
    check_error(dimcheck != ntar,"Error in theodolite file '%s': Theodolite file has %zu targets. It should be %zu.",filename.c_str(),dimcheck,ntar);
    netcdf_check(&nc,dimcheck = nc.ncid->getDim("nhv").getSize());
    check_error(dimcheck != nhv,"Error in theodolite file '%s': Theodolite file has %zu-sized H/V dimension. It should be %zu.",filename.c_str(),dimcheck,nhv);
    vector<int> rail(nmeas*ntar);
    vector<int> tar_full(nmeas*ntar); // Targets of all measurements.
    // All measurements must have consistent targets. At the end, one set
    // of targets is returned.
    tar.resize(ntar);
    ang.resize(nmeas*ntar*nhv);
    netcdf_check(&nc,nc.ncid->getVar("stage_act_angle").getVar(act.data()));
    netcdf_check(&nc,nc.ncid->getVar("stage_alt_angle").getVar(alt.data()));
    netcdf_check(&nc,nc.ncid->getVar("rail_number").getVar(rail.data()));
    netcdf_check(&nc,nc.ncid->getVar("target").getVar(tar_full.data()));
    netcdf_check(&nc,nc.ncid->getVar("theodolite_angles").getVar(ang.data()));
    // Convert degrees to radians.
    for (size_t imeas=0 ; imeas<nmeas ; imeas++) {
        act[imeas] *= DEGREES;
        alt[imeas] *= DEGREES;
    }
    for (size_t iel=0 ; iel<nmeas*ntar*nhv ; iel++) ang[iel] *= DEGREES;
    // Put all measurements on rail 1, by shifting H for all measurements
    // on rail 2.
    for (size_t iel=0 ; iel<nmeas*ntar ; iel++) {
        if (rail[iel] == 2) {
            ang[iel*nhv] += 90.0*DEGREES;
        } else {
            check_error(rail[iel] != 1,"Error: Unrecognized rail number: %d",rail[iel]);
        }
    }
    // Verify that every measurement has the same targets.
    for (size_t imeas=1 ; imeas<nmeas ; imeas++) {
        check_error(tar_full[imeas*ntar] != tar_full[0],"Error: Inconsistent targets: %s and %s for first target in theodolite file.",targetnames[tar_full[imeas*ntar]].c_str(),targetnames[tar_full[0]].c_str());
        check_error(tar_full[imeas*ntar+1] != tar_full[1],"Error: Inconsistent targets: %s and %s for second target in theodolite file.",targetnames[tar_full[imeas*ntar+1]].c_str(),targetnames[tar_full[1]].c_str());
    }
    tar = {tar_full[0],tar_full[1]}; // Keep what we need. Throw away the rest.

    // Optional: The star stimulus.
    NcVar var_stst;
    netcdf_check(&nc,var_stst = nc.ncid->getVar("star_stimulus"));
    if (!var_stst.isNull()) {
        stst.resize(nhv);
        netcdf_check(&nc,var_stst.getVar(stst.data()));
        // Also read the rail number. If there is a star stimulus, there must also be rail
        // associated to it.
        int ststrail;
        netcdf_check(&nc,nc.ncid->getVar("star_stimulus_rail_number").getVar(&ststrail));
        // Convert degrees to radians.
        for (size_t ihv=0 ; ihv<nhv ; ihv++) {
            stst[ihv] *= DEGREES;
        }
        // Correct for rail number.
        if (ststrail == 2) {
            stst[0] += 90.0*DEGREES;
        } else {
            check_error(ststrail != 1,"Error: Unrecognized rail number for star stimulus: %d",ststrail);
        }
    }

    return 0;

} // }}}

// Turn theodolite angles into orientation Quaternion.
int Theodolite::theo_to_quat( // {{{
    double *theo_angs,
    Quaternion *quat
)
{

    // This is the matrix, but then as vector of Vectors.
    vector<Vector> basevectors(DIM_VEC);

    vector<size_t> axes(ntar); // Save target axes for post-processing.
    for (size_t itar=0 ; itar<ntar ; itar++) {
        // Convert six possible targets to an axis and a negative flag.
        int ax = tar[itar] / 2; // 0=X, 1=Y, 2=Z.
        double fact = tar[itar] - 2*ax == 0?-1.0:1.0; // Toggle negative and positive.
        basevectors[ax] = Vector(
            // X = cos(V).
            fact*cos(theo_angs[itar*nhv+1]), // The 1 stands for V.
            // Y = cos(H)*sin(v).
            fact*cos(theo_angs[itar*nhv])*sin(theo_angs[itar*nhv+1]),
            // Z = sin(H)*sin(v).
            fact*sin(theo_angs[itar*nhv])*sin(theo_angs[itar*nhv+1])
        );

        // Write down axis.
        axes[itar] = ax;

    }

    // Verify that the two theodolite are not looking at the same or
    // excactly at the opposite side of the reference cube.
    check_error(axes[0] == axes[1],"Error: Cannot interpret theodolite measurements if the two targets (%s and %s) are the same axis.",targetnames[tar[0]].c_str(),targetnames[tar[1]].c_str());

    // I assume that third axis = first target <cross> second target.
    // With axes turned positive. This is true if
    // First target: X (0), second target: Y (1).
    // First target: Y (1), second target: Z (2).
    // First target: Z (2), second target: X (0).
    int axdiff = static_cast<int>(axes[1])-static_cast<int>(axes[0]);
    double fact = (axdiff == 1 || axdiff == -2)?1.0:-1.0;
    size_t ax_third = 3-axes[0]-axes[1]; // Remaining axes.
    // The two calculated base vectors may not be exactly perpendicular.
    // This is fixed by adapting the second target.
    basevectors[ax_third] = fact*basevectors[axes[0]].cross(basevectors[axes[1]]);
    basevectors[ax_third].normalize(); // This does something only if there is a imperfection.
    // Progress to second target. With the same factor, we can assume that.
    // second target = third axis <cross> first target.
    basevectors[axes[1]] = fact*basevectors[ax_third].cross(basevectors[axes[0]]);

    // Last step is to convert the base vectors into a Quaternion. It may
    // be a bit ambiguous in what direction we define this quaternion.
    // Are we rotating the axes or the vector inside the axes, resulting
    // into conjugate results. Check carefully if the world is not turned
    // upside down.
    double realpart = 0.5*sqrt(basevectors[0][0]+basevectors[1][1]+basevectors[2][2]+1.0);
    *quat = Quaternion(
        (basevectors[1][2]-basevectors[2][1]) / (4.0*realpart),
        (basevectors[2][0]-basevectors[0][2]) / (4.0*realpart),
        (basevectors[0][1]-basevectors[1][0]) / (4.0*realpart),
        realpart
    );

    return 0;

} // }}}

// Do the exact opposite: Turn orientation quaternion into theodolite angles.
int Theodolite::quat_to_theo( // {{{
    Quaternion *quat, // Input orientation quaternion.
    double *theo_angs // Output theodolite angles.
)
{

    for (size_t itar=0 ; itar<ntar ; itar++) {
        int ax = tar[itar] / 2; // 0=X, 1=Y, 2=Z.
        Vector vec;
        handle(quat_to_unitvec(*quat,ax,vec));
        if (tar[itar] - 2*ax == 0) vec *= -1.0;

        // Convert them into theodolite angles.
        // X = cos(V).
        // Y = cos(H)*sin(v).
        // Z = sin(H)*sin(v).
        theo_angs[itar*nhv] = atan2(vec[2],vec[1]); // H.
        theo_angs[itar*nhv+1] = acos(vec[0]); // V.

    }

    return 0;

} // }}}

int Theodolite::quat_to_unitvec( // {{{
    Quaternion quat,
    int ax, // 0=X, 1=Y, 2=Z.
    Vector &vec
)
{

    // Get other axes in right-handed order (X - Y - Z - X).
    size_t next_ax;
    size_t third_ax;
    next_ax = ax + 1;
    if (next_ax == DIM_VEC) next_ax = 0;
    third_ax = next_ax + 1;
    if (third_ax == DIM_VEC) third_ax = 0;
    // Diagonal term: 1 - squares of imaginary quaternion elements other than ax.
    vec[ax] = (1.0 - 2.0*pow(quat[next_ax],2.0) - 2.0*pow(quat[third_ax],2.0));
    // Next-ax, terms with plus sign.
    vec[next_ax] = 2.0*(quat[ax]*quat[next_ax] + quat[third_ax]*quat[DIM_VEC]);
    // Third-ax, terms with minus sign.
    vec[third_ax] = 2.0*(quat[ax]*quat[third_ax] - quat[next_ax]*quat[DIM_VEC]);

    return 0;

} // }}}

// Acquire thedolite angles that correspond to a certain stage position.
int Theodolite::stage_to_theo( // {{{
    double act, // Across-track stage angle.
    double alt, // Along-track stage angle.
    double *theo_angs // Output theodolite angles.
)
{

    // 1: Evaluate the ALT model to convert the ALT angle into an orientation quaternion
    // for ACT = 0.
    // 2: Convert this quaternion into H and V.
    // 3: Use ACT slopes to include nonzero ACT influence.
    // 4: Add offset in H and V from spot check file.
    Quaternion altmodel_rot(sin(alt/2.0)*(*altmodel_rotax),cos(alt/2.0));
    Quaternion quat = altmodel_rot * (*altmodel_base);
    handle(quat_to_theo(&quat,theo_angs));
    // Search best viewport (for as long as ACT-derivatives are indeed
    // calculated per viewport).
    size_t ivp_best = 0;
    double diff_alt = abs(alt-actual_alt_angles[0]);
    for (size_t ivp=1 ; ivp<dim_vp ; ivp++) {
        double newdiff = abs(alt-actual_alt_angles[ivp]);
        if (newdiff < diff_alt) {
            diff_alt = newdiff;
            ivp_best = ivp;
        }
    }
    // Apply ACT derivative from best viewport and offset.
    for (size_t iel=0 ; iel<ntar*nhv ; iel++) {
        theo_angs[iel] += act*diff_theo_angs_act[ivp_best*ntar*nhv+iel] + offset[iel];
    }

    return 0;

} // }}}

// Inversion routine for the ALT-model.
// Constructor.
Altmodel_Inversion::Altmodel_Inversion( // {{{
    Theodolite *creator, // Creating theodolite instance.
    size_t a_ivec_fixed, // Dimension of rotation axis to set to one (before normalization).
    bool neg, // True if it is minus one instead of plus one.
    vector<double> a_alt_angles // ALT stage angles at measurement points.
) : Inversion(
    creator, // The theodolite is also the logger creator.
    5, // This model always has five state parameters.
    creator->getDimVp()*DIM_QUAT // All orientation quaternions will be compared. That is one quaternion (DIM_QUAT values) per viewport.
)
{

    // Set mapping indices for fourth and fifth state parameters.
    ivec_fixed = a_ivec_fixed;
    ivec_rotax_first = ivec_fixed == 0?1:0;
    ivec_rotax_second = ivec_fixed == 2?1:2;
    val_fixed = neg?-1.0:1.0;
    // Copy ALT-angles (abscissa). This is a true copy of dim_vp values.
    alt_angles = a_alt_angles;

    // Save number of viewports. Yes, it is nmeas/DIM_QUAT, but that is clumsy.
    dim_vp = creator->getDimVp();

} // }}}

// Empty destructor.
Altmodel_Inversion::~Altmodel_Inversion(){}

// Forward function for the ALT-model.
int Altmodel_Inversion::fwd( // {{{
    double *state, // State vector.
    double *signal, // Modelled signal (outout).
    double *jacobian // Modelled jacobian (outout).
)
{

    // Base quaternion. Define it with real part 1.
    Quaternion base(state[0],state[1],state[2],1.0);
    Quaternion diff_base_st0(1.0,0.0,0.0,0.0);
    Quaternion diff_base_st1(0.0,1.0,0.0,0.0);
    Quaternion diff_base_st2(0.0,0.0,1.0,0.0);
    // Magnitude of the quaternion.
    double mag = sqrt(
        pow(state[0],2.0) +
        pow(state[1],2.0) +
        pow(state[2],2.0) +
        1.0
    );
    double diff_mag_st0 = state[0] / mag;
    double diff_mag_st1 = state[1] / mag;
    double diff_mag_st2 = state[2] / mag;
    // Normalized base quaternion.
    Quaternion base_nrm = base / mag;
    Quaternion diff_base_nrm_st0 = (diff_base_st0 * mag - diff_mag_st0 * base) / pow(mag,2.0);
    Quaternion diff_base_nrm_st1 = (diff_base_st1 * mag - diff_mag_st1 * base) / pow(mag,2.0);
    Quaternion diff_base_nrm_st2 = (diff_base_st2 * mag - diff_mag_st2 * base) / pow(mag,2.0);

    // Rotation axis.
    Vector ax;
    ax[ivec_fixed] = val_fixed;
    ax[ivec_rotax_first] = state[3];
    ax[ivec_rotax_second] = state[4];
    Vector diff_ax_st3(0.0,0.0,0.0);
    diff_ax_st3[ivec_rotax_first] = 1.0;
    Vector diff_ax_st4(0.0,0.0,0.0);
    diff_ax_st4[ivec_rotax_second] = 1.0;
    // Magnitude of rotation axis.
    double axmag = sqrt(
        1.0 +
        pow(state[3],2.0) +
        pow(state[4],2.0)
    );
    double diff_axmag_st3 = state[3] / mag;
    double diff_axmag_st4 = state[4] / mag;
    // Normalized rotation axis.
    Vector ax_nrm = ax / axmag;
    Vector diff_ax_nrm_st3 = (diff_ax_st3 * axmag - diff_axmag_st3 * ax) / pow(axmag,2.0);
    Vector diff_ax_nrm_st4 = (diff_ax_st4 * axmag - diff_axmag_st4 * ax) / pow(axmag,2.0);

    // For each viewport.
    for (size_t ivp=0 ; ivp<dim_vp ; ivp++) {
        // Rotation quaternion for each viewport.
        double sinterm = sin(alt_angles[ivp]/2.0);
        double costerm = cos(alt_angles[ivp]/2.0);
        Quaternion rotquat(sinterm*ax_nrm,costerm);
        Quaternion diff_rotquat_st3(sinterm*diff_ax_nrm_st3,0.0);
        Quaternion diff_rotquat_st4(sinterm*diff_ax_nrm_st4,0.0);
        // And the final result.
        Quaternion res = rotquat * base_nrm;
        Quaternion diff_res_st0 = rotquat * diff_base_nrm_st0;
        Quaternion diff_res_st1 = rotquat * diff_base_nrm_st1;
        Quaternion diff_res_st2 = rotquat * diff_base_nrm_st2;
        Quaternion diff_res_st3 = diff_rotquat_st3 * base;
        Quaternion diff_res_st4 = diff_rotquat_st4 * base;
        // Write results into signal and jacobian.
        res.get(&signal[ivp*DIM_QUAT]);
        diff_res_st0.get(&jacobian[ivp*DIM_QUAT]);
        diff_res_st1.get(&jacobian[nmeas + ivp*DIM_QUAT]);
        diff_res_st2.get(&jacobian[2*nmeas + ivp*DIM_QUAT]);
        diff_res_st3.get(&jacobian[3*nmeas + ivp*DIM_QUAT]);
        diff_res_st4.get(&jacobian[4*nmeas + ivp*DIM_QUAT]);
    }

    return 0;

} // }}}

int Altmodel_Inversion::translate( // {{{
    double *state, // State vector.
    Quaternion *base, // Base quaternion.
    Vector *rotax // Rotation axis.
)
{

    // Define base quaternion and rotation axis as in the forward model.
    // Do not differentiate, so the function normalize() can be used.
    (*base)[0] = state[0];
    (*base)[1] = state[1];
    (*base)[2] = state[2];
    (*base)[3] = 1.0;
    base->normalize();
    (*rotax)[ivec_fixed] = val_fixed;
    (*rotax)[ivec_rotax_first] = state[3];
    (*rotax)[ivec_rotax_second] = state[4];
    rotax->normalize();

    return 0;

} // }}}

// DEBUG.
void Theodolite::logReferenceStuff(){
    writelog(log_debug,"dim_vp: %zu",dim_vp);
    writelog(log_debug,"tar.size(): %zu",tar.size());
    for (size_t itar=0 ; itar<tar.size() ; itar++) {
        writelog(log_debug,"tar[%zu]: %d",itar,tar[itar]);
    }
    writelog(log_debug,"actual_alt_angles.size(): %zu",actual_alt_angles.size());
    for (size_t ivp=0 ; ivp<actual_alt_angles.size() ; ivp++) {
        writelog(log_debug,"actual_alt_angles[%zu]: %.12f",ivp,actual_alt_angles[ivp]/DEGREES);
    }
    writelog(log_debug,"diff_theo_angs_act.size(): %zu",diff_theo_angs_act.size());
    for (size_t iel=0 ; iel<diff_theo_angs_act.size() ; iel++) {
        writelog(log_debug,"diff_theo_angs_act[%zu]: %.12f",iel,diff_theo_angs_act[iel]);
    }
    for (size_t iel=0 ; iel<DIM_QUAT ; iel++) {
        writelog(log_debug,"altmodel_base[%zu]: %.12f",iel,(*altmodel_base)[iel]);
    }
    for (size_t iel=0 ; iel<DIM_VEC ; iel++) {
        writelog(log_debug,"altmodel_rotax[%zu]: %.12f",iel,(*altmodel_rotax)[iel]);
    }
}

void Theodolite::logSpotcheckStuff(){
    for (size_t iel=0 ; iel<offset.size() ; iel++) {
        writelog(log_debug,"offset[%zu]: %.12f",iel,offset[iel]/DEGREES);
    }
    for (size_t iel=0 ; iel<DIM_VEC ; iel++) {
        writelog(log_debug,"star_stimulus[%zu]: %.12f",iel,(*star_stimulus)[iel]);
    }
}

