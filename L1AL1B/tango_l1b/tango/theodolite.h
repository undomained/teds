// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef THEODOLITE_H
#define THEODOLITE_H

#include "header.h"
#include "logger.h"
#include "inversion.h"

// Forward declaration.
class Vector;
class Quaternion;

// Some extra constants for logging.
const vector<string> targetnames = {
    "-X",
    "+X",
    "-Y",
    "+Y",
    "-Z",
    "+Z"
};
const vector<string> hvnames = {"H","V"};
// Static dimension sizes. It is nice to see if a '2' is
// a number of targets or a number of H/V theodolite angles.
const size_t ntar = 2;
const size_t nhv = 2;

class Theodolite : public Logger {

    public:
    // Constructor.
    Theodolite(
        Logger *creator
    );
    ~Theodolite();

    // Do what you can do with the reference file.
    int calculate_reference(
        size_t dim_vp,
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
    );

    // Add information from spot check file.
    int calculate_spotcheck(
        string spotcheck_filename
    );

    // Get swath vector from stage angle using the theodolite measurements.
    int calculate_swath_vector(
        double act,
        double alt,
        Vector &swathhvector
    );

    // Getters for private data members.
    size_t getDimVp();

    // Intermediate results are stored into member variables.
    private:
    size_t dim_vp; // Number of viewports.
    vector<int> tar; // First and second target. For now, this must be consistent all over the files.
    // Results from reference calculation.
    vector<double> actual_alt_angles; // ALT-angles per viewport, to recognize viewport.
    vector<double> diff_theo_angs_act; // Derivatives of theodolite angles with respect to ACT stage angle, for each target and viewport.
    unique_ptr<Quaternion> altmodel_base; // Base quaternion for the ALT angle dependence model.
    unique_ptr<Vector> altmodel_rotax; // Rotation axis for ALT angle dependence model.
    // Results from spotcheck calculation.
    vector<double> offset;
    unique_ptr<Vector> star_stimulus;

    // Read a theodolite file.
    int read(
        string filename,
        size_t &nmeas,
        vector<double> &act,
        vector<double> &alt,
        vector<int> &tar,
        vector<double> &ang,
        vector<double> &stst
    );

    // Turn theodolite anlges into orientation quaternion.
    int theo_to_quat(
        double *theo_angs, // Input theodolite angles.
        Quaternion *quat // Output orientation quaternion.
    );
    // Do the exact opposite: Turn orientation quaternion into theodolite angles.
    int quat_to_theo(
        Quaternion *quat, // Input orientation quaternion.
        double *theo_angs // Output theodolite angles.
    );
    // Get a unit vector from the matrix that belongs to a quaternion.
    int quat_to_unitvec(
        Quaternion quat,
        int ax,
        Vector &vec
    );
    // Acquire thedolite angles that correspond to a certain stage position.
    int stage_to_theo(
        double act, // Across-track stage angle.
        double alt, // Along-track stage angle.
        double *theo_angs // Output theodolite angles.
    );

    public:
    // DEBUG.
    void logReferenceStuff();
    void logSpotcheckStuff();

};

// The non-linear inversion for the ALT-stage model.
class Altmodel_Inversion : public Inversion {

    public:
    // Constructor.
    Altmodel_Inversion(
        Theodolite *creator, // Creating theodolite instance.
        size_t a_ivec_fixed, // Dimension of rotation axis to set to one (before normalization).
        bool neg, // True if it is minus one instead of plus one.
        vector<double> a_alt_angles // ALT stage angles at measurement points.
    );
    ~Altmodel_Inversion();

    // Overwrite forward function.
    int fwd(
        double *state, // State vector.
        double *signal, // Modelled signal (outout).
        double *jacobian // Modelled jacobian (outout).
    ) override;

    // Translate state function into information for the theodolite.
    int translate(
        double *state, // State vector.
        Quaternion *base, // Base (ALT=0) orientation quaternion.
        Vector *rotax // Rotation axis.
    );

    private:
    // Mapping parameters for fourth and fifth state parameters.
    size_t ivec_fixed;
    size_t ivec_rotax_first;
    size_t ivec_rotax_second;
    double val_fixed;
    // Abscissa.
    vector<double> alt_angles;
    // Copy of number of viewports.
    size_t dim_vp;

};


#endif
