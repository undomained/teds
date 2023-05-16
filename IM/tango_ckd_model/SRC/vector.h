// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef VECTOR_H
#define VECTOR_H

#include "header.h"

struct Vector {
    private:
    // Member.
    vector<double> v;
    public:
    // Constructors.
    Vector();
    Vector(double *input);
    Vector(double x,double y, double z);
    // Elements.
    double& operator[](int idx);
    // Multiplication (vector and scalar, commutative).
    Vector operator*(double rhs);
    friend Vector operator*(double lhs, Vector rhs);
    Vector& operator*=(double rhs);
    // Division (vector divided by scalar).
    Vector operator/(double rhs);
    Vector& operator/=(double rhs);
    // Addition, subtraction and unary minus sign.
    Vector& operator+=(Vector rhs);
    Vector operator+(Vector rhs);
    Vector operator-=(Vector rhs);
    Vector operator-(Vector rhs);
    Vector operator-();
    // Dot and cross product.
    double dot(Vector other);
    Vector cross(Vector other);
    // Normalization.
    void normalize();
    // Transform to double array.
    void get(double *res);
};
struct Quaternion {
    private:
    // Members.
    double r;
    Vector v;
    public:
    // Constructors.
    Quaternion();
    Quaternion(double *input);
    Quaternion(Vector v, double r);
    Quaternion(double x, double y, double z, double r);
    // Elements.
    double& operator[](int idx);
    // Multiplication (quaternion and scalar, commutative).
    Quaternion operator*(double rhs);
    friend Quaternion operator*(double lhs, Quaternion rhs);
    Quaternion& operator*=(double rhs);
    // Division (quaternion divided by scalar).
    Quaternion operator/(double rhs);
    Quaternion& operator/=(double rhs);
    // Addition, subtraction and unary minus sign.
    Quaternion& operator+=(Quaternion rhs);
    Quaternion operator+(Quaternion rhs);
    Quaternion operator-=(Quaternion rhs);
    Quaternion operator-(Quaternion rhs);
    Quaternion operator-();
    // Conjugate.
    Quaternion conjug();
    // Rotation.
    Vector quaternion_rotate(Vector orig);
    // Quaternion multiplication.
    Quaternion operator*(Quaternion rhs);
    // Quaternion division.
    Quaternion operator/(Quaternion rhs);
    // Raising a power.
    Quaternion power(double p);
    // Quaternion interpolation.
    Quaternion interpolate(Quaternion next, double w);
    // Normalization.
    void normalize();
    // Make the real part positive by eventually negating it.
    void intuitive_representation();
    // Transform to double array.
    void get(double *res);
};

#endif
