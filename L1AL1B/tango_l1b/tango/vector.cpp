// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "vector.h"

// Constructors.
Vector::Vector() {
    v.resize(DIM_VEC);
    for (size_t ivec=0 ; ivec<DIM_VEC ; ivec++) v[ivec] = 0.0;
}
Vector::Vector(double *input) {
    v.resize(DIM_VEC);
    memcpy(v.data(),input,DIM_VEC*sizeof(double));
}
Vector::Vector(double x,double y, double z) {
    v.resize(DIM_VEC);
    v[0] = x;
    v[1] = y;
    v[2] = z;
}
// Elements.
double& Vector::operator[](int idx) {
    return v[idx];
}
// Multiplication (vector and scalar, commutative).
Vector Vector::operator*(double rhs) {
    Vector res = *this;
    res *= rhs;
    return res;
}
Vector operator*(double lhs, Vector rhs) {
    return rhs * lhs;
}
Vector& Vector::operator*=(double rhs) {
    for (size_t ivec=0 ; ivec<DIM_VEC ; ivec++) v[ivec] *= rhs;
    return *this;
}
// Division (vector divided by scalar).
Vector Vector::operator/(double rhs) {
    Vector res = *this;
    res /= rhs;
    return res;
}
Vector& Vector::operator/=(double rhs) {
    for (size_t ivec=0 ; ivec<DIM_VEC ; ivec++) v[ivec] /= rhs;
    return *this;
}
// Addition, subtraction and unary minus sign.
Vector& Vector::operator+=(Vector rhs) {
    for (size_t ivec=0 ; ivec<DIM_VEC ; ivec++) v[ivec] += rhs.v[ivec];
    return *this;
}
Vector Vector::operator+(Vector rhs) {
    Vector res = *this;
    res += rhs;
    return res;
}
Vector Vector::operator-=(Vector rhs) {
    for (size_t ivec=0 ; ivec<DIM_VEC ; ivec++) v[ivec] -= rhs.v[ivec];
    return *this;
}
Vector Vector::operator-(Vector rhs) {
    Vector res = *this;
    res -= rhs;
    return res;
}
Vector Vector::operator-() {
    Vector res = *this;
    res *= -1.0;
    return res;
}
// Dot and cross product.
double Vector::dot(Vector other) {
    double res = 0.0;
    for (size_t ivec=0 ; ivec<DIM_VEC ; ivec++) res += v[ivec]*other[ivec];
    return res;
}
Vector Vector::cross(Vector other) {
    Vector res;
    for (size_t ivec=0 ; ivec<DIM_VEC ; ivec++) res[ivec] = v[(ivec+1)%DIM_VEC]*other[(ivec+2)%DIM_VEC] - v[(ivec+2)%DIM_VEC]*other[(ivec+1)%DIM_VEC];
    return res;
}
// Normalization.
void Vector::normalize() {
    double mag = sqrt(pow(v[0],2.0)+pow(v[1],2.0)+pow(v[2],2.0));
    if (mag != 0.0) this->operator/=(mag);
}
// Transform to double array.
void Vector::get(double *res) {
    memcpy(res,v.data(),DIM_VEC*sizeof(double));
}

// Constructors.
Quaternion::Quaternion() {
    v = Vector();
    r = 0.0;
}
Quaternion::Quaternion(double *input) {
    v = Vector(input);
    r = input[DIM_VEC];
}
Quaternion::Quaternion(Vector v, double r) {
    this->v = Vector(v);
    this->r = r;
}
Quaternion::Quaternion(double x, double y, double z, double r) {
    v = Vector(x,y,z);
    this->r = r;
}
// Elements
double& Quaternion::operator[](int idx) {
    return idx<(int)DIM_VEC?v[idx]:r;
}
// Multiplication (quaternion and scalar, commutative).
Quaternion Quaternion::operator*(double rhs) {
    Quaternion res = *this;
    res *= rhs;
    return res;
}
Quaternion operator*(double lhs, Quaternion rhs) {
    return rhs * lhs;
}
Quaternion& Quaternion::operator*=(double rhs) {
    v *= rhs;
    r *= rhs;
    return *this;
}
// Division (quaternion divided by scalar).
Quaternion Quaternion::operator/(double rhs) {
    Quaternion res = *this;
    res /= rhs;
    return res;
}
Quaternion& Quaternion::operator/=(double rhs) {
    v /= rhs;
    r /= rhs;
    return *this;
}
// Addition, subtraction and unary minus sign.
Quaternion& Quaternion::operator+=(Quaternion rhs) {
    v += rhs.v;
    r += rhs.r;
    return *this;
}
Quaternion Quaternion::operator+(Quaternion rhs) {
    Quaternion res = *this;
    res += rhs;
    return res;
}
Quaternion Quaternion::operator-=(Quaternion rhs) {
    v -= rhs.v;
    r -= rhs.r;
    return *this;
}
Quaternion Quaternion::operator-(Quaternion rhs) {
    Quaternion res = *this;
    res -= rhs;
    return res;
}
Quaternion Quaternion::operator-() {
    Quaternion res = *this;
    res *= -1.0;
    return res;
}
// Conjugate.
Quaternion Quaternion::conjug() {
    return Quaternion(-this->v,this->r);
}
// Rotation.
Vector Quaternion::quaternion_rotate(Vector orig) {
    return (pow(r,2.0)-v.dot(v))*orig + 2.0*r*v.cross(orig) + 2.0*v*v.dot(orig);
}
// Quaternion multiplication.
Quaternion Quaternion::operator*(Quaternion rhs) {
    return Quaternion(r*rhs.v+rhs.r*v+v.cross(rhs.v),r*rhs.r - v.dot(rhs.v));
}
// Quaternion division.
Quaternion Quaternion::operator/(Quaternion rhs) {
    return (*this)*rhs.conjug() / (pow(rhs.r,2.0)+rhs.v.dot(rhs.v));
}
// Raising a power.
Quaternion Quaternion::power(double p) {
    double norm = sqrt(pow(r,2.0) + v.dot(v));
    double c = r/norm;
    Vector dir = v/norm;
    double s = sqrt(dir.dot(dir));
    dir /= s;
    double mag = pow(norm,p);
    double ang = atan2(s,c) * p;
    return mag*Quaternion(dir*sin(ang),cos(ang));
}
// Quaternion interpolation.
Quaternion Quaternion::interpolate(Quaternion next, double w) {
    Quaternion diff = next / (*this);
    return diff.power(w) * (*this);
}
// Normalization.
void Quaternion::normalize() {
    double mag = sqrt(pow(r,2.0) + v.dot(v));
    *this /= mag;
}
// Make the real part positive by eventually negating it.
void Quaternion::intuitive_representation() {
    if (r < 0) *this = -*this;
}
// Transform to double array.
void Quaternion::get(double *res) {
    v.get(res);
    res[3] = r;
}

