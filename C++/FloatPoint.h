
#pragma once


#include "Utils.h"


/**
 * Value returned when trying to decide orientation of 3 points.
 */
enum class TrianglePointOrientation : uint8_t {

    /**
     * Points are in clockwise orientation.
     */
    Clockwise,


    /**
     * Points are in counter-clockwise orientation.
     */
    CounterClockwise,


    /**
     * Points are collinear and orientation cannot be determined.
     */
    Collinear
};


/**
 * Represents a point in a two-dimensional plane defined as two floating point
 * values X and Y.
 */
class FloatPoint final {
public:

    /**
     * Constructs point with X and Y values set to zero.
     */
    constexpr FloatPoint() {
    }


    /**
     * Constructs point with given X and Y values.
     */
    constexpr FloatPoint(const double x, const double y)
    :   X(x),
        Y(y)
    {
    }


    /**
     * Constructs point with X and Y set to a given value.
     */
    constexpr explicit FloatPoint(const double value)
    :   X(value),
        Y(value)
    {
    }

public:

    /**
     * Returns greater than zero if three given points are clockwise. Returns
     * less than zero if points are counter-clockwise. And returns zero if
     * points are collinear.
     */
    static double Turn(const FloatPoint &p0, const FloatPoint &p1,
        const FloatPoint &p2);


    /**
     * Determines orientation of triangle defined by three given points.
     */
    static TrianglePointOrientation DetermineTriangleOrientation(
        const FloatPoint &p0, const FloatPoint &p1, const FloatPoint &p2);


    /**
     * Returns true if triangle defined by three given points is clockwise.
     * Returns false if triangle is counter-clockwise or if points are
     * collinear.
     */
    static bool IsTriangleClockwise(const FloatPoint &p0,
        const FloatPoint &p1, const FloatPoint &p2);

public:

    /**
     * Returns true if this point is equal to another point.
     *
     * @param point Another point to compare with.
     *
     * @param epsilon Comparison epsilon. This value will be used when
     * comparing both X and Y components separately. Must be at least zero.
     */
    bool IsEqual(const FloatPoint &point, const double epsilon = DBL_EPSILON) const;


    /**
     * Returns distance from this point to another point.
     */
    double DistanceTo(const FloatPoint &point) const;


    /**
     * Returns squared distance from this point to another point.
     */
    double DistanceToSquared(const FloatPoint &point) const;


    /**
     * Returns length of vector defined by X and Y components of this point.
     */
    double Length() const;


    /**
     * Returns squared length of vector defined by X and Y components of this
     * point.
     */
    double LengthSquared() const;


    /**
     * Returns normalized version of this vector. If this vector has length of
     * zero, vector with both components set to zero will be returned.
     */
    FloatPoint UnitVector() const;


    /**
     * Returns vector which has direction perpendicular to the direction of
     * this vector.
     */
    FloatPoint NormalVector() const;


    /**
     * Returns vector which has direction perpendicular to the direction of
     * this vector. Before returning, resulting vector is normalized.
     */
    FloatPoint UnitNormalVector() const;


    /**
     * Returns cross product of two 2D vectors (this × point).
     *
     * Since both 2D vectors lie on the same XY plane, the only meaningful
     * return value is Z component of cross product. This method returns that
     * and does not calculate anything else.
     */
    double Cross(const FloatPoint &point) const;


    /**
     * Returns dot product of this vector and a given vector.
     */
    double Dot(const FloatPoint &point) const;


    /**
     * Rotates this vector by 90 degrees counter-clockwise and return rotated
     * version of this vector.
     */
    FloatPoint Rotated90CCW() const;

public:
    void operator+=(const FloatPoint &p);
    void operator-=(const FloatPoint &p);
    void operator*=(const double s);
    void operator*=(const FloatPoint &p);
    void operator/=(const double s);
    bool operator==(const FloatPoint &p) const;
    bool operator!=(const FloatPoint &p) const;
public:
    double X = 0;
    double Y = 0;
};


FORCE_INLINE FloatPoint operator+(const FloatPoint &p0, const FloatPoint &p1) {
    return FloatPoint(p0.X + p1.X, p0.Y + p1.Y);
}


FORCE_INLINE FloatPoint operator-(const FloatPoint &p0, const FloatPoint &p1) {
    return FloatPoint(p0.X - p1.X, p0.Y - p1.Y);
}


FORCE_INLINE FloatPoint operator*(const FloatPoint &p, const double s) {
    return FloatPoint(p.X * s, p.Y * s);
}


FORCE_INLINE FloatPoint operator*(const double s, const FloatPoint &p) {
    return FloatPoint(p.X * s, p.Y * s);
}


FORCE_INLINE FloatPoint operator*(const FloatPoint &p0, const FloatPoint &p1) {
    return FloatPoint(p0.X * p1.X, p0.Y * p1.Y);
}


FORCE_INLINE FloatPoint operator/(const FloatPoint &p, const double s) {
    return FloatPoint(p.X / s, p.Y / s);
}


FORCE_INLINE FloatPoint operator/(const FloatPoint &p0, const FloatPoint &p1) {
    return FloatPoint(p0.X / p1.X, p0.Y / p1.Y);
}


FORCE_INLINE FloatPoint operator-(const FloatPoint &p) {
    return FloatPoint(-p.X, -p.Y);
}


FORCE_INLINE double FloatPoint::Turn(const FloatPoint &point1,
    const FloatPoint &point2, const FloatPoint &point3)
{
    return (point2 - point1).Cross(point3 - point1);
}


FORCE_INLINE bool FloatPoint::IsTriangleClockwise(const FloatPoint &p0, const FloatPoint &p1, const FloatPoint &p2) {
    return DetermineTriangleOrientation(p0, p1, p2) == TrianglePointOrientation::Clockwise;
}


FORCE_INLINE double FloatPoint::DistanceTo(const FloatPoint &point) const {
    return (point - *this).Length();
}


FORCE_INLINE double FloatPoint::DistanceToSquared(const FloatPoint &point) const {
    return (point - *this).LengthSquared();
}


FORCE_INLINE double FloatPoint::Length() const {
    return Sqrt(LengthSquared());
}


FORCE_INLINE double FloatPoint::LengthSquared() const {
    return X * X + Y * Y;
}


FORCE_INLINE FloatPoint FloatPoint::NormalVector() const {
    return FloatPoint(Y, -X);
}


FORCE_INLINE FloatPoint FloatPoint::UnitNormalVector() const {
    return NormalVector().UnitVector();
}


FORCE_INLINE double FloatPoint::Cross(const FloatPoint &point) const {
    return (X * point.Y) - (Y * point.X);
}


FORCE_INLINE double FloatPoint::Dot(const FloatPoint &point) const {
    return (X * point.X) + (Y * point.Y);
}


FORCE_INLINE FloatPoint FloatPoint::Rotated90CCW() const {
    return FloatPoint(Y, -X);
}


FORCE_INLINE void FloatPoint::operator+=(const FloatPoint &p) {
    X += p.X;
    Y += p.Y;
}


FORCE_INLINE void FloatPoint::operator-=(const FloatPoint &p) {
    X -= p.X;
    Y -= p.Y;
}


FORCE_INLINE void FloatPoint::operator*=(const double s) {
    X *= s;
    Y *= s;
}


FORCE_INLINE void FloatPoint::operator*=(const FloatPoint &p) {
    X *= p.X;
    Y *= p.Y;
}


FORCE_INLINE void FloatPoint::operator/=(const double s) {
    X /= s;
    Y /= s;
}


FORCE_INLINE bool FloatPoint::operator==(const FloatPoint &p) const {
    return FuzzyIsEqual(X, p.X) and FuzzyIsEqual(Y, p.Y);
}


FORCE_INLINE bool FloatPoint::operator!=(const FloatPoint &p) const {
    return !FuzzyIsEqual(X, p.X) or !FuzzyIsEqual(Y, p.Y);
}
