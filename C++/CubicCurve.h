
#pragma once


#include "FloatPoint.h"
#include "FloatLine.h"


/**
 * Represents a two dimensional cubic curve. It is defined by start point, end
 * point and two control points.
 */
class CubicCurve final {
public:

    /**
     * Constructs cubic curve with all coordinates set to zero.
     */
    constexpr CubicCurve() {
    }


    /**
     * Constructs cubic curve with given points.
     */
    constexpr CubicCurve(const double p0x, const double p0y, const double p1x,
        const double p1y, const double p2x, const double p2y,
        const double p3x, const double p3y)
    :   P0(p0x, p0y),
        P1(p1x, p1y),
        P2(p2x, p2y),
        P3(p3x, p3y)
    {
    }


    /**
     * Constructs cubic curve with given points.
     *
     * @param p0 Start point of curve.
     *
     * @param p1 First control point of curve.
     *
     * @param p2 Second control point of curve.
     *
     * @param p3 End point of curve.
     */
    constexpr CubicCurve(const FloatPoint &p0, const FloatPoint &p1,
        const FloatPoint &p2, const FloatPoint &p3)
    :   P0(p0),
        P1(p1),
        P2(p2),
        P3(p3)
    {
    }


    /**
     * Constructs cubic curve from quadratic curve parameters.
     */
    CubicCurve(const FloatPoint &p0, const FloatPoint &p1,
        const FloatPoint &p2);


    /**
     * Constructs cubic curve from a line. First control point of constructed
     * curve will be on line from p0 and p1 at 1/3rd of distance. Second
     * control point will be at 2/3rds of distance on line from p0 to p1.
     */
    CubicCurve(const FloatPoint &p0, const FloatPoint &p1);

public:

    /**
     * Returns true if all points of this curve are at the same place within
     * given tolerance.
     *
     * @param tolerance Tolerance value to test curve flatness against. Must
     * be greater than 0.0. Larger value means less accurate test.
     */
    bool IsPoint(const double tolerance = DBL_EPSILON) const;


    /**
     * Returns true if this curve is a straight line. Curve is straight line
     * when both control points lie on line segment between the first and the
     * last point of the curve.
     */
    bool IsStraight() const;


    /**
     * Returns start tangent of this curve. Start tangent is a line going from
     * P0 to P1, P2 or P3, depending which one is the first found to be not
     * equal to P0.
     *
     * @param epsilon Point comparison tolerance.
     */
    FloatLine StartTangent(const double epsilon = DBL_EPSILON) const;


    /**
     * Returns end tangent of this curve. Start tangent is a line going from
     * P3 to P2, P1 or P0, depending which one is the first found to be not
     * equal to P3.
     *
     * @param epsilon Point comparison tolerance.
     */
    FloatLine EndTangent(const double epsilon = DBL_EPSILON) const;


    /**
     * Returns point on curve at a given t.
     *
     * @param t Value from 0 to 1 to evaluate curve at.
     */
    FloatPoint PointAt(const double t) const;


    /**
     * Returns normal vector of this curve at a given t. Returned vector is
     * not normalized.
     *
     * @param t Value from 0 to 1 to evaluate curve at.
     */
    FloatPoint NormalVector(const double t) const;


    /**
     * Returns normal vector of this curve at a given t. Length of returned
     * vector is 1.
     *
     * @param t Value from 0 to 1 to evaluate curve at.
     */
    FloatPoint UnitNormalVector(const double t) const;


    /**
     * Returns first derivative at a given t.
     */
    FloatPoint DerivedAt(const double t) const;


    /**
     * Returns second derivative at a given t.
     */
    FloatPoint SecondDerivedAt(const double t) const;


    /**
     * Returns a sub-curve of this curve. If start and end positions are
     * equal, a curve with all its points set to the same point is returned.
     *
     * @param t0 Start position on this curve. Must be less than or equal to
     * t1.
     *
     * @param t1 End position on this curve. Must be greater than or equal to
     * t0.
     */
    CubicCurve GetSubcurve(const double t0, const double t1) const;


    /**
     * Finds points of maximum curvature and returns the number of points
     * found. Up to 3 points are returned.
     *
     * @param tValues Array of t values. Upon return, up to 3 values are set
     * in this array. So it must be large enough to hold 3 values. All values
     * written to this array are within 0-1 range, including 0 and 1.
     */
    int FindMaxCurvature(double tValues[3]) const;


    /**
     * Finds inflection points of this curve and returns a number of points
     * found.
     *
     * @param tValues Array of t values. Upon return, up to 2 values are set
     * in this array. So it must be large enough to hold 2 values. All values
     * written to this array are within 0-1 range.
     */
    int FindInflections(double tValues[2]) const;


    /**
     * Split this curve in half and put result into two output parameters.
     *
     * @param l Will receive the first half of this curve. The first point of
     * curve l will be the first point of this curve.
     *
     * @param r Will receive the second half of this curve. The last point of
     * curve r will be the last point of this curve.
     */
    void Split(CubicCurve &l, CubicCurve &r) const;


    /**
     * Finds intersections between this curve and line defined by two points.
     * Intersections found outside of line segment going from the first to the
     * second point will be included in the result. This method returns the
     * amount of intersections found. There can be up to 3 intersections.
     *
     * @param linePointA First point of line.
     *
     * @param linePointB Second point of line.
     *
     * @param t Array where up to 3 values will be stored. These values are
     * positions on this curve where intersections are found. All of them are
     * in range from 0 to 1 including 0 and 1.
     */
    int FindRayIntersections(const FloatPoint &linePointA,
        const FloatPoint &linePointB, double t[3]) const;

public:

    /**
     * Start point of the curve.
     */
    FloatPoint P0;


    /**
     * First control point of the curve.
     */
    FloatPoint P1;


    /**
     * Second control point of the curve.
     */
    FloatPoint P2;


    /**
     * End point.
     */
    FloatPoint P3;
};
