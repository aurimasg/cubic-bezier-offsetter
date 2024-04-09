
#pragma once


#include "FloatPoint.h"


/**
 * Describes type of intersection between two line segments.
 */
enum class LineIntersectionKind : uint8_t {

    /**
     * No intersections found. Line segments are either zero length or they
     * are collinear.
     */
    None = 0,


    /**
     * Intersection was found within line segments.
     */
    Bounded,


    /**
     * Intersection was found, but beyong line segments.
     */
    Unbounded
};


/**
 * Describes intersection between two line segments.
 */
struct LineIntersection final {

    constexpr LineIntersection() {
    }


    constexpr LineIntersection(const LineIntersectionKind kind,
        const FloatPoint &intersectionPoint)
    :   Kind(kind),
        IntersectionPoint(intersectionPoint)
    {
    }


    /**
     * What kind of intersection was found.
     */
    LineIntersectionKind Kind = LineIntersectionKind::None;


    /**
     * Intersection point. Only valid if intersection kind is something other
     * than None.
     */
    FloatPoint IntersectionPoint;
};


/**
 * Describes intersection between two line segments.
 */
struct LineIntersectionSimple final {

    constexpr LineIntersectionSimple() {
    }


    constexpr LineIntersectionSimple(const FloatPoint &intersectionPoint)
    :   Found(true),
        IntersectionPoint(intersectionPoint)
    {
    }


    /**
     * Indicates wether intersection point was found.
     */
    bool Found = false;


    /**
     * Intersection point. Only valid if intersection kind is something other
     * than None.
     */
    FloatPoint IntersectionPoint;
};


/**
 * Represents a line segment defined by two points with floating-point
 * components.
 */
class FloatLine final {
public:

    /**
     * Constructs a line segment with both points set as zero points.
     */
    constexpr FloatLine() {
    }


    /**
     * Constructs a line segment between two given points.
     */
    constexpr FloatLine(const FloatPoint &p0, const FloatPoint &p1)
    :   P0(p0),
        P1(p1)
    {
    }


    /**
     * Constructs a line segment between two given points.
     */
    constexpr FloatLine(const double x0, const double y0, const double x1,
        const double y1)
    :   P0(x0, y0),
        P1(x1, y1)
    {
    }


    /**
     * Constructs a line segment between two given points.
     */
    constexpr FloatLine(const double x0, const double y0, const FloatPoint &p1)
    :   P0(x0, y0),
        P1(p1)
    {
    }


    /**
     * Constructs a line segment between two given points.
     */
    constexpr FloatLine(const FloatPoint &p0, const double x1, const double y1)
    :   P0(p0),
        P1(x1, y1)
    {
    }

public:

    /**
     * Return X coordinate of the first point.
     */
    double X0() const;


    /**
     * Return Y coordinate of the first point.
     */
    double Y0() const;


    /**
     * Return X coordinate of the second point.
     */
    double X1() const;


    /**
     * Return Y coordinate of the second point.
     */
    double Y1() const;


    /**
     * Returns difference between X component of point 2 and point 1 of this
     * line.
     */
    double Dx() const;


    /**
     * Returns difference between Y component of point 2 and point 1 of this
     * line.
     */
    double Dy() const;


    /**
     * Returns true if both points of this line are at the same place within
     * given error.
     */
    bool IsPoint(const double epsilon = DBL_EPSILON) const;


    /**
     * Returns a positive number of degrees from this linne to a given line.
     * Returned value is always positive and never exceeds 180.
     */
    double GetDegreesToLine(const FloatLine &l) const;


    /**
     * Returns a positive number of radians from this linne to a given line.
     * Returned value is always positive and never exceeds M_PI.
     */
    double GetRadiansToLine(const FloatLine &l) const;


    /**
     * Returns angle of this line in degrees. Angle is measured
     * counter-clock-wise from x axis.
     */
    double Angle() const;


    /**
     * Returns length of this line.
     */
    double Length() const;


    /**
     * Returns squared length of this line.
     */
    double LengthSquared() const;


    /**
     * Returns FloatLine with points reversed.
     */
    FloatLine Reversed() const;


    /**
     * Returns normalized vector which points to the same direction as this
     * line segment.
     */
    FloatPoint UnitVector() const;


    /**
     * Returns normal vector of this line. If both points of this line are in
     * the same place, returned vector will have both components set to zero.
     */
    FloatPoint NormalVector() const;


    /**
     * Returns normal vector of this line. Before returning, result is
     * normalized and will have length of 1. If both points of this line are
     * in the same place, returned vector will have both components set to
     * zero and its length will be zero.
     */
    FloatPoint UnitNormalVector() const;


    /**
     * Returns line intersection between this and other line segment.
     */
    LineIntersection Intersect(const FloatLine &l) const;


    /**
     * Returns simple line intersection between this and other line. Use this
     * method if you are not interested in finding out if intersection is
     * within line segment or not.
     */
    LineIntersectionSimple IntersectSimple(const FloatLine &l) const;


    /**
     * Returns a copy of this line with both points translated by a given
     * amount.
     */
    FloatLine Translated(const FloatPoint &p) const;


    /**
     * Returns point which is in the middle of this line segment.
     */
    FloatPoint MidPoint() const;


    /**
     * Extend this line segment by a given length. The first point of the line
     * will remain the same while the second point will be moved to the
     * direction of this line by a given distance. If points of this line
     * segment are equal or length value is zero, this method does nothing.
     *
     * @param length The distance by which to push the second point of this
     * line segment. Negative values are allowed.
     */
    void ExtendByLengthFront(const double length);


    /**
     * Extend this line segment by a given length. The second point of the
     * line will remain the same while the first point will be moved to the
     * direction opposite to the direction of this line by a given distance.
     * If points of this line segment are equal or length value is zero, this
     * method does nothing.
     *
     * @param length The distance by which to push the first point of this
     * line segment. Negative values are allowed.
     */
    void ExtendByLengthBack(const double length);


    /**
     * Returns true if a given point lies on this line segment.
     *
     * @param point Point to test.
     *
     * @param epsilon Tolerance for number comparison.
     */
    bool IsPointOnLineSegment(const FloatPoint &point,
        const double epsilon = DBL_EPSILON) const;


    /**
     * Returns true if a given point lies on this line.
     *
     * @param point Point to test.
     *
     * @param epsilon Tolerance for number comparison.
     */
    bool IsPointOnLine(const FloatPoint &point,
        const double epsilon = DBL_EPSILON) const;

public:
    FloatPoint P0;
    FloatPoint P1;
};


FORCE_INLINE double FloatLine::X0() const {
    return P0.X;
}


FORCE_INLINE double FloatLine::Y0() const {
    return P0.Y;
}


FORCE_INLINE double FloatLine::X1() const {
    return P1.X;
}


FORCE_INLINE double FloatLine::Y1() const {
    return P1.Y;
}


FORCE_INLINE double FloatLine::Dx() const {
    return P1.X - P0.X;
}


FORCE_INLINE double FloatLine::Dy() const {
    return P1.Y - P0.Y;
}


FORCE_INLINE FloatLine FloatLine::Reversed() const {
    return FloatLine(P1, P0);
}


FORCE_INLINE FloatPoint FloatLine::UnitVector() const {
    return (P1 - P0).UnitVector();
}


FORCE_INLINE FloatPoint FloatLine::NormalVector() const {
    return FloatPoint(Dy(), -Dx());
}


FORCE_INLINE FloatPoint FloatLine::UnitNormalVector() const {
    return NormalVector().UnitVector();
}


FORCE_INLINE FloatLine FloatLine::Translated(const FloatPoint &p) const {
    return FloatLine(P0 + p, P1 + p);
}


FORCE_INLINE FloatPoint FloatLine::MidPoint() const {
    return FloatPoint((P0.X + P1.X) * 0.5, (P0.Y + P1.Y) * 0.5);
}
