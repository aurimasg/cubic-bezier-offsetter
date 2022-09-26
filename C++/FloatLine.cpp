
#include "FloatLine.h"


bool FloatLine::IsPoint(const double epsilon) const
{
    return (IsEqualWithEpsilon(P1.X, P2.X, epsilon) and
        IsEqualWithEpsilon(P1.Y, P2.Y, epsilon));
}


double FloatLine::GetDegreesToLine(const FloatLine &l) const
{
    return Rad2Deg(GetRadiansToLine(l));
}


double FloatLine::GetRadiansToLine(const FloatLine &l) const
{
    if (IsPoint() or l.IsPoint()) {
        return 0;
    }

    const double c = (Dx() * l.Dx() + Dy() * l.Dy()) / (Length() * l.Length());

    // FLT_EPSILON instead of DBL_EPSILON is used deliberately.
    static constexpr double kMinRange = -1.0 - FLT_EPSILON;
    static constexpr double kMaxRange = 1.0 + FLT_EPSILON;

    // Return 0 instead of PI if c is outside range.
    if (c >= kMinRange and c <= kMaxRange) {
        return Acos(Clamp(c, -1.0, 1.0));
    }

    return 0;
}


double FloatLine::Angle() const
{
    const double dx = P2.X - P1.X;
    const double dy = P2.Y - P1.Y;
    const double theta = Rad2Deg(Atan2(-dy, dx));
    const double thetaNormalized = theta < 0 ? theta + 360 : theta;

    if (FuzzyIsEqual(thetaNormalized, 360.0)) {
        return 0;
    }

    if (FuzzyIsZero(thetaNormalized)) {
        // In case we have -0, return positive zero.
        return 0;
    }

    return thetaNormalized;
}


double FloatLine::Length() const
{
    const double x = P2.X - P1.X;
    const double y = P2.Y - P1.Y;

    return Sqrt(x * x + y * y);
}


double FloatLine::LengthSquared() const
{
    const double x = P2.X - P1.X;
    const double y = P2.Y - P1.Y;

    return x * x + y * y;
}


LineIntersection FloatLine::Intersect(const FloatLine &l) const
{
    const FloatPoint a = P2 - P1;
    const FloatPoint b = l.P1 - l.P2;
    const double denominator = a.Y * b.X - a.X * b.Y;

    if (denominator == 0) {
        return LineIntersection();
    }

    const FloatPoint c = P1 - l.P1;
    const double reciprocal = 1.0 / denominator;
    const double na = (b.Y * c.X - b.X * c.Y) * reciprocal;

    const FloatPoint point = P1 + a * na;

    if (na < 0 or na > 1) {
        return LineIntersection(LineIntersectionKind::Unbounded, point);
    }

    const double nb = (a.X * c.Y - a.Y * c.X) * reciprocal;

    if (nb < 0 or nb > 1) {
        return LineIntersection(LineIntersectionKind::Unbounded, point);
    }

    return LineIntersection(LineIntersectionKind::Bounded, point);
}


LineIntersectionSimple FloatLine::IntersectSimple(const FloatLine &l) const
{
    const FloatPoint a = P2 - P1;
    const FloatPoint b = l.P1 - l.P2;
    const double denominator = a.Y * b.X - a.X * b.Y;

    if (denominator == 0) {
        return LineIntersectionSimple();
    }

    const FloatPoint c = P1 - l.P1;
    const double reciprocal = 1.0 / denominator;
    const double na = (b.Y * c.X - b.X * c.Y) * reciprocal;

    return LineIntersectionSimple(P1 + a * na);
}


void FloatLine::ExtendByLengthFront(const double length)
{
    if (IsPoint() or FuzzyIsZero(length)) {
        return;
    }

    const FloatPoint v = UnitVector();

    P2 = FloatPoint(P2.X + (v.X * length), P2.Y + (v.Y * length));
}


void FloatLine::ExtendByLengthBack(const double length)
{
    if (IsPoint() or FuzzyIsZero(length)) {
        return;
    }

    const FloatPoint v = UnitVector();

    P1 = FloatPoint(P1.X - (v.X * length), P1.Y - (v.Y * length));
}


bool FloatLine::IsPointOnLineSegment(const FloatPoint &point,
    const double epsilon) const
{
    const double cross = (point.Y - P1.Y) * (P2.X - P1.X) -
        (point.X - P1.X) * (P2.Y - P1.Y);

    if (Abs(cross) > epsilon) {
        return false;
    }

    const double dot = (point.X - P1.X) * (P2.X - P1.X) +
        (point.Y - P1.Y) * (P2.Y - P1.Y);

    if (dot < 0) {
        return false;
    }

    const double sql = (P2.X - P1.X) * (P2.X - P1.X) +
        (P2.Y - P1.Y) * (P2.Y - P1.Y);

    return dot <= sql;
}


bool FloatLine::IsPointOnLine(const FloatPoint &point,
    const double epsilon) const
{
    const double cross = (point.Y - P1.Y) * (P2.X - P1.X) -
        (point.X - P1.X) * (P2.Y - P1.Y);

    return Abs(cross) <= epsilon;
}
