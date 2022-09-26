
#include "FloatPoint.h"


TrianglePointOrientation FloatPoint::DetermineTriangleOrientation(
    const FloatPoint &p1, const FloatPoint &p2, const FloatPoint &p3)
{
    const double turn = Turn(p1, p2, p3);

    if (FuzzyIsZero(turn)) {
        return TrianglePointOrientation::Collinear;
    } else if (turn > 0.0) {
        return TrianglePointOrientation::Clockwise;
    }

    return TrianglePointOrientation::CounterClockwise;
}


bool FloatPoint::IsEqual(const FloatPoint &point, const double epsilon) const
{
    return IsEqualWithEpsilon(X, point.X, epsilon) and
        IsEqualWithEpsilon(Y, point.Y, epsilon);
}


FloatPoint FloatPoint::UnitVector() const
{
    const double mag2 = LengthSquared();

    if (mag2 != 0.0 and mag2 != 1.0) {
        const double length = Sqrt(mag2);

        return FloatPoint(X / length, Y / length);
    }

    return *this;
}
