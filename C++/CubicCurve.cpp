
#include "CubicCurve.h"
#include "Quicksort.h"


CubicCurve::CubicCurve(const FloatPoint &p0, const FloatPoint &p1, const FloatPoint &p2)
:   P0(p0),
    P3(p2)
{
    P1.X = p0.X + (2.0 / 3.0) * (p1.X - p0.X);
    P1.Y = p0.Y + (2.0 / 3.0) * (p1.Y - p0.Y);
    P2.X = p1.X + (1.0 / 3.0) * (p2.X - p1.X);
    P2.Y = p1.Y + (1.0 / 3.0) * (p2.Y - p1.Y);
}


CubicCurve::CubicCurve(const FloatPoint &p0, const FloatPoint &p1)
:   P0(p0),
    P3(p1)
{
    P1.X = p0.X + (1.0 / 3.0) * (p1.X - p0.X);
    P1.Y = p0.Y + (1.0 / 3.0) * (p1.Y - p0.Y);
    P2.X = p0.X + (2.0 / 3.0) * (p1.X - p0.X);
    P2.Y = p0.Y + (2.0 / 3.0) * (p1.Y - p0.Y);
}


bool CubicCurve::IsPoint(const double tolerance) const
{
    return P0.IsEqual(P3, tolerance) and
        P0.IsEqual(P1, tolerance) and
        P0.IsEqual(P2, tolerance);
}


bool CubicCurve::IsStraight() const
{
    const double minx = Min(P0.X, P3.X);
    const double miny = Min(P0.Y, P3.Y);
    const double maxx = Max(P0.X, P3.X);
    const double maxy = Max(P0.Y, P3.Y);

    return
        // Is P1 located between P0 and P3?
        minx <= P1.X and
        miny <= P1.Y and
        maxx >= P1.X and
        maxy >= P1.Y and
        // Is P2 located between P0 and P3?
        minx <= P2.X and
        miny <= P2.Y and
        maxx >= P2.X and
        maxy >= P2.Y and
        // Are all points collinear?
        FuzzyIsZero(FloatPoint::Turn(P0, P1, P3)) and
        FuzzyIsZero(FloatPoint::Turn(P0, P2, P3));
}


FloatLine CubicCurve::StartTangent(const double epsilon) const
{
    if (P0.IsEqual(P1, epsilon)) {
        if (P0.IsEqual(P2, epsilon)) {
            return FloatLine(P0, P3);
        }

        return FloatLine(P0, P2);
    }

    return FloatLine(P0, P1);
}


FloatLine CubicCurve::EndTangent(const double epsilon) const
{
    if (P3.IsEqual(P2, epsilon)) {
        if (P3.IsEqual(P1, epsilon)) {
            return FloatLine(P3, P0);
        }

        return FloatLine(P3, P1);
    }

    return FloatLine(P3, P2);
}


FloatPoint CubicCurve::PointAt(const double t) const
{
    const double it = 1.0 - t;

    const FloatPoint a0 = P0 * it + P1 * t;
    const FloatPoint b0 = P1 * it + P2 * t;
    const FloatPoint c0 = P2 * it + P3 * t;

    const FloatPoint a1 = a0 * it + b0 * t;
    const FloatPoint b1 = b0 * it + c0 * t;

    return a1 * it + b1 * t;
}


FloatPoint CubicCurve::NormalVector(const double t) const
{
    if (FuzzyIsZero(t)) {
        if (P0 == P1) {
            if (P0 == P2) {
                return FloatLine(P0, P3).NormalVector();
            } else {
                return FloatLine(P0, P2).NormalVector();
            }
        }
    } else if (FuzzyIsEqual(t, 1.0)) {
        if (P2 == P3) {
            if (P1 == P3) {
                return FloatLine(P0, P3).NormalVector();
            } else {
                return FloatLine(P1, P3).NormalVector();
            }
        }
    }

    const FloatPoint d = DerivedAt(t);

    return FloatPoint(d.Y, -d.X);
}


FloatPoint CubicCurve::UnitNormalVector(const double t) const
{
    const FloatPoint n = NormalVector(t);

    return n.UnitVector();
}


FloatPoint CubicCurve::DerivedAt(const double t) const
{
    const double it = 1.0 - t;
    const double d = t * t;
    const double a = -it * it;
    const double b = 1.0 - 4.0 * t + 3.0 * d;
    const double c = 2.0 * t - 3.0 * d;

    return 3.0 * FloatPoint(a * P0.X + b * P1.X + c * P2.X + d * P3.X,
        a * P0.Y + b * P1.Y + c * P2.Y + d * P3.Y);
}


FloatPoint CubicCurve::SecondDerivedAt(const double t) const
{
    const double a = 2.0 - 2.0 * t;
    const double b = -4.0 + 6.0 * t;
    const double c = 2.0 - 6.0 * t;
    const double d = 2.0 * t;

    return 3.0 * FloatPoint(a * P0.X + b * P1.X + c * P2.X + d * P3.X,
        a * P0.Y + b * P1.Y + c * P2.Y + d * P3.Y);
}


CubicCurve CubicCurve::GetSubcurve(const double t0, const double t1) const
{
    ASSERT(t0 <= t1);

    if (FuzzyIsEqual(t0, t1)) {
        const FloatPoint p = PointAt(t0);

        return CubicCurve(p, p, p, p);
    }

    if (t0 <= DBL_EPSILON) {
        if (t1 >= (1.0 - DBL_EPSILON)) {
            // Both parameters are 0.0 to 1.0.
            return *this;
        }

        // Cut at t1 only.
        const FloatPoint ab = InterpolateLinear(P0, P1, t1);
        const FloatPoint bc = InterpolateLinear(P1, P2, t1);
        const FloatPoint cd = InterpolateLinear(P2, P3, t1);
        const FloatPoint abc = InterpolateLinear(ab, bc, t1);
        const FloatPoint bcd = InterpolateLinear(bc, cd, t1);
        const FloatPoint abcd = InterpolateLinear(abc, bcd, t1);

        return CubicCurve(P0, ab, abc, abcd);
    }

    if (t1 >= (1.0 - DBL_EPSILON)) {
        // Cut at t0 only.
        const FloatPoint ab = InterpolateLinear(P0, P1, t0);
        const FloatPoint bc = InterpolateLinear(P1, P2, t0);
        const FloatPoint cd = InterpolateLinear(P2, P3, t0);
        const FloatPoint abc = InterpolateLinear(ab, bc, t0);
        const FloatPoint bcd = InterpolateLinear(bc, cd, t0);
        const FloatPoint abcd = InterpolateLinear(abc, bcd, t0);

        return CubicCurve(abcd, bcd, cd, P3);
    }

    // Cut at both t0 and t1.
    const FloatPoint ab0 = InterpolateLinear(P0, P1, t1);
    const FloatPoint bc0 = InterpolateLinear(P1, P2, t1);
    const FloatPoint cd0 = InterpolateLinear(P2, P3, t1);
    const FloatPoint abc0 = InterpolateLinear(ab0, bc0, t1);
    const FloatPoint bcd0 = InterpolateLinear(bc0, cd0, t1);
    const FloatPoint abcd0 = InterpolateLinear(abc0, bcd0, t1);

    const double m = t0 / t1;

    const FloatPoint ab1 = InterpolateLinear(P0, ab0, m);
    const FloatPoint bc1 = InterpolateLinear(ab0, abc0, m);
    const FloatPoint cd1 = InterpolateLinear(abc0, abcd0, m);
    const FloatPoint abc1 = InterpolateLinear(ab1, bc1, m);
    const FloatPoint bcd1 = InterpolateLinear(bc1, cd1, m);
    const FloatPoint abcd1 = InterpolateLinear(abc1, bcd1, m);

    return CubicCurve(abcd1, bcd1, cd1, abcd0);
}


static int AcceptRoot(double *t, const double root)
{
    if (root < -DBL_EPSILON) {
        return 0;
    } else if (root > (1.0 + DBL_EPSILON)) {
        return 0;
    }

    t[0] = Clamp(root, 0.0, 1.0);

    return 1;
}


/*
 * Roots must not be nullptr. Returns 0, 1 or 2.
 */
static int FindQuadraticRoots(const double a, const double b, const double c,
    double roots[2])
{
    ASSERT(roots != nullptr);

    const double delta = b * b - 4.0 * a * c;

    if (delta < 0.0) {
        return 0;
    }

    if (delta > 0.0) {
        const double d = Sqrt(delta);
        const double q = -0.5 * (b + (b < 0.0 ? -d : d));
        const double rv0 = q / a;
        const double rv1 = c / q;

        if (FuzzyIsEqual(rv0, rv1)) {
            return AcceptRoot(roots, rv0);
        }

        if (rv0 < rv1) {
            int n = AcceptRoot(roots, rv0);

            n += AcceptRoot(roots + n, rv1);

            return n;
        } else {
            int n = AcceptRoot(roots, rv1);

            n += AcceptRoot(roots + n, rv0);

            return n;
        }
    }

    if (a != 0) {
        return AcceptRoot(roots, -0.5 * b / a);
    }

    return 0;
}


static double CubeRoot(const double x)
{
    return pow(x, 1.0 / 3.0);
}


static bool DoubleArrayContainsValue(const double *array, const int count,
    const double value)
{
    for (int i = 0; i < count; i++) {
        if (FuzzyIsEqual(array[i], value)) {
            return true;
        }
    }

    return false;
}


static int DeduplicateDoubleArray(double *array, const int currentCount)
{
    int newCount = 0;

    for (int i = 0; i < currentCount; i++) {
        const double value = array[i];

        if (DoubleArrayContainsValue(array, newCount, value)) {
            continue;
        }

        array[newCount++] = value;
    }

    return newCount;
}


/**
 * This function is based on Numerical Recipes.
 * 5.6 Quadratic and Cubic Equations.
 */
static int FindCubicRoots(const double coe0, const double coe1, const double coe2,
    const double coe3, double roots[3])
{
    ASSERT(roots != nullptr);

    if (FuzzyIsZero(coe0)) {
        return FindQuadraticRoots(coe1, coe2, coe3, roots);
    }

    const double inva = 1.0 / coe0;

    const double a = coe1 * inva;
    const double b = coe2 * inva;
    const double c = coe3 * inva;

    const double Q = (a * a - b * 3.0) / 9.0;
    const double R = (2.0 * a * a * a - 9.0 * a * b + 27.0 * c) / 54.0;

    const double R2 = R * R;
    const double Q3 = Q * Q * Q;
    const double R2subQ3 = R2 - Q3;
    const double adiv3 = a / 3.0;

    if (R2subQ3 < 0.0) {
        // If Q and R are real (always true when a, b, c are real) and R2 < Q3,
        // then the cubic equation has three real roots.
        const double theta = Acos(Clamp(R / Sqrt(Q3), -1.0, 1.0));
        const double negative2RootQ = -2.0 * Sqrt(Q);

        const double x0 = negative2RootQ * Cos(theta / 3.0) - adiv3;
        const double x1 = negative2RootQ * Cos((theta + 2.0 * M_PI) / 3.0) - adiv3;
        const double x3 = negative2RootQ * Cos((theta - 2.0 * M_PI) / 3.0) - adiv3;

        int n = AcceptRoot(roots, x0);
        n += AcceptRoot(roots + n, x1);
        n += AcceptRoot(roots + n, x3);

        Quicksort(roots, n, [](const double a, const double b) -> bool {
            return a < b;
        });

        return DeduplicateDoubleArray(roots, n);
    }

    double A = Abs(R) + Sqrt(R2subQ3);

    A = CubeRoot(A);

    if (R > 0.0) {
        A = -A;
    }

    if (A != 0.0) {
        A += Q / A;
    }

    return AcceptRoot(roots, A - adiv3);
}


int CubicCurve::FindMaxCurvature(double *tValues) const
{
    const double axx = P1.X - P0.X;
    const double bxx = P2.X - 2.0 * P1.X + P0.X;
    const double cxx = P3.X + 3.0 * (P1.X - P2.X) - P0.X;

    const double cox0 = cxx * cxx;
    const double cox1 = 3.0 * bxx * cxx;
    const double cox2 = 2.0 * bxx * bxx + cxx * axx;
    const double cox3 = axx * bxx;

    const double ayy = P1.Y - P0.Y;
    const double byy = P2.Y - 2.0 * P1.Y + P0.Y;
    const double cyy = P3.Y + 3.0 * (P1.Y - P2.Y) - P0.Y;

    const double coy0 = cyy * cyy;
    const double coy1 = 3.0 * byy * cyy;
    const double coy2 = 2.0 * byy * byy + cyy * ayy;
    const double coy3 = ayy * byy;

    const double coe0 = cox0 + coy0;
    const double coe1 = cox1 + coy1;
    const double coe2 = cox2 + coy2;
    const double coe3 = cox3 + coy3;

    return FindCubicRoots(coe0, coe1, coe2, coe3, tValues);
}


int CubicCurve::FindInflections(double *tValues) const
{
    const double ax = P1.X - P0.X;
    const double ay = P1.Y - P0.Y;
    const double bx = P2.X - 2.0 * P1.X + P0.X;
    const double by = P2.Y - 2.0 * P1.Y + P0.Y;
    const double cx = P3.X + 3.0 * (P1.X - P2.X) - P0.X;
    const double cy = P3.Y + 3.0 * (P1.Y - P2.Y) - P0.Y;

    return FindQuadraticRoots(bx * cy - by * cx, ax * cy - ay * cx,
        ax * by - ay * bx, tValues);
}


void CubicCurve::Split(CubicCurve &l, CubicCurve &r) const
{
    const FloatPoint c = (P1 + P2) * 0.5;

    const FloatPoint aP2 = (P0 + P1) * 0.5;
    const FloatPoint bP3 = (P2 + P3) * 0.5;
    const FloatPoint aP3 = (aP2 + c) * 0.5;
    const FloatPoint bP2 = (bP3 + c) * 0.5;
    const FloatPoint m = (aP3 + bP2) * 0.5;

    l.P0 = P0;
    l.P1 = aP2;
    l.P2 = aP3;
    l.P3 = m;

    r.P0 = m;
    r.P1 = bP2;
    r.P2 = bP3;
    r.P3 = P3;
}


int CubicCurve::FindRayIntersections(const FloatPoint &linePointA,
    const FloatPoint &linePointB, double *t) const
{
    const FloatPoint v = linePointB - linePointA;

    const double ax = (P0.Y - linePointA.Y) * v.X - (P0.X - linePointA.X) * v.Y;
    const double bx = (P1.Y - linePointA.Y) * v.X - (P1.X - linePointA.X) * v.Y;
    const double cx = (P2.Y - linePointA.Y) * v.X - (P2.X - linePointA.X) * v.Y;
    const double dx = (P3.Y - linePointA.Y) * v.X - (P3.X - linePointA.X) * v.Y;

    const double a = dx;
    const double b = cx * 3;
    const double c = bx * 3;

    const double D = ax;
    const double A = a - (D - c + b);
    const double B = b + (3 * D - 2 * c);
    const double C = c - (3 * D);

    return FindCubicRoots(A, B, C, D, t);
}
