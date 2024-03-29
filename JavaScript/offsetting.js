
/**
 * A maximum number of iterations for searching closest point to cusp that has
 * the first derivative long enough for finding start or end points of
 * circular arc for cusp.
 *
 * Smaller value means faster search, but worse accuracy when handling
 * cusp-like points of the curve.
 */
const NearCuspPointSearchMaxIterationCount = 18;


/**
 * After an attempt to find offset curve is made, squared lengths of all edges
 * of polygon enclosing curve is calculated and added together. If this length
 * is equal to or less that this number, resulting curve will be discarded
 * immediadely without attempt to add it to the output.
 *
 * Smaller value means smaller curves will be accepted for output.
 */
const MaximumTinyCurvePolygonPerimeterSquared = 1e-7;


/**
 * If a good circular arc approximation of a curve is found, but its radius is
 * very close to offset amount, scaled arc can collapse to a point or almost a
 * point. This is epsilon for testing if arc is large enough. Arcs with radius
 * smaller than this value will not be added to the output.
 *
 * Smaller value means smaller arcs will be accepted for output.
 */
const MinimumArcRadius = 1e-8;


/**
 * An upper limit of arc radius. Circular arcs with calculated radius greater
 * than this value will not be considered as accepted approximations of curve
 * segments.
 */
const MaximumArcRadius = 1e+6;


/**
 * Offsetter does not attempt to find exact cusp locations and does not
 * consider cusp only to be where derivative vector length is exactly zero.
 *
 * Smaller values means that sharper curve edges are considered cusps.
 */
const CuspDerivativeLengthSquared = 1.5e-4;


/**
 * If X and Y components of all points are equal when compared with this
 * epsilon, curve is considered a point.
 */
const CurvePointClumpTestEpsilon = 1e-14;


/**
 * Epsilon used to compare coordinates of circular arc centers to see if they
 * can be merged into a single circular arc.
 */
const ArcCenterComparisonEpsilon = 1e-8;


/**
 * When testing if curve is almost straight, cross products of unit vectors
 * are calculated as follows
 *
 *     Turn1 = (P1 → P2) ⨯ (P1 → P4)
 *     Turn2 = (P2 → P3) ⨯ (P1 → P4)
 *
 * Where P1, P2, P3 and P4 are curve points and (X → Y) are unit vectors going
 * from X to a direction of Y.
 *
 * Then these values are compared with zero. If they both are close to zero,
 * curve is considered approximately straight. This is the epsilon used for
 * comparing Turn1 and Turn2 values to zero.
 *
 * Bigger value means less straight curves are considered approximately
 * straight.
 */
const ApproximatelyStraightCurveTestApsilon = 1e-5;


/**
 * The logic is the same as for ApproximatelyStraightCurveTestApsilon value.
 * This value is used to determine if curve is completely straight, not just
 * approximately straight.
 *
 * Bigger value means less straight curves are considered completely straight.
 * This value should be smaller than ApproximatelyStraightCurveTestApsilon.
 */
const CompletelyStraightCurveTestApsilon = 1e-15;


/**
 * A list of positions on curve for testing if circular arc approximation of a
 * curve is good enough. All values must be from 0 to 1. Positions 0 and 1
 * should not be included here because distances at 0 and 1 are either
 * guaranteed to be exactly on the curve being approximated (on the sides of
 * the input curve) or they are tested before using another methods (ray cast
 * from arc center to triangle incenter point).
 */
const ArcProbePositions = [
    0.2,
    0.4,
    0.6,
    0.8
];


/**
 * A list of positions on curve for testing if candidate offset curve is good
 * enough. From points on original curve at each of these positions a ray is
 * cast to normal direction and intersection is found with offset candidate.
 * Then distance is checked to see if it is within maximum error. 0 and 1
 * should not be added in this list because candidate points at 0 and 1 is
 * guaranteed to be at the right place already.
 *
 * Note that testing involves cubic root finding which is not very cheap
 * operation. Adding more testing points increases precision, but also
 * increases the time spent for testing if candidate is good.
 */
const SimpleOffsetProbePositions = [
    0.25,
    0.5,
    0.85
];


/**
 * Returns true if two given numbers are considered equal.
 */
function fuzzyIsEqual(a, b) {
    return Math.abs(a - b) < Number.EPSILON;
}


/**
 * Returns true if a number can be considered being equal to zero.
 */
function fuzzyIsZero(a) {
    return Math.abs(a) < Number.EPSILON;
}


/**
 * Returns true if two numbers differ not more than a given epsilon.
 */
function isEqualWithEpsilon(a, b, epsilon) {
    return Math.abs(a - b) < epsilon;
}


/**
 * Returns true if a number differs from zero by less than a given epsilon.
 */
function isZeroWithEpsilon(a, epsilon) {
    return Math.abs(a) < epsilon;
}


/**
 * Convert degrees to radians.
 */
function deg2Rad(x) {
    return x * (Math.PI / 180.0);
}


/**
 * Convert radians to degrees.
 */
function rad2Deg(x) {
    return x * (180.0 / Math.PI);
}


/**
 * Returns value clamped to range between minimum and maximum values.
 */
function clamp(val, min, max) {
    return val > max ? max : val < min ? min : val;
}


/**
 * Linearly interpolate between A and B.
 * If t is 0, returns A.
 * If t is 1, returns B.
 * If t is something else, returns value linearly interpolated between A and B.
 */
function interpolateLinear(A, B, t) {
    //ASSERT(t >= 0);
    //ASSERT(t <= 1);
    return A + ((B - A) * t);
}


/**
 * Finds the greatest of the three values.
 */
function max3(a, b, c) {
    return Math.max(a, Math.max(b, c));
}


/**
 * Finds the smallest of the three values.
 */
function min3(a, b, c) {
    return Math.min(a, Math.min(b, c));
}


/**
 * Finds the greatest of the four values.
 */
function max4(a, b, c, d) {
    return Math.max(a, Math.max(b, Math.max(c, d)));
}


/**
 * Finds the smallest of the four values.
 */
function min4(a, b, c, d) {
    return Math.min(a, Math.min(b, Math.min(c, d)));
}


/**
 * Returns array containing a single number if it is accepted as result of
 * root finder or empty array if anumber is not accepted. Root is accepted if
 * it is between -EPSILON to 1.0 + EPSILON, allowing small deviation from
 * range of 0.0 to 1.0. All returned values are clamped to 0.0 - 1.0 range.
 */
function acceptRoot(root) {
    if (root < -Number.EPSILON) {
        return [];
    } else if (root > (1.0 + Number.EPSILON)) {
        return [];
    }

    return [ clamp(root, 0.0, 1.0) ];
}


function findQuadraticRoots(a, b, c) {
    const delta = b * b - 4.0 * a * c;

    if (delta < 0) {
        return [];
    }

    if (delta > 0.0) {
        const d = Math.sqrt(delta);
        const q = -0.5 * (b + (b < 0.0 ? -d : d));
        const rv0 = q / a;
        const rv1 = c / q;

        if (fuzzyIsEqual(rv0, rv1)) {
            return acceptRoot(roots, rv0);
        }

        if (rv0 < rv1) {
            var roots = acceptRoot(rv0);

            roots.concat(acceptRoot(rv1));

            return roots;
        } else {
            var roots = acceptRoot(rv1);

            roots.concat(acceptRoot(rv0));

            return roots;
        }
    }

    if (a != 0) {
        return acceptRoot(-0.5 * b / a);
    }

    return [];
}


function arrayContainsFloat(array, value) {
    for (let i = 0; i < array.length; i++) {
        if (fuzzyIsEqual(array[i], value)) {
            return true;
        }
    }

    return false;
}


function deduplicateFloatArray(array) {
    var rv = [];

    for (let i = 0; i < array.length; i++) {
        if (!arrayContainsFloat(rv, array[i])) {
            rv.push(array[i]);
        }
    }

    return rv;
}


/**
 * This function is based on Numerical Recipes.
 * 5.6 Quadratic and Cubic Equations.
 */
function findCubicRoots(coe0, coe1, coe2, coe3) {
    if (fuzzyIsZero(coe0)) {
        return findQuadraticRoots(coe1, coe2, coe3);
    }

    const inva = 1.0 / coe0;

    const a = coe1 * inva;
    const b = coe2 * inva;
    const c = coe3 * inva;

    const Q = (a * a - b * 3.0) / 9.0;
    const R = (2.0 * a * a * a - 9.0 * a * b + 27.0 * c) / 54.0;

    const R2 = R * R;
    const Q3 = Q * Q * Q;
    const R2subQ3 = R2 - Q3;
    const adiv3 = a / 3.0;

    if (R2subQ3 < 0.0) {
        // If Q and R are real (always true when a, b, c are real) and R2 < Q3,
        // then the cubic equation has three real roots.
        const theta = Math.acos(clamp(R / Math.sqrt(Q3), -1.0, 1.0));
        const negative2RootQ = -2.0 * Math.sqrt(Q);

        const x1 = negative2RootQ * Math.cos(theta / 3.0) - adiv3;
        const x2 = negative2RootQ * Math.cos((theta + 2.0 * Math.PI) / 3.0) - adiv3;
        const x3 = negative2RootQ * Math.cos((theta - 2.0 * Math.PI) / 3.0) - adiv3;

        var roots = acceptRoot(x1);
        roots = roots.concat(acceptRoot(x2));
        roots = roots.concat(acceptRoot(x3));

        roots.sort();

        return deduplicateFloatArray(roots);
    }

    var A = Math.abs(R) + Math.sqrt(R2subQ3);

    A = Math.pow(A, 1.0 / 3.0);

    if (R > 0.0) {
        A = -A;
    }

    if (A != 0.0) {
        A += Q / A;
    }

    return acceptRoot(A - adiv3);
}


/**
 * Value returned when trying to decide orientation of 3 points.
 */
const TrianglePointOrientation = {

    /**
     * Points are in clockwise orientation.
     */
    Clockwise: 0,


    /**
     * Points are in counter-clockwise orientation.
     */
    CounterClockwise: 1,


    /**
     * Points are collinear and orientation cannot be determined.
     */
    Collinear: 2
}


/**
 * Represents a point in a two-dimensional plane defined as two floating point
 * values X and Y.
 */
class FloatPoint {
    X = 0;
    Y = 0;


    /**
     * Constructs point with given X and Y values.
     */
    constructor(x, y) {
        this.X = x;
        this.Y = y;
    }


    /**
     * Returns greater than zero if three given points are clockwise. Returns
     * less than zero if points are counter-clockwise. And returns zero if
     * points are collinear.
     */
    static turn(p1, p2, p3) {
        return p2.minus(p1).cross(p3.minus(p1));
    }


    /**
     * Determines orientation of triangle defined by three given points.
     */
    static determineTriangleOrientation(p1, p2, p3) {
        const t = FloatPoint.turn(p1, p2, p3);

        if (fuzzyIsZero(t)) {
            return TrianglePointOrientation.Collinear;
        } else if (t > 0.0) {
            return TrianglePointOrientation.Clockwise;
        }

        return TrianglePointOrientation.CounterClockwise;
    }


    /**
     * Returns true if triangle defined by three given points is clockwise.
     * Returns false if triangle is counter-clockwise or if points are
     * collinear.
     */
    static isTriangleClockwise(p1, p2, p3) {
        return FloatPoint.determineTriangleOrientation(p1, p2, p3) == TrianglePointOrientation.Clockwise;
    }


    /**
     * Returns true if this point is equal to another point.
     *
     * @param point Another point to compare with.
     */
    equals(point) {
        return fuzzyIsEqual(point.X, this.X) && fuzzyIsEqual(point.Y, this.Y);
    }


    /**
     * Returns true if this point is equal to another point.
     *
     * @param point Another point to compare with.
     *
     * @param epsilon Comparison epsilon. This value will be used when
     * comparing both X and Y components separately. Must be at least zero.
     */
    equalsWithEpsilon(point, epsilon) {
        return isEqualWithEpsilon(point.X, this.X, epsilon) &&
          isEqualWithEpsilon(point.Y, this.Y, epsilon);
    }


    /**
     * Returns distance from this point to another point.
     */
    distanceTo(point) {
        return point.minus(this).length();
    }


    /**
     * Returns squared distance from this point to another point.
     */
    distanceToSquared(point) {
        return point.minus(this).lengthSquared();
    }


    /**
     * Returns length of vector defined by X and Y components of this point.
     */
    length() {
        return Math.sqrt(this.lengthSquared());
    }


    /**
     * Returns squared length of vector defined by X and Y components of this
     * point.
     */
    lengthSquared() {
        return this.X * this.X + this.Y * this.Y;
    }


    /**
     * Returns normalized version of this vector. If this vector has length of
     * zero, vector with both components set to zero will be returned.
     */
    unitVector() {
        const mag2 = this.lengthSquared();

        if (mag2 != 0.0 && mag2 != 1.0) {
            const length = Math.sqrt(mag2);

            return new FloatPoint(this.X / length, this.Y / length);
        }

        return this;
    }


    /**
     * Returns vector which has direction perpendicular to the direction of
     * this vector.
     */
    normalVector() {
        return new FloatPoint(this.Y, -this.X);
    }


    /**
     * Returns vector which has direction perpendicular to the direction of
     * this vector. Before returning, resulting vector is normalized.
     */
    unitNormalVector() {
        return this.normalVector().unitVector();
    }


    /**
     * Returns cross product of two 2D vectors (this × point).
     *
     * Since both 2D vectors lie on the same XY plane, the only meaningful
     * return value is Z component of cross product. This method returns that
     * and does not calculate anything else.
     */
    cross(point) {
        return (this.X * point.Y) - (this.Y * point.X);
    }


    /**
     * Returns dot product of this vector and a given vector.
     */
    dot(point) {
        return (this.X * point.X) + (this.Y * point.Y);
    }


    /**
     * Rotates this vector by 90 degrees counter-clockwise and return rotated
     * version of this vector.
     */
    rotated90CCW() {
        return new FloatPoint(this.Y, -this.X);
    }


    /**
     * Returns point which is a result of linear interpolation between this
     * point and at given point.
     */
    lerp(to, t) {
        //ASSERT(t >= 0);
        //ASSERT(t <= 1);
        return this.plus(to.minus(this).multiplyScalar(t));
    }


    /**
     * Binary addition operator. Adds a given point to this point.
     *
     * Returns result as new point.
     */
    plus(point) {
        return new FloatPoint(this.X + point.X, this.Y + point.Y);
    }


    /**
     * Binary subtraction operator. Subtracts given point from this point.
     *
     * Returns result as new point.
     */
    minus(point) {
        return new FloatPoint(this.X - point.X, this.Y - point.Y);
    }


    /**
     * Binary multiplication operator. Multiplies given point by this point.
     *
     * Returns result as new point.
     */
    multiply(point) {
        return new FloatPoint(this.X * point.X, this.Y * point.Y);
    }


    /**
     * Binary scalar multiplication operator. Multiplies this point by a given
     * scalar.
     *
     * Returns result as new point.
     */
    multiplyScalar(v) {
        return new FloatPoint(this.X * v, this.Y * v);
    }


    /**
     * Binary division operator. Divides this point by a given point.
     *
     * Returns result as new point.
     */
    divide(point) {
        return new FloatPoint(this.X / point.X, this.Y / point.Y);
    }


    /**
     * Binary scalar division operator. Divides this point by a given scalar
     * value.
     *
     * Returns result as new point.
     */
    divideScalar(v) {
        return new FloatPoint(this.X / v, this.Y / v);
    }
}


/**
 * Describes type of intersection between two line segments.
 */
const LineIntersectionKind = {

    /**
     * No intersections found. Line segments are either zero length or they
     * are collinear.
     */
    None: 0,


    /**
     * Intersection was found within line segments.
     */
    Bounded: 1,


    /**
     * Intersection was found, but beyong line segments.
     */
    Unbounded: 2
}


/**
 * Describes intersection between two line segments.
 */
class LineIntersection {

    constructor(kind, intersectionPoint) {
        this.Kind = kind;
        this.IntersectionPoint = intersectionPoint;
    }


    /**
     * What kind of intersection was found.
     */
    Kind = LineIntersectionKind.None;


    /**
     * Intersection point. Only valid if intersection kind is something other
     * than None.
     */
    IntersectionPoint = new FloatPoint(0, 0);
}


/**
 * Represents a line segment defined by two points with floating-point
 * components.
 */
class FloatLine {
    P1 = new FloatPoint(0, 0);
    P2 = new FloatPoint(0, 0);


    /**
     * Constructs a line segment between two given points.
     */
    constructor(p1, p2) {
        this.P1 = p1;
        this.P2 = p2;
    }


    /**
     * Return X coordinate of the first point.
     */
    x1() {
        return this.P1.X;
    }


    /**
     * Return Y coordinate of the first point.
     */
    y1() {
        return this.P1.Y;
    }


    /**
     * Return X coordinate of the second point.
     */
    x2() {
        return this.P2.X;
    }


    /**
     * Return Y coordinate of the second point.
     */
    y2() {
        return this.P2.Y;
    }


    /**
     * Returns difference between X component of point 2 and point 1 of this
     * line.
     */
    dx() {
        return this.P2.X - this.P1.X;
    }


    /**
     * Returns difference between Y component of point 2 and point 1 of this
     * line.
     */
    dy() {
        return this.P2.Y - this.P1.Y;
    }


    /**
     * Returns true if both points of this line are at the same place within
     * given error.
     */
    isPoint() {
        return (fuzzyIsEqual(this.P1.X, this.P2.X) &&
            fuzzyIsEqual(this.P1.Y, this.P2.Y));
    }


    /**
     * Returns true if both points of this line are at the same place within
     * given error.
     */
    isPointWithEpsilon() {
        return (isEqualWithEpsilon(this.P1.X, this.P2.X, epsilon) &&
            isEqualWithEpsilon(this.P1.Y, this.P2.Y, epsilon));
    }


    /**
     * Returns a positive number of degrees from this linne to a given line.
     * Returned value is always positive and never exceeds 180.
     */
    getDegreesToLine(l) {
        return rad2Deg(this.getRadiansToLine(l));
    }


    /**
     * Returns a positive number of radians from this linne to a given line.
     * Returned value is always positive and never exceeds M_PI.
     */
    getRadiansToLine(l) {
        if (this.isPoint() || l.isPoint()) {
            return 0;
        }

        const c = (this.dx() * l.dx() + this.dy() * l.dy()) / (this.length() * l.length());

        const epsilon = Number.EPSILON * 8.0;

        // Return 0 instead of PI if c is outside range.
        if (c >= (-1.0 - epsilon) && c <= (1.0 + epsilon)) {
            return Math.acos(clamp(c, -1.0, 1.0));
        }

        return 0;
    }


    /**
     * Returns angle of this line in degrees. Angle is measured
     * counter-clock-wise from x axis.
     */
    angle() {
        const dx = this.P2.X - this.P1.X;
        const dy = this.P2.Y - this.P1.Y;
        const theta = rad2Deg(Math.atan2(-dy, dx));
        const thetaNormalized = theta < 0 ? theta + 360 : theta;

        if (fuzzyIsEqual(thetaNormalized, 360.0)) {
            return 0;
        }

        if (fuzzyIsZero(thetaNormalized)) {
            // In case we have -0, return positive zero.
            return 0;
        }

        return thetaNormalized;
    }


    /**
     * Returns length of this line.
     */
    length() {
        const x = this.P2.X - this.P1.X;
        const y = this.P2.Y - this.P1.Y;

        return Math.sqrt(x * x + y * y);
    }


    /**
     * Returns squared length of this line.
     */
    lengthSquared() {
        const x = this.P2.X - this.P1.X;
        const y = this.P2.Y - this.P1.Y;

        return x * x + y * y;
    }


    /**
     * Returns FloatLine with points reversed.
     */
    reversed() {
        return new FloatLine(this.P2, this.P1);
    }


    /**
     * Returns normalized vector which points to the same direction as this
     * line segment.
     */
    unitVector() {
        return this.P2.minus(this.P1).unitVector();
    }


    /**
     * Returns normal vector of this line. If both points of this line are in
     * the same place, returned vector will have both components set to zero.
     */
    normalVector() {
        return new FloatPoint(this.dy(), -this.dx());
    }


    /**
     * Returns normal vector of this line. Before returning, result is
     * normalized and will have length of 1. If both points of this line are
     * in the same place, returned vector will have both components set to
     * zero and its length will be zero.
     */
    unitNormalVector() {
        return this.normalVector().unitVector();
    }


    /**
     * Returns line intersection between this and other line segment.
     */
    intersect(l) {
        const a = this.P2.minus(this.P1);
        const b = l.P1.minus(l.P2);
        const denominator = a.Y * b.X - a.X * b.Y;

        if (denominator == 0) {
            return new LineIntersection();
        }

        const c = this.P1.minus(l.P1);
        const reciprocal = 1.0 / denominator;
        const na = (b.Y * c.X - b.X * c.Y) * reciprocal;

        const point = this.P1.plus(a.multiplyScalar(na));

        if (na < 0 || na > 1) {
            return new LineIntersection(LineIntersectionKind.Unbounded, point);
        }

        const nb = (a.X * c.Y - a.Y * c.X) * reciprocal;

        if (nb < 0 || nb > 1) {
            return new LineIntersection(LineIntersectionKind.Unbounded, point);
        }

        return LineIntersection(LineIntersectionKind.Bounded, point);
    }


    /**
     * Returns simple line intersection between this and other line. Use this
     * method if you are not interested in finding out if intersection is
     * within line segment or not.
     */
    intersectSimple(l) {
        const a = this.P2.minus(this.P1);
        const b = l.P1.minus(l.P2);
        const denominator = a.Y * b.X - a.X * b.Y;

        if (denominator == 0) {
            return null;
        }

        const c = this.P1.minus(l.P1);
        const reciprocal = 1.0 / denominator;
        const na = (b.Y * c.X - b.X * c.Y) * reciprocal;

        return this.P1.plus(a.multiplyScalar(na));
    }


    /**
     * Returns a copy of this line with both points translated by a given
     * amount.
     */
    translated(p) {
        return new FloatLine(this.P1.plus(p), this.P2.plus(p));
    }


    /**
     * Returns point which is in the middle of this line segment.
     */
    midPoint() {
        return new FloatPoint((this.P1.X + this.P2.X) * 0.5, (this.P1.Y + this.P2.Y) * 0.5);
    }


    /**
     * Extend this line segment by a given length. The first point of the line
     * will remain the same while the second point will be moved to the
     * direction of this line by a given distance. If points of this line
     * segment are equal or length value is zero, this method does nothing.
     *
     * @param length The distance by which to push the second point of this
     * line segment. Negative values are allowed.
     */
    extendByLengthFront(length) {
        if (this.isPoint() || fuzzyIsZero(length)) {
            return;
        }

        const v = this.unitVector();

        this.P2 = new FloatPoint(this.P2.X + (v.X * length),
            this.P2.Y + (v.Y * length));
    }


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
    extendByLengthBack(length) {
        if (this.isPoint() || fuzzyIsZero(length)) {
            return;
        }

        const v = this.unitVector();

        this.P1 = new FloatPoint(this.P1.X - (v.X * length),
            this.P1.Y - (v.Y * length));
    }


    /**
     * Returns true if a given point lies on this line segment.
     *
     * @param point Point to test.
     *
     * @param epsilon Tolerance for number comparison.
     */
    isPointOnLineSegmentWithEpsilon(point, epsilon) {
        const cross = (point.Y - this.P1.Y) * (this.P2.X - this.P1.X) -
            (this.point.X - this.P1.X) * (this.P2.Y - this.P1.Y);

        if (Math.abs(cross) > epsilon) {
            return false;
        }

        const dot = (point.X - this.P1.X) * (this.P2.X - this.P1.X) +
            (point.Y - this.P1.Y) * (this.P2.Y - this.P1.Y);

        if (dot < 0) {
            return false;
        }

        const sql = (this.P2.X - this.P1.X) * (this.P2.X - this.P1.X) +
            (this.P2.Y - this.P1.Y) * (this.P2.Y - this.P1.Y);

        return dot <= sql;
    }


    /**
     * Returns true if a given point lies on this line segment.
     *
     * @param point Point to test.
     */
    isPointOnLineSegment(point) {
        return this.isPointOnLineSegmentWithEpsilon(point, Number.EPSILON);
    }


    /**
     * Returns true if a given point lies on this line.
     *
     * @param point Point to test.
     *
     * @param epsilon Tolerance for number comparison.
     */
    isPointOnLineWithEpsilon(point, epsilon) {
        const cross = (point.Y - this.P1.Y) * (this.P2.X - this.P1.X) -
            (point.X - this.P1.X) * (this.P2.Y - this.P1.Y);

        return Math.abs(cross) <= epsilon;
    }


    /**
     * Returns true if a given point lies on this line.
     *
     * @param point Point to test.
     */
    isPointOnLine(point) {
        return this.isPointOnLineWithEpsilon(point, Number.EPSILON);
    }
}


/**
 * Represents a two dimensional cubic curve. It is defined by start point, end
 * point and two control points.
 */
class CubicCurve {

    /**
     * Start point of the curve.
     */
    P1 = new FloatPoint(0, 0);


    /**
     * First control point of the curve.
     */
    P2 = new FloatPoint(0, 0);


    /**
     * Second control point of the curve.
     */
    P3 = new FloatPoint(0, 0);


    /**
     * End point.
     */
    P4 = new FloatPoint(0, 0);


    /**
     * Constructs cubic curve with given points.
     *
     * @param p1 Start point of curve.
     *
     * @param p2 First control point of curve.
     *
     * @param p3 Second control point of curve.
     *
     * @param p4 End point of curve.
     */
    constructor(p1, p2, p3, p4) {
        this.P1 = p1;
        this.P2 = p2;
        this.P3 = p3;
        this.P4 = p4;
    }


    /**
     * Constructs cubic curve from quadratic curve parameters.
     */
    static createFromQuadratic(p1, p2, p3) {
        return new CubicCurve(p1,
            new FloatPoint(p1.X + (2.0 / 3.0) * (p2.X - p1.X),
                p1.Y + (2.0 / 3.0) * (p2.Y - p1.Y)),
            new FloatPoint(p2.X + (1.0 / 3.0) * (p3.X - p2.X),
                p2.Y + (1.0 / 3.0) * (p3.Y - p2.Y)), p4);
    }


    /**
     * Constructs cubic curve from a line. First control point of constructed
     * curve will be on line from p1 and p2 at 1/3rd of distance. Second
     * control point will be at 2/3rds of distance on line from p1 to p2.
     */
    static createFromLine(p1, p2) {
        return new CubicCurve(p1,
            new FloatPoint(p1.X + (1.0 / 3.0) * (p2.X - p1.X),
                p1.Y + (1.0 / 3.0) * (p2.Y - p1.Y)),
            new FloatPoint(p1.X + (2.0 / 3.0) * (p2.X - p1.X),
                p1.Y + (2.0 / 3.0) * (p2.Y - p1.Y)), p2);
    }


    /**
     * Returns true if all points of this curve are at the same place within
     * given tolerance.
     *
     * @param epsilon Tolerance value to test curve flatness against. Must
     * be greater than 0.0. Larger value means less accurate test.
     */
    isPointWithEpsilon(epsilon) {
        return this.P1.equalsWithEpsilon(this.P4, epsilon) &&
            this.P1.equalsWithEpsilon(this.P2, epsilon) &&
            this.P1.equalsWithEpsilon(this.P3, epsilon);
    }


    /**
     * Returns true if all points of this curve are at the same place within
     * given tolerance.
     */
    isPoint() {
        return this.isPointWithEpsilon(Number.EPSILON);
    }


    /**
     * Returns true if this curve is a straight line. Curve is straight line
     * when both control points lie on line segment between the first and the
     * last point of the curve.
     */
    isStraight() {
        const minx = Math.min(this.P1.X, this.P4.X);
        const miny = Math.min(this.P1.Y, this.P4.Y);
        const maxx = Math.max(this.P1.X, this.P4.X);
        const maxy = Math.max(this.P1.Y, this.P4.Y);

        // Is P2 located between P1 and P4?
        return minx <= this.P2.X &&
            miny <= this.P2.Y &&
            maxx >= this.P2.X &&
            maxy >= this.P2.Y &&
            // Is P3 located between P1 and P4?
            minx <= this.P3.X &&
            miny <= this.P3.Y &&
            maxx >= this.P3.X &&
            maxy >= this.P3.Y &&
            // Are all points collinear?
            fuzzyIsZero(FloatPoint.turn(this.P1, this.P2, this.P4)) &&
            fuzzyIsZero(FloatPoint.turn(this.P1, this.P3, this.P4));
    }


    /**
     * Returns start tangent of this curve. Start tangent is a line going from
     * P1 to P2, P3 or P4, depending which one is the first found to be not
     * equal to P1.
     *
     * @param epsilon Point comparison tolerance.
     */
    startTangentWithEpsilon(epsilon) {
        if (this.P1.equalsWithEpsilon(this.P2, epsilon)) {
            if (this.P1.equalsWithEpsilon(this.P3, epsilon)) {
                return new FloatLine(this.P1, this.P4);
            }

            return new FloatLine(this.P1, this.P3);
        }

        return new FloatLine(this.P1, this.P2);
    }


    /**
     * Returns start tangent of this curve. Start tangent is a line going from
     * P1 to P2, P3 or P4, depending which one is the first found to be not
     * equal to P1.
     */
    startTangent() {
        return this.startTangentWithEpsilon(Number.EPSILON);
    }


    /**
     * Returns end tangent of this curve. Start tangent is a line going from
     * P4 to P3, P2 or P1, depending which one is the first found to be not
     * equal to P4.
     *
     * @param epsilon Point comparison tolerance.
     */
    endTangentWithEpsilon(epsilon) {
        if (this.P4.equalsWithEpsilon(this.P3, epsilon)) {
            if (this.P4.equalsWithEpsilon(this.P2, epsilon)) {
                return new FloatLine(this.P4, this.P1);
            }

            return new FloatLine(this.P4, this.P2);
        }

        return new FloatLine(this.P4, this.P3);
    }


    /**
     * Returns end tangent of this curve. Start tangent is a line going from
     * P4 to P3, P2 or P1, depending which one is the first found to be not
     * equal to P4.
     */
    endTangent() {
        return this.endTangentWithEpsilon(Number.EPSILON);
    }


    /**
     * Returns point on curve at a given t.
     *
     * @param t Value from 0 to 1 to evaluate curve at.
     */
    pointAt(t) {
        const it = 1.0 - t;

        const a0 = this.P1.multiplyScalar(it).plus(this.P2.multiplyScalar(t));
        const b0 = this.P2.multiplyScalar(it).plus(this.P3.multiplyScalar(t));
        const c0 = this.P3.multiplyScalar(it).plus(this.P4.multiplyScalar(t));

        const a1 = a0.multiplyScalar(it).plus(b0.multiplyScalar(t));
        const b1 = b0.multiplyScalar(it).plus(c0.multiplyScalar(t));

        return a1.multiplyScalar(it).plus(b1.multiplyScalar(t));
    }


    /**
     * Returns normal vector of this curve at a given t. Returned vector is
     * not normalized.
     *
     * @param t Value from 0 to 1 to evaluate curve at.
     */
    normalVector(t) {
        if (fuzzyIsZero(t)) {
            if (this.P1.equals(this.P2)) {
                if (this.P1.equals(this.P3)) {
                    return new FloatLine(this.P1, this.P4).normalVector();
                } else {
                    return new FloatLine(this.P1, this.P3).normalVector();
                }
            }
        } else if (fuzzyIsEqual(t, 1.0)) {
            if (this.P3.equals(this.P4)) {
                if (this.P2.equals(this.P4)) {
                    return new FloatLine(this.P1, this.P4).normalVector();
                } else {
                    return new FloatLine(this.P2, this.P4).normalVector();
                }
            }
        }

        const d = this.derivedAt(t);

        return new FloatPoint(d.Y, -d.X);
    }


    /**
     * Returns normal vector of this curve at a given t. Length of returned
     * vector is 1.
     *
     * @param t Value from 0 to 1 to evaluate curve at.
     */
    unitNormalVector(t) {
        const n = this.normalVector(t);

        return n.unitVector();
    }


    /**
     * Returns first derivative at a given t.
     */
    derivedAt(t) {
        const it = 1.0 - t;
        const d = t * t;
        const a = -it * it;
        const b = 1.0 - 4.0 * t + 3.0 * d;
        const c = 2.0 * t - 3.0 * d;

        return new FloatPoint(a * this.P1.X + b * this.P2.X + c * this.P3.X + d * this.P4.X,
            a * this.P1.Y + b * this.P2.Y + c * this.P3.Y + d * this.P4.Y).multiplyScalar(3.0);
    }


    /**
     * Returns second derivative at a given t.
     */
    secondDerivedAt(t) {
        const a = 2.0 - 2.0 * t;
        const b = -4.0 + 6.0 * t;
        const c = 2.0 - 6.0 * t;
        const d = 2.0 * t;

        return new FloatPoint(a * this.P1.X + b * this.P2.X + c * this.P3.X + d * this.P4.X,
            a * this.P1.Y + b * this.P2.Y + c * this.P3.Y + d * this.P4.Y).multiplyScalar(3.0);
    }


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
    getSubcurve(t0, t1) {
        //ASSERT(t0 <= t1);

        if (fuzzyIsEqual(t0, t1)) {
            const p = this.pointAt(t0);

            return new CubicCurve(p, p, p, p);
        }

        if (t0 <= Number.EPSILON) {
            if (t1 >= (1.0 - Number.EPSILON)) {
                // Both parameters are 0.0 to 1.0.
                return this;
            }

            // Cut at t1 only.
            const ab = this.P1.lerp(this.P2, t1);
            const bc = this.P2.lerp(this.P3, t1);
            const cd = this.P3.lerp(this.P4, t1);
            const abc = ab.lerp(bc, t1);
            const bcd = bc.lerp(cd, t1);
            const abcd = abc.lerp(bcd, t1);

            return new CubicCurve(this.P1, ab, abc, abcd);
        }

        if (t1 >= (1.0 - Number.EPSILON)) {
            // Cut at t0 only.
            const ab = this.P1.lerp(this.P2, t0);
            const bc = this.P2.lerp(this.P3, t0);
            const cd = this.P3.lerp(this.P4, t0);
            const abc = ab.lerp(bc, t0);
            const bcd = bc.lerp(cd, t0);
            const abcd = abc.lerp(bcd, t0);

            return new CubicCurve(abcd, bcd, cd, this.P4);
        }

        // Cut at both t0 and t1.
        const ab0 = this.P1.lerp(this.P2, t1);
        const bc0 = this.P2.lerp(this.P3, t1);
        const cd0 = this.P3.lerp(this.P4, t1);
        const abc0 = ab0.lerp(bc0, t1);
        const bcd0 = bc0.lerp(cd0, t1);
        const abcd0 = abc0.lerp(bcd0, t1);

        const m = t0 / t1;

        const ab1 = this.P1.lerp(ab0, m);
        const bc1 = ab0.lerp(abc0, m);
        const cd1 = abc0.lerp(abcd0, m);
        const abc1 = ab1.lerp(bc1, m);
        const bcd1 = bc1.lerp(cd1, m);
        const abcd1 = abc1.lerp(bcd1, m);

        return new CubicCurve(abcd1, bcd1, cd1, abcd0);
    }


    /**
     * Finds points of maximum curvature and returns the number of points
     * found. Up to 3 points are returned.
     */
    findMaxCurvature() {
        const axx = this.P2.X - this.P1.X;
        const bxx = this.P3.X - 2.0 * this.P2.X + this.P1.X;
        const cxx = this.P4.X + 3.0 * (this.P2.X - this.P3.X) - this.P1.X;

        const cox0 = cxx * cxx;
        const cox1 = 3.0 * bxx * cxx;
        const cox2 = 2.0 * bxx * bxx + cxx * axx;
        const cox3 = axx * bxx;

        const ayy = this.P2.Y - this.P1.Y;
        const byy = this.P3.Y - 2.0 * this.P2.Y + this.P1.Y;
        const cyy = this.P4.Y + 3.0 * (this.P2.Y - this.P3.Y) - this.P1.Y;

        const coy0 = cyy * cyy;
        const coy1 = 3.0 * byy * cyy;
        const coy2 = 2.0 * byy * byy + cyy * ayy;
        const coy3 = ayy * byy;

        const coe0 = cox0 + coy0;
        const coe1 = cox1 + coy1;
        const coe2 = cox2 + coy2;
        const coe3 = cox3 + coy3;

        return findCubicRoots(coe0, coe1, coe2, coe3);
    }


    /**
     * Finds inflection points of this curve and returns a number of points
     * found.
     */
    findInflections() {
        const ax = this.P2.X - this.P1.X;
        const ay = this.P2.Y - this.P1.Y;
        const bx = this.P3.X - 2.0 * this.P2.X + this.P1.X;
        const by = this.P3.Y - 2.0 * this.P2.Y + this.P1.Y;
        const cx = this.P4.X + 3.0 * (this.P2.X - this.P3.X) - this.P1.X;
        const cy = this.P4.Y + 3.0 * (this.P2.Y - this.P3.Y) - this.P1.Y;

        return findQuadraticRoots(bx * cy - by * cx, ax * cy - ay * cx,
            ax * by - ay * bx);
    }


    /**
     * Split this curve in half and return array containing two curves.
     */
    split() {
        const c = this.P2.plus(this.P3).multiplyScalar(0.5);

        const aP2 = this.P1.plus(this.P2).multiplyScalar(0.5);
        const bP3 = this.P3.plus(this.P4).multiplyScalar(0.5);
        const aP3 = aP2.plus(c).multiplyScalar(0.5);
        const bP2 = bP3.plus(c).multiplyScalar(0.5);
        const m = aP3.plus(bP2).multiplyScalar(0.5);

        return [
            new CubicCurve(this.P1, aP2, aP3, m),
            new CubicCurve(m, bP2, bP3, this.P4)
        ];
    }


    /**
     * Finds intersections between this curve and line defined by two points.
     * Intersections found outside of line segment going from the first to the
     * second point will be included in the result. This method returns the
     * amount of intersections found. There can be up to 3 intersections.
     *
     * @param linePointA First point of line.
     *
     * @param linePointB Second point of line.
     */
    findRayIntersections(linePointA, linePointB) {
        const v = linePointB.minus(linePointA);

        const ax = (this.P1.Y - linePointA.Y) * v.X - (this.P1.X - linePointA.X) * v.Y;
        const bx = (this.P2.Y - linePointA.Y) * v.X - (this.P2.X - linePointA.X) * v.Y;
        const cx = (this.P3.Y - linePointA.Y) * v.X - (this.P3.X - linePointA.X) * v.Y;
        const dx = (this.P4.Y - linePointA.Y) * v.X - (this.P4.X - linePointA.X) * v.Y;

        const a = dx;
        const b = cx * 3;
        const c = bx * 3;

        const D = ax;
        const A = a - (D - c + b);
        const B = b + (3 * D - 2 * c);
        const C = c - (3 * D);

        return findCubicRoots(A, B, C, D);
    }
}


/**
 * Used for parallel curve construction.
 */
class CubicCurveBuilder {
    segments = [];

    /**
     * Adds line.
     */
    addLine(p1, p2) {
        this.segments.push(CubicCurve.createFromLine(p1, p2));
    }


    /**
     * Adds cubic curve.
     */
    addCubic(p1, cp1, cp2, to) {
        this.segments.push(new CubicCurve(p1, cp1, cp2, to));
    }


    /**
     * Returns the first point.
     */
    getFirstPoint() {
        return this.segments[0].P1;
    }


    /**
     * Returns the last point.
     */
    getLastPoint() {
        return this.segments[this.segments.length - 1].P4;
    }


    /**
     * Returns start tangent.
     */
    getStartTangent() {
        return this.segments[0].startTangent();
    }


    /**
     * Returns end tangent.
     */
    getEndTangent() {
        return this.segments[this.segments.length - 1].endTangent();
    }


    /**
     * Returns the number of curves.
     */
    getSegmentCount() {
        return this.segments.length;
    }


    /**
     * Returns output segment at a given index.
     *
     * @param index Segment index. Must be equal or greater than zero and less
     * than the value returned by getSegmentCount.
     */
    getSegmentAt(index) {
        return this.segments[index];
    }


    /**
     * Clears all segments in this builder.
     */
    reset() {
        this.segments.length = 0;
    }
};


/**
 * Keeps data needed to generate a set of output segments.
 */
class OutputBuilder {
    builder = null;
    previousPoint = new FloatPoint();
    previousPointT = new FloatPoint();
    cuspPoint = new FloatPoint();
    needsCuspArc = false;
    cuspArcClockwise = false;
    scale = 1;
    translation = new FloatPoint();

    constructor(builder, scale, translation) {
        this.builder = builder;
        this.scale = scale;
        this.translation = translation;
    }
};


/**
 * Called once when the first point of output is calculated.
 */
function moveTo(builder, to) {
    builder.previousPoint = to;
    builder.previousPointT = to.multiplyScalar(builder.scale).plus(builder.translation);
}


/**
 * Called when a new line needs to be added to the output. Line starts at the
 * last point of previously added segment or point set by a call to MoveTo.
 */
function lineTo(builder, to) {
    const previous = builder.previousPoint;

    if (!previous.equals(to)) {
        const t = to.multiplyScalar(builder.scale).plus(builder.translation);

        builder.builder.addLine(builder.previousPointT, t);

        builder.previousPoint = to;
        builder.previousPointT = t;
    }
}


/**
 * Called when a new cubic curve needs to be added to the output. Curve starts
 * at the last point of previously added segment or point set by a call to
 * MoveTo.
 */
function cubicTo(builder, cp1, cp2, to) {
    const previous = builder.previousPoint;

    if (!previous.equals(cp1) || !previous.equals(cp2) || !previous.equals(to)) {
        const c1 = cp1.multiplyScalar(builder.scale).plus(builder.translation);
        const c2 = cp2.multiplyScalar(builder.scale).plus(builder.translation);
        const t = to.multiplyScalar(builder.scale).plus(builder.translation);

        builder.builder.addCubic(builder.previousPointT, c1, c2, t);

        builder.previousPoint = to;
        builder.previousPointT = t;
    }
}


/**
 * Returns unit cubic curve for given circular arc parameters. Arc center is
 * assumed to be at 0, 0.
 *
 * @param p1 Starting point of circular arc. Both components must be in range
 * from -1 to 1.
 *
 * @param p4 End point of circular arc. Both components must be in range from
 * -1 to 1.
 */
function findUnitCubicCurveForArc(p1, p4) {
    const ax = p1.X;
    const ay = p1.Y;
    const bx = p4.X;
    const by = p4.Y;
    const q1 = ax * ax + ay * ay;
    const q2 = q1 + ax * bx + ay * by;
    const k2 = (4.0 / 3.0) * (Math.sqrt(2.0 * q1 * q2) - q2) /
        (ax * by - ay * bx);
    const x2 = p1.X - k2 * p1.Y;
    const y2 = p1.Y + k2 * p1.X;
    const x3 = p4.X + k2 * p4.Y;
    const y3 = p4.Y - k2 * p4.X;

    return new CubicCurve(p1, new FloatPoint(x2, y2), new FloatPoint(x3, y3), p4);
}


/**
 * Called when a circular arc needs to be added to the output. Arc starts at
 * the last point of previously added segment or point set by a call to MoveTo
 * and goes to a given end point.
 */
function arcTo(builder, center, to, clockwise) {
    const arcFrom = builder.previousPoint;

    const arcRadius = center.distanceTo(arcFrom);

    if (arcRadius < MinimumArcRadius) {
        return;
    }

    const centerToCurrentPoint = new FloatLine(center, arcFrom);
    const centerToEndPoint = new FloatLine(center, to);
    const startAngle = deg2Rad(centerToCurrentPoint.angle());

    var sweepAngle = centerToCurrentPoint.getRadiansToLine(centerToEndPoint);

    if (fuzzyIsZero(sweepAngle)) {
        return;
    }

    const determinedOrientation =
        FloatPoint.determineTriangleOrientation(center, arcFrom, to);

    if (determinedOrientation != TrianglePointOrientation.Collinear) {
        // If our three points are not collinear, we check if they are
        // clockwise. If we see that their orientation is opposite of what we
        // are told to draw, we draw large arc.
        const determinedClockwise =
            determinedOrientation == TrianglePointOrientation.Clockwise;

        if (determinedClockwise != clockwise) {
            sweepAngle = (Math.PI * 2.0) - sweepAngle;
        }
    }

    const nSteps = Math.ceil(sweepAngle / (Math.PI / 2.0));
    const step = sweepAngle / nSteps * (clockwise ? -1.0 : 1.0);

    var s = -Math.sin(startAngle);
    var c = Math.cos(startAngle);

    for (let i = 1; i <= nSteps; i++) {
        const a1 = startAngle + step * i;

        const s1 = -Math.sin(a1);
        const c1 = Math.cos(a1);

        const unitCurve = findUnitCubicCurveForArc(
            new FloatPoint(c, s), new FloatPoint(c1, s1));

        const p2 = unitCurve.P2.multiplyScalar(arcRadius).plus(center);
        const p3 = unitCurve.P3.multiplyScalar(arcRadius).plus(center);

        if (i < nSteps) {
            const p4 = unitCurve.P4.multiplyScalar(arcRadius).plus(center);

            cubicTo(builder, p2, p3, p4);
        } else {
            // Last point. Make sure we end with it. This is quite important
            // thing to do.
            cubicTo(builder, p2, p3, to);
        }

        s = s1;
        c = c1;
    }
}


function maybeAddCuspArc(builder, toPoint) {
    if (builder.needsCuspArc) {
        builder.needsCuspArc = false;

        arcTo(builder, builder.cuspPoint, toPoint, builder.cuspArcClockwise);

        builder.cuspPoint = new FloatPoint();
        builder.cuspArcClockwise = false;
    }
}


/**
 * Returns true if the curve is close enough to be considered parallel to the
 * original curve.
 *
 * @param original The original curve.
 *
 * @param parallel Candidate parallel curve to be tested.
 *
 * @param offset Offset from original curve to candidate parallel curve.
 *
 * @param maximumError Maximum allowed error.
 */
function acceptOffset(original, parallel, offset, maximumError) {
    // Using shape control method, sometimes output curve becomes completely
    // off in some situations involving start and end tangents being almost
    // parallel. These two checks are to prevent accepting such curves as good.
    if (FloatPoint.isTriangleClockwise(original.P1, original.P2, original.P4) !=
        FloatPoint.isTriangleClockwise(parallel.P1, parallel.P2, parallel.P4)) {
        return false;
    }

    if (FloatPoint.isTriangleClockwise(original.P1, original.P3, original.P4) !=
        FloatPoint.isTriangleClockwise(parallel.P1, parallel.P3, parallel.P4)) {
        return false;
    }

    for (let i = 0; i < SimpleOffsetProbePositions.length; i++) {
        const t = SimpleOffsetProbePositions[i];
        const p0 = original.pointAt(t);
        const n = original.normalVector(t);

        const roots = parallel.findRayIntersections(p0, p0.plus(n));

        if (roots.length != 1) {
            return false;
        }

        const p1 = parallel.pointAt(roots[0]);
        const d = p0.distanceTo(p1);
        const error = Math.abs(d - Math.abs(offset));

        if (error > maximumError) {
            return false;
        }
    }

    return true;
}


function arcOffset(builder, offset, center, from, to, clockwise) {
    var l1 = new FloatLine(center, from);
    var l2 = new FloatLine(center, to);

    if (clockwise) {
        l1.extendByLengthFront(offset);
        l2.extendByLengthFront(offset);
    } else {
        l1.extendByLengthFront(-offset);
        l2.extendByLengthFront(-offset);
    }

    maybeAddCuspArc(builder, l1.P2);

    // Determine if it is clockwise again since arc orientation may have
    // changed if arc radius was smaller than offset.
    //
    // Also it is important to use previous point to determine orientation
    // instead of the point we just calculated as the start of circular arc
    // because for small arcs a small numeric error can result in incorrect
    // arc orientation.
    arcTo(builder, center, l2.P2, FloatPoint.isTriangleClockwise(center,
        builder.previousPoint, l2.P2));
}


function unitTurn(point1, point2, point3) {
    return point2.minus(point1).unitVector().cross(point3.minus(point1).unitVector());
}


/**
 * Represents curve tangents as two line segments and some precomputed data.
 */
class CurveTangentData {
    startTangent = null;
    endTangent = null;
    turn1 = 0;
    turn2 = 0;
    startUnitNormal = null;
    endUnitNormal = null;

    constructor(curve) {
        this.startTangent = curve.startTangent();
        this.endTangent = curve.endTangent();
        this.turn1 = unitTurn(this.startTangent.P1, this.startTangent.P2,
            this.endTangent.P1);
        this.turn2 = unitTurn(this.startTangent.P1, this.endTangent.P2,
            this.endTangent.P1);
        this.startUnitNormal = this.startTangent.unitNormalVector();
        this.endUnitNormal = this.endTangent.unitNormalVector();
    }
};


/**
 * Returns true if an attempt to approximate a curve with given tangents
 * should be made.
 */
function canTryArcOffset(d) {
    // Arc approximation is only attempted if curve is not considered
    // approximately straight. But it can be attemped for curves which have
    // their control points on the different sides of line connecting points
    // P1 and P4.
    //
    // We need to make sure we don't try to do arc approximation for these S
    // type curves because such curves cannot be approximated by arcs in such
    // cases.

    const P = ApproximatelyStraightCurveTestApsilon;
    const N = -P;

    return (d.turn1 >= P && d.turn2 >= P) || (d.turn1 <= N && d.turn2 <= N);
}


/**
 * Returns true if an attempt to use simple offsetting for a curve with given
 * tangents should be made.
 */
function canTrySimpleOffset(d) {
    // Arc approximation is only attempted if curve is not considered
    // approximately straight. But it can be attemped for curves which have
    // their control points on the different sides of line connecting points
    // P1 and P4.
    //
    // We need to make sure we don't try to do arc approximation for these S
    // type curves because the shape control method behaves really badly with
    // S shape curves.

    return (d.turn1 >= 0 && d.turn2 >= 0) || (d.turn1 <= 0 && d.turn2 <= 0);
}


/**
 * Returns true if curve is considered too small to be added to offset output.
 */
function curveIsTooTiny(curve) {
    const lengthsSquared =
        curve.P1.distanceToSquared(curve.P2) +
        curve.P2.distanceToSquared(curve.P3) +
        curve.P3.distanceToSquared(curve.P4);

    return lengthsSquared <= MaximumTinyCurvePolygonPerimeterSquared;
}


/**
 * Attempts to perform simple curve offsetting and returns true if it succeeds
 * to generate good enough parallel curve.
 */
function trySimpleCurveOffset(curve, d, builder, offset, maximumError) {
    if (!canTrySimpleOffset(d)) {
        return false;
    }

    const d1 = curve.P2.minus(curve.P1);
    const d2 = curve.P3.minus(curve.P4);
    const div = d1.cross(d2);

    if (fuzzyIsZero(div)) {
        return false;
    }

    // Start point.
    const p1 = d.startTangent.P1.plus(d.startTangent.unitNormalVector().multiplyScalar(offset));

    // End point.
    const p4 = d.endTangent.P1.minus(d.endTangent.unitNormalVector().multiplyScalar(offset));

    // Middle point.
    const mp = curve.pointAt(0.5);
    const mpN = curve.unitNormalVector(0.5);
    const p = mp.plus(mpN.multiplyScalar(offset));

    const bxby = p.minus(p1.plus(p4).multiplyScalar(0.5)).multiplyScalar(8.0 / 3.0);

    const a = bxby.cross(d2) / div;
    const b = d1.cross(bxby) / div;

    const p2 = new FloatPoint(p1.X + a * d1.X, p1.Y + a * d1.Y);
    const p3 = new FloatPoint(p4.X + b * d2.X, p4.Y + b * d2.Y);

    const candidate = new CubicCurve(p1, p2, p3, p4);

    if (curveIsTooTiny(candidate)) {
        // If curve is too tiny, tell caller there was a great success.
        return true;
    }

    if (!acceptOffset(curve, candidate, offset, maximumError)) {
        return false;
    }

    maybeAddCuspArc(builder, candidate.P1);

    cubicTo(builder, candidate.P2, candidate.P3, candidate.P4);

    return true;
}


function arrayContainsCurvePosition(array, value, epsilon) {
    for (let i = 0; i < array.length; i++) {
        const v = array[i];

        if (isEqualWithEpsilon(value, v, epsilon)) {
            return true;
        }
    }

    return false;
}


function mergeCurvePositions(array, na, epsilon) {
    const va = array.concat(na);

    var gx = []

    for (let i = 0; i < va.length; i++) {
        const v = va[i];

        if (isZeroWithEpsilon(v, epsilon)) {
            continue;
        }

        if (isEqualWithEpsilon(v, 1.0, epsilon)) {
            continue;
        }

        if (arrayContainsCurvePosition(gx, v)) {
            continue;
        }

        gx.push(v);
    }

    return gx;
}


/**
 * Returns true if a given line segment intersects with circle. Only
 * intersection within line segment is considered.
 *
 * @param line Line segment.
 *
 * @param circleCenter Position of the circle center.
 *
 * @param circleRadius Circle radius. Must not be negative.
 */
function lineCircleIntersect(line, circleCenter, circleRadius) {
    //ASSERT(circleRadius >= 0);

    const d = line.P2.minus(line.P1);
    const g = line.P1.minus(circleCenter);
    const a = d.dot(d);
    const b = 2.0 * g.dot(d);
    const crSquared = circleRadius * circleRadius;
    const c = g.dot(g) - crSquared;
    const discriminant = b * b - 4.0 * a * c;

    if (discriminant > 0) {
        const dsq = Math.sqrt(discriminant);
        const a2 = a * 2.0;
        const t1 = (-b - dsq) / a2;
        const t2 = (-b + dsq) / a2;

        return (t1 >= 0.0 && t1 <= 1.0) || (t2 >= 0.0 && t2 <= 1.0);
    }

    return false;
}


/**
 * Returns true if circular arc with given parameters approximate curve close
 * enough.
 *
 * @param arcCenter Point where arc center is located.
 *
 * @param arcRadius Radius of arc.
 *
 * @param curve Curve being approximated with arc.
 *
 * @param maximumError Maximum allowed error.
 */
function goodArc(arcCenter, arcRadius, curve, maximumError, tFrom, tTo) {
    if (arcRadius > MaximumArcRadius) {
        return false;
    }

    const e = Math.min(maximumError, arcRadius / 3.0);

    // Calculate value equal to slightly more than half of maximum error.
    // Slightly more to minimize false negatives due to finite precision in
    // circle-line intersection test.
    const me = (e * (0.5 + 1e-4));

    for (let i = 0; i < ArcProbePositions.length; i++) {
        const t = ArcProbePositions[i];

        // Find t on a given curve.
        const curveT = interpolateLinear(t, tFrom, tTo);

        // Find point and normal at this position.
        const point = curve.pointAt(curveT);
        const n = curve.unitNormalVector(curveT);

        // Create line segment which has its center at curve on point and
        // extends by half of maximum allowed error to both directions from
        // curve point along normal.
        const segment = new FloatLine(point.plus(n.multiplyScalar(me)), point.minus(n.multiplyScalar(me)));

        // Test if intersection exists.
        if (!lineCircleIntersect(segment, arcCenter, arcRadius)) {
            return false;
        }
    }

    return true;
}


/**
 * Attempts to use circular arc offsetting method on a given curve.
 */
function tryArcApproximation(curve, d, builder, offset, maximumError) {
    if (!canTryArcOffset(d)) {
        return false;
    }

    // Cast ray from curve end points to start and end tangent directions.
    const vectorFrom = d.startTangent.unitVector();
    const vectorTo = d.endTangent.unitVector();
    const denom = vectorTo.X * vectorFrom.Y - vectorTo.Y * vectorFrom.X;

    // Should not happen as we already elliminated parallel case.
    if (fuzzyIsZero(denom)) {
        return false;
    }

    const asv = d.startTangent.P1;
    const bsv = d.endTangent.P1;
    const u = ((bsv.Y - asv.Y) * vectorTo.X - (bsv.X - asv.X) * vectorTo.Y) / denom;
    const v = ((bsv.Y - asv.Y) * vectorFrom.X - (bsv.X - asv.X) * vectorFrom.Y) / denom;

    if (u < 0.0 || v < 0.0) {
        // Intersection is on the wrong side.
        return false;
    }

    const V = asv.plus(vectorFrom.multiplyScalar(u));

    // If start or end tangents extend too far beyond intersection, return
    // early since it will not result in good approximation.
    if (curve.P1.distanceToSquared(V) < (d.startTangent.lengthSquared() * 0.25) ||
        curve.P4.distanceToSquared(V) < (d.endTangent.lengthSquared() * 0.25)) {
        return false;
    }

    const P2VDistance = curve.P4.distanceTo(V);
    const P1VDistance = curve.P1.distanceTo(V);
    const P1P4Distance = curve.P1.distanceTo(curve.P4);
    const G = curve.P1.multiplyScalar(P2VDistance).plus(curve.P4.multiplyScalar(P1VDistance).plus(V.multiplyScalar(P1P4Distance))).divideScalar(P2VDistance + P1VDistance + P1P4Distance);

    const P1G = new FloatLine(curve.P1, G);
    const GP4 = new FloatLine(G, curve.P4);

    const E = new FloatLine(P1G.midPoint(), P1G.midPoint().minus(P1G.normalVector()));
    const E1 = new FloatLine(d.startTangent.P1, d.startTangent.P1.minus(d.startTangent.normalVector()));
    const C1 = E.intersectSimple(E1);

    if (C1 == null) {
        return false;
    }

    const roots = curve.findRayIntersections(C1, G);

    if (roots.length != 1) {
        return false;
    }

    const tG = roots[0];
    const dist0 = G.distanceTo(curve.pointAt(tG));

    if (dist0 > maximumError) {
        return false;
    }

    const F = new FloatLine(GP4.midPoint(), GP4.midPoint().minus(GP4.normalVector()));
    const F1 = new FloatLine(d.endTangent.P1, d.endTangent.P1.plus(d.endTangent.normalVector()));
    const C2 = F.intersectSimple(F1);

    if (C2 == null) {
        return false;
    }

    if (C1.equalsWithEpsilon(C2, ArcCenterComparisonEpsilon)) {
        const radius = C1.distanceTo(curve.P1);

        if (goodArc(C1, radius, curve, maximumError, 0, 1)) {
            const clockwise = FloatPoint.isTriangleClockwise(curve.P1, V,
                curve.P4);

            arcOffset(builder, offset, C1, curve.P1, curve.P4, clockwise);

            return true;
        }
    } else {
        const radius1 = C1.distanceTo(curve.P1);

        if (!goodArc(C1, radius1, curve, maximumError, 0, tG)) {
            return false;
        }

        const radius2 = C2.distanceTo(curve.P4);

        if (!goodArc(C2, radius2, curve, maximumError, tG, 1)) {
            return false;
        }

        const clockwise = FloatPoint.isTriangleClockwise(curve.P1, V,
            curve.P4);

        arcOffset(builder, offset, C1, curve.P1, G, clockwise);
        arcOffset(builder, offset, C2, G, curve.P4, clockwise);

        return true;
    }

    return false;
}


function isCurveApproximatelyStraight(d) {
    const minx = Math.min(d.startTangent.x1(), d.endTangent.x1());
    const miny = Math.min(d.startTangent.y1(), d.endTangent.y1());
    const maxx = Math.max(d.startTangent.x1(), d.endTangent.x1());
    const maxy = Math.max(d.startTangent.y1(), d.endTangent.y1());

    const x2 = d.startTangent.x2();
    const y2 = d.startTangent.y2();
    const x3 = d.endTangent.x2();
    const y3 = d.endTangent.y2();

    // Is P2 located between P1 and P4?
    return minx <= x2 &&
        miny <= y2 &&
        maxx >= x2 &&
        maxy >= y2 &&
        // Is P3 located between P1 and P4?
        minx <= x3 &&
        miny <= y3 &&
        maxx >= x3 &&
        maxy >= y3 &&
        // Are all points collinear?
        isZeroWithEpsilon(d.turn1,
            ApproximatelyStraightCurveTestApsilon) &&
        isZeroWithEpsilon(d.turn2,
            ApproximatelyStraightCurveTestApsilon);
}


function curveIsCompletelyStraight(d) {
    return isZeroWithEpsilon(d.turn1, CompletelyStraightCurveTestApsilon) &&
        isZeroWithEpsilon(d.turn2, CompletelyStraightCurveTestApsilon);
}


/**
 * Main function for approximating offset of a curve without cusps.
 */
function approximateBezier(curve, d, builder, offset, maximumError) {
    if (!curve.isPointWithEpsilon(CurvePointClumpTestEpsilon)) {
        if (isCurveApproximatelyStraight(d)) {
            if (curveIsCompletelyStraight(d)) {
                // Curve is extremely close to being straight.
                const line = new FloatLine(curve.P1, curve.P2);
                const normal = line.unitNormalVector();

                maybeAddCuspArc(builder, line.P1.plus(normal.multiplyScalar(offset)));

                lineTo(builder, line.P2.plus(normal.multiplyScalar(offset)));
            } else {
                const p1o = d.startTangent.P1.plus(d.startUnitNormal.multiplyScalar(offset));
                const p2o = d.startTangent.P2.plus(d.startUnitNormal.multiplyScalar(offset));
                const p3o = d.endTangent.P2.minus(d.endUnitNormal.multiplyScalar(offset));
                const p4o = d.endTangent.P1.minus(d.endUnitNormal.multiplyScalar(offset));

                maybeAddCuspArc(builder, p1o);

                cubicTo(builder, p2o, p3o, p4o);
            }
        } else {
            if (!trySimpleCurveOffset(curve, d, builder, offset, maximumError)) {
                if (!tryArcApproximation(curve, d, builder, offset, maximumError)) {
                    // Split in half and continue.
                    const curves = curve.split();

                    for (let i = 0; i < curves.length; i++) {
                        const da = new CurveTangentData(curves[i]);

                        approximateBezier(curves[i], da, builder, offset, maximumError);
                    }
                }
            }
        }
    }
}


function findPositionOnCurveWithLargeEnoughDerivative(curve, previousT, currentT) {
    // ASSERT(currentT > previousT);

    const kPrecision = CuspDerivativeLengthSquared * 2.0;

    var t = Math.max(interpolateLinear(previousT, currentT, 0.8), currentT - 0.05);

    for (let i = 0; i < NearCuspPointSearchMaxIterationCount; i++) {
        const derivative = curve.derivedAt(t);
        const lengthSquared = derivative.lengthSquared();

        if (lengthSquared < kPrecision) {
            return t;
        }

        const a = t + currentT;

        t = a / 2.0;
    }

    return t;
}


function findPositionOnCurveWithLargeEnoughDerivativeStart(curve, currentT, nextT) {
    // ASSERT(currentT < nextT);

    const kPrecision = CuspDerivativeLengthSquared * 2.0;

    var t = Math.min(interpolateLinear(currentT, nextT, 0.2), currentT + 0.05);

    for (let i = 0; i < NearCuspPointSearchMaxIterationCount; i++) {
        const derivative = curve.derivedAt(t);
        const lengthSquared = derivative.lengthSquared();

        if (lengthSquared < kPrecision) {
            return t;
        }

        const a = currentT + t;

        t = a / 2.0;
    }

    return t;
}


/**
 * If all points of the curve are collinear, a shortcut must be made because
 * general offsetting algorithm does not handle such curves very well. In case
 * where are points are collinear, lines between cusps are offset to direction
 * of their normals and at the points where curve has a cusps, semi-circles
 * are added to the output.
 */
function offsetLinearCuspyCurve(curve, builder, offset, maximumCurvaturePoints) {
    const startTangent = curve.startTangent();
    const normal = startTangent.unitNormalVector();

    var previousPoint = startTangent.P1;
    var previousOffsetPoint = previousPoint.plus(normal.multiplyScalar(offset));

    moveTo(builder, previousOffsetPoint);

    for (let i = 0; i < maximumCurvaturePoints.length; i++) {
        // Skip 0 and 1!
        const t = maximumCurvaturePoints[i];
        const derived = curve.derivedAt(t);
        const lengthSquared = derived.lengthSquared();

        if (lengthSquared <= 1e-9) {
            // Cusp. Since we know all curve points are on the same line, some
            // of maximum curvature points will have nearly zero length
            // derivative vectors.
            const pointAtCusp = curve.pointAt(t);

            // Draw line from previous point to point at cusp.
            const l = new FloatLine(previousPoint, pointAtCusp);
            const n = l.unitNormalVector();
            const to = pointAtCusp.plus(n.multiplyScalar(offset));

            lineTo(builder, to);

            // Draw semi circle at cusp.
            const arcTo = pointAtCusp.minus(n.multiplyScalar(offset));

            arcTo(builder, pointAtCusp, arcTo, FloatPoint.isTriangleClockwise(
                previousPoint, previousOffsetPoint, pointAtCusp));

            previousPoint = pointAtCusp;
            previousOffsetPoint = arcTo;
        }
    }

    const endTangent = curve.endTangent();
    const normal2 = endTangent.unitNormalVector();

    lineTo(builder, endTangent.P1.minus(normal2.multiplyScalar(offset)));
}


function doApproximateBezier(curve, d, builder, offset, maximumError) {
    // First find maximum curvature positions.
    const maximumCurvaturePositions = curve.findMaxCurvature();

    // Handle special case where the input curve is a straight line, but
    // control points do not necessary lie on line segment between curve
    // points P1 and P4.
    if (curveIsCompletelyStraight(curve)) {
        offsetLinearCuspyCurve(curve, builder, offset,
            maximumCurvaturePositions);
        return;
    }

    // Now find inflection point positions.
    const inflections = curve.findInflections();

    // Merge maximum curvature and inflection point positions.
    const t = mergeCurvePositions(maximumCurvaturePositions, inflections);

    t.sort();

    if (t.length == 0) {
        // No initial subdivision suggestions.
        approximateBezier(curve, d, builder, offset, maximumError);
    } else {
        var previousT = 0;

        for (let i = 0; i < t.length; i++) {
            const T = t[i];
            const derivative = curve.derivedAt(T);
            const lengthSquared = derivative.lengthSquared();

            if (lengthSquared < CuspDerivativeLengthSquared) {
                // Squared length of derivative becomes tiny. This is where
                // the cusp is. The goal here is to find a spon on curve,
                // located before T which has large enough derivative and draw
                // circular arc to the next point on curve with large enough
                // derivative.

                const t1 = findPositionOnCurveWithLargeEnoughDerivative(
                    curve, previousT, T);

                // ASSERT(t1 < T);

                const k = curve.getSubcurve(previousT, t1);
                const nd = new CurveTangentData(k);

                approximateBezier(k, nd, builder, offset, maximumError);

                const t2 = findPositionOnCurveWithLargeEnoughDerivativeStart(
                    curve, T, i == (t.length - 1) ? 1.0 : t[i + 1]);

                // ASSERT(t2 > T);

                builder.cuspPoint = curve.pointAt(T);
                builder.needsCuspArc = true;
                builder.cuspArcClockwise = FloatPoint.isTriangleClockwise(
                    k.P4, builder.cuspPoint, curve.pointAt(t2));

                previousT = t2;
            } else {
                // Easy, feed subcurve between previous and current t values
                // to offset approximation function.

                const k = curve.getSubcurve(previousT, T);
                const nd = new CurveTangentData(k);

                approximateBezier(k, nd, builder, offset, maximumError);

                previousT = T;
            }
        }

        // ASSERT(previousT < 1.0);

        const k = curve.getSubcurve(previousT, 1.0);
        const nd = new CurveTangentData(k);

        approximateBezier(k, nd, builder, offset, maximumError);
    }
}


/**
 * Flattens ends of curves if control points are too close to end points.
 */
function fixRedundantTangents(curve) {
    var p2 = curve.P2;
    var p3 = curve.P3;

    if (curve.P1.distanceToSquared(p2) < 1e-12) {
        p2 = curve.P1;
    }

    if (curve.P4.distanceToSquared(p3) < 1e-12) {
        p3 = curve.P4;
    }

    return new CubicCurve(curve.P1, p2, p3, curve.P4);
}


function offsetCurve(curve, offset, maximumError, builder) {
    builder.reset();

    const minx = min4(curve.P1.X, curve.P2.X, curve.P3.X, curve.P4.X);
    const maxx = max4(curve.P1.X, curve.P2.X, curve.P3.X, curve.P4.X);
    const miny = min4(curve.P1.Y, curve.P2.Y, curve.P3.Y, curve.P4.Y);
    const maxy = max4(curve.P1.Y, curve.P2.Y, curve.P3.Y, curve.P4.Y);

    const dx = maxx - minx;
    const dy = maxy - miny;

    if (dx < CurvePointClumpTestEpsilon && dy < CurvePointClumpTestEpsilon) {
        return;
    }

    // Select bigger of width and height.
    const m = Math.max(dx, dy) / 2.0;

    // Calculate scaled offset.
    const so = offset / m;

    if (fuzzyIsZero(so)) {
        builder.addCubic(curve.P1, curve.P2, curve.P3, curve.P4);
        return;
    }

    // Calculate "normalized" curve which kind of fits into range from -1 to 1.
    const tx = (minx + maxx) / 2.0;
    const ty = (miny + maxy) / 2.0;
    const t = new FloatPoint(tx, ty);

    const p1 = curve.P1.minus(t);
    const p2 = curve.P2.minus(t);
    const p3 = curve.P3.minus(t);
    const p4 = curve.P4.minus(t);

    const sc = new CubicCurve(p1.divideScalar(m), p2.divideScalar(m),
        p3.divideScalar(m), p4.divideScalar(m));

    const c = fixRedundantTangents(sc);

    b = new OutputBuilder(builder, m, t);

    const d = new CurveTangentData(c);

    if (isCurveApproximatelyStraight(d)) {
        if (curveIsCompletelyStraight(d)) {
            // Curve is extremely close to being straight, use simple line
            // translation.
            const line = new FloatLine(c.P1, c.P4);
            const normal = line.unitNormalVector();
            const translated = line.translated(normal.multiplyScalar(so));

            moveTo(b, translated.P1);

            lineTo(b, translated.P2);
        } else {
            // Curve is almost straight. Translate start and end tangents
            // separately and generate a cubic curve.
            const p1o = d.startTangent.P1.plus(d.startUnitNormal.multiplyScalar(so));
            const p2o = d.startTangent.P2.plus(d.startUnitNormal.multiplyScalar(so));
            const p3o = d.endTangent.P2.minus(d.endUnitNormal.multiplyScalar(so));
            const p4o = d.endTangent.P1.minus(d.endUnitNormal.multiplyScalar(so));

            moveTo(b, p1o);

            cubicTo(b, p2o, p3o, p4o);
        }
    } else {
        // Arbitrary curve.
        moveTo(b, d.startTangent.P1.plus(d.startUnitNormal.multiplyScalar(so)));

        // Try arc approximation first in case this curve was intended to
        // approximate circle. If that is indeed true, we avoid a lot of
        // expensive calculations like finding inflection and maximum
        // curvature points.
        if (!tryArcApproximation(c, d, b, so, maximumError)) {
            doApproximateBezier(c, d, b, so, maximumError);
        }
    }
}
