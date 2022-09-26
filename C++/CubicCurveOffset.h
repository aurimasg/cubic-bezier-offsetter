
#pragma once


#include "CubicCurve.h"


/**
 * Used for parallel curve construction. Users need to override methods to
 * accept offsetter output.
 */
class CubicCurveBuilder {
public:
    virtual ~CubicCurveBuilder() {
    }
public:

    /**
     * Called exactly once for each offsetter invocation when called with
     * curve which is not a point.
     *
     * @param point Point to move to.
     */
    virtual void MoveTo(const FloatPoint &point) = 0;


    /**
     * Draw line from previous to a given point.
     *
     * @param point End point of a line.
     */
    virtual void LineTo(const FloatPoint &point) = 0;


    /**
     * Draw cubic curve from previous point to a given point.
     *
     * @param cp1 The first control point of a cubic curve.
     *
     * @param cp2 The second control point of a cubic curve.
     *
     * @param to End point of a cubic curve.
     */
    virtual void CubicTo(const FloatPoint &cp1, const FloatPoint &cp2,
        const FloatPoint &to) = 0;

public:
    FloatPoint PreviousPoint;
    FloatPoint CuspPoint;
    bool NeedsCuspArc = false;
    bool CuspArcClockwise = false;
};


/**
 * Find a set of segments that approximate parallel curve.
 *
 * @param curve Input curve.
 *
 * @param offset Offset amount. If it is zero, resulting curve will be
 * identical to input curve. Can be negative.
 *
 * @param maximumError Maximum error. Lower value means better precision and
 * more output segments. Larger value means worse precision, but fewer output
 * segments.
 *
 * @param builder Output receiver.
 */
extern void OffsetCurve(const CubicCurve &curve, const double offset,
    const double maximumError, CubicCurveBuilder &builder);
