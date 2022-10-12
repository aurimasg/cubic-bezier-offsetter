
#pragma once


#include "CubicCurve.h"
#include <new>
#include <cstring>


/**
 * Used for parallel curve construction.
 */
class CubicCurveBuilder final {
public:
    CubicCurveBuilder() {
    }
public:
   ~CubicCurveBuilder();
public:

    /**
     * Adds line.
     */
    void AddLine(const FloatPoint &p1, const FloatPoint &p2);


    /**
     * Adds cubic curve.
     */
    void AddCubic(const FloatPoint &p1, const FloatPoint &cp1,
        const FloatPoint &cp2, const FloatPoint &to);


    /**
     * Returns the first point.
     */
    FloatPoint GetFirstPoint() const;


    /**
     * Returns the last point.
     */
    FloatPoint GetLastPoint() const;


    /**
     * Returns start tangent.
     */
    FloatLine GetStartTangent() const;


    /**
     * Returns end tangent.
     */
    FloatLine GetEndTangent() const;


    /**
     * Returns the number of curves.
     */
    int GetSegmentCount() const;


    /**
     * Returns pointer to the output segment at a given index.
     *
     * @param index Segment index. Must be equal or greater than zero and less
     * than the value returned by GetSegmentCount.
     */
    const CubicCurve *GetSegmentAt(const int index) const;


    /**
     * Clears all segments in this builder.
     */
    void Reset();

private:
    void MakeRoomForCurve();
private:
    static constexpr int EmbeddedCurveCount = 32;

    char mEmbeddedMemory[SIZE_OF(CubicCurve) * EmbeddedCurveCount];
    CubicCurve *mCurves = reinterpret_cast<CubicCurve *>(mEmbeddedMemory);

    // How many curves there are.
    int mCurveCount = 0;

    // For how many curves memory was allocated.
    int mCurveCapacity = EmbeddedCurveCount;
private:
    DISABLE_COPY_AND_ASSIGN(CubicCurveBuilder);
};


FORCE_INLINE CubicCurveBuilder::~CubicCurveBuilder() {
    if (mCurves != reinterpret_cast<CubicCurve *>(mEmbeddedMemory)) {
        free(mCurves);
    }
}


FORCE_INLINE void CubicCurveBuilder::AddLine(const FloatPoint &p1, const FloatPoint &p2) {
    MakeRoomForCurve();
    void *p = mCurves + mCurveCount;
    new (p) CubicCurve(p1, p2);
    mCurveCount++;
}


FORCE_INLINE void CubicCurveBuilder::AddCubic(const FloatPoint &p1, const FloatPoint &cp1, const FloatPoint &cp2, const FloatPoint &to) {
    MakeRoomForCurve();
    void *p = mCurves + mCurveCount;
    new (p) CubicCurve(p1, cp1, cp2, to);
    mCurveCount++;
}


FORCE_INLINE FloatPoint CubicCurveBuilder::GetFirstPoint() const {
    ASSERT(mCurveCount > 0);

    return mCurves->P1;
}


FORCE_INLINE FloatPoint CubicCurveBuilder::GetLastPoint() const {
    ASSERT(mCurveCount > 0);

    return mCurves[mCurveCount - 1].P4;
}


FORCE_INLINE FloatLine CubicCurveBuilder::GetStartTangent() const {
    ASSERT(mCurveCount > 0);

    return mCurves->StartTangent();
}


FORCE_INLINE FloatLine CubicCurveBuilder::GetEndTangent() const {
    ASSERT(mCurveCount > 0);

    return mCurves[mCurveCount - 1].EndTangent();
}


FORCE_INLINE void CubicCurveBuilder::MakeRoomForCurve() {
    if (mCurveCapacity == mCurveCount) {
        const int newCapacity = mCurveCapacity * 2;

        CubicCurve *curves = reinterpret_cast<CubicCurve *>(
            malloc(SIZE_OF(CubicCurve) * newCapacity));

        memcpy(curves, mCurves, SIZE_OF(CubicCurve) * mCurveCount);

        if (mCurves != reinterpret_cast<CubicCurve *>(mEmbeddedMemory)) {
            free(mCurves);
        }

        mCurves = curves;
        mCurveCapacity = newCapacity;
    }
}


FORCE_INLINE int CubicCurveBuilder::GetSegmentCount() const {
    return mCurveCount;
}


FORCE_INLINE const CubicCurve *CubicCurveBuilder::GetSegmentAt(const int index) const {
    ASSERT(index >= 0);
    ASSERT(index < mCurveCount);

    return mCurves + index;
}


FORCE_INLINE void CubicCurveBuilder::Reset() {
    mCurveCount = 0;
}


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
