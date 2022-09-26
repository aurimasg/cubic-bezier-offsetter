
#pragma once


#include <cassert>
#include <cfloat>
#include <cmath>


#define FORCE_INLINE inline

#ifdef DEBUG
#define ASSERT(x) assert(x)
#else
#define ASSERT(x)
#endif

#define ARRAY_SIZE(a) \
    ((sizeof(a) / sizeof(*(a))) / \
      static_cast<size_t>(!(sizeof(a) % sizeof(*(a)))))


/**
 * Returns absolute of a given value.
 */
template <typename T>
constexpr T Abs(const T t) {
    return t >= 0 ? t : -t;
}


/**
 * Returns greater of two given values.
 */
template <typename T>
constexpr T Max(const T a, const T b) {
    return a < b ? b : a;
}


/**
 * Returns smaller of two given values.
 */
template <typename T>
constexpr T Min(const T a, const T b) {
    return a < b ? a : b;
}


/**
 * Returns value clamped to range between minimum and maximum values.
 */
template <typename T>
constexpr T Clamp(const T val, const T min, const T max) {
    return val > max ? max : val < min ? min : val;
}


/**
 * Convert degrees to radians.
 */
constexpr double Deg2Rad(const double x) {
    // pi / 180
    return x * 0.01745329251994329576923690768489;
}


/**
 * Convert radians to degrees.
 */
constexpr double Rad2Deg(const double x) {
    // 180 / pi.
    return x * 57.295779513082320876798154814105;
}


/**
 * Linearly interpolate between A and B.
 * If t is 0, returns A.
 * If t is 1, returns B.
 * If t is something else, returns value linearly interpolated between A and B.
 */
template <typename T, typename V>
static FORCE_INLINE T InterpolateLinear(const T A, const T B, const V t) {
    ASSERT(t >= 0);
    ASSERT(t <= 1);
    return A + ((B - A) * t);
}


/**
 * Returns true if two numbers differ not more than a given epsilon.
 */
static FORCE_INLINE bool IsEqualWithEpsilon(const double a, const double b, const double epsilon) {
    return (fabs(a - b) < epsilon);
}


/**
 * Returns true if two numbers differ not more than a given epsilon.
 */
static FORCE_INLINE bool IsEqualWithEpsilon(const float a, const float b, const float epsilon) {
    return (fabsf(a - b) < epsilon);
}


/**
 * Returns true if a number differs from zero by less than a given epsilon.
 */
static FORCE_INLINE bool IsZeroWithEpsilon(const double a, const double epsilon) {
    return (fabs(a) < epsilon);
}


/**
 * Returns true if a number differs from zero by less than a given epsilon.
 */
static FORCE_INLINE bool IsZeroWithEpsilon(const float a, const float epsilon) {
    return (fabsf(a) < epsilon);
}


/**
 * Returns true if two given numbers are considered equal.
 */
static FORCE_INLINE bool FuzzyIsEqual(const double a, const double b) {
    return (fabs(a - b) < DBL_EPSILON);
}


/**
 * Returns true if two given numbers are considered equal.
 */
static FORCE_INLINE bool FuzzyIsEqual(const float a, const float b) {
    return (fabsf(a - b) < FLT_EPSILON);
}


/**
 * Returns true if a number can be considered being equal to zero.
 */
static FORCE_INLINE bool FuzzyIsZero(const double d) {
    return fabs(d) < DBL_EPSILON;
}


/**
 * Returns true if a number can be considered being equal to zero.
 */
static FORCE_INLINE bool FuzzyIsZero(const float f) {
    return fabsf(f) < FLT_EPSILON;
}


/**
 * Returns square root of a given number.
 */
static FORCE_INLINE float Sqrt(const float v) {
    return sqrtf(v);
}


/**
 * Returns square root of a given number.
 */
static FORCE_INLINE double Sqrt(const double v) {
    return sqrt(v);
}


/**
 * Calculates sine of a given number.
 */
static FORCE_INLINE float Sin(const float v) {
    return sinf(v);
}


/**
 * Calculates sine of a given number.
 */
static FORCE_INLINE double Sin(const double v) {
    return sin(v);
}


/**
 * Calculate cosine of a given number.
 */
static FORCE_INLINE float Cos(const float v) {
    return cosf(v);
}


/**
 * Calculate cosine of a given number.
 */
static FORCE_INLINE double Cos(const double v) {
    return cos(v);
}


/**
 * Returns the angle whose cosine is the specified number. Returned angle is
 * in radians.
 *
 * @param v A number representing a cosine which must be greater than or
 * equal to -1 and less than or equal to 1.
 */
static FORCE_INLINE float Acos(const float v) {
    ASSERT(v <= 1.0f);
    ASSERT(v >= -1.0f);
    return acosf(v);
}


/**
 * Returns the angle whose cosine is the specified number. Returned angle is
 * in radians.
 *
 * @param v A number representing a cosine which must be greater than or
 * equal to -1 and less than or equal to 1.
 */
static FORCE_INLINE double Acos(const double v) {
    ASSERT(v <= 1.0);
    ASSERT(v >= -1.0);
    return acos(v);
}


/**
 * Returns tangent of a given number.
 */
static FORCE_INLINE float Tan(const float v) {
    return tanf(v);
}


/**
 * Returns tangent of a given number.
 */
static FORCE_INLINE double Tan(const double v) {
    return tan(v);
}


/**
 * Returns the angle expressed in radians between the positive X-axis and the
 * ray cast from 0, 0 to the point x, y.
 */
static FORCE_INLINE float Atan2(const float y, const float x) {
    return atan2f(y, x);
}


/**
 * Returns the angle expressed in radians between the positive X-axis and the
 * ray cast from 0, 0 to the point x, y.
 */
static FORCE_INLINE double Atan2(const double y, const double x) {
    return atan2(y, x);
}


/**
 * Returns true if a given floating point value is not a number.
 */
constexpr bool IsNaN(const float x) {
    return x != x;
}


/**
 * Returns true if a given double precision floating point value is not a
 * number.
 */
constexpr bool IsNaN(const double x) {
    return x != x;
}


/**
 * Returns true if a given double precision floating point number is finite.
 */
static FORCE_INLINE bool DoubleIsFinite(const double x) {
    // 0 × finite → 0
    // 0 × infinity → NaN
    // 0 × NaN → NaN
    const double p = x * 0;

    return !IsNaN(p);
}


/**
 * Rounds-up a given number.
 */
static FORCE_INLINE float Ceil(const float v) {
    return ceilf(v);
}


/**
 * Rounds-up a given number.
 */
static FORCE_INLINE double Ceil(const double v) {
    return ceil(v);
}


/**
 * Rounds-down a given number.
 */
static FORCE_INLINE float Floor(const float v) {
    return floorf(v);
}


/**
 * Rounds-down a given number.
 */
static FORCE_INLINE double Floor(const double v) {
    return floor(v);
}
