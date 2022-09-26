
#pragma once


#include <sys/types.h>


// In-place quicksort implementation.
template <typename T>
static inline void __quicksort_swap(T *array, size_t first, size_t second) {
    T old = array[first];
    array[first] = array[second];
    array[second] = old;
}


template <typename T, typename Compare>
static void __quicksort_it(T *array, size_t beginning, size_t end,
    const Compare compare)
{
    if (end > (beginning + 1)) {
        const T &piv = array[beginning];
        size_t l = beginning + 1;
        size_t r = end;
        while (l < r) {
            if (compare(array[l], piv)) {
                l++;
            } else {
                __quicksort_swap(array, l, --r);
            }
        }

        __quicksort_swap(array, --l, beginning);
        __quicksort_it(array, beginning, l, compare);
        __quicksort_it(array, r, end, compare);
    }
}


template <typename T, typename Compare>
static inline void Quicksort(T *array, size_t size, const Compare compare) {
    __quicksort_it(array, 0, size, compare);
}
