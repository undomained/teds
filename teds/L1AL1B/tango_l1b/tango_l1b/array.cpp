// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "header.h"
#include "functions.h"
#include "array.h"

namespace tango {

template<class T>
void Array<T>::partition( // {{{
    int lo,
    int hi,
    T *arr,
    int &p,
    size_t *indexarray
)
{

    T pivot = arr[(lo+hi)/2];
    // Iterators from low and high sides.
    int ilo = lo-1;
    int ihi = hi+1;
    while (true) {
        // Progress i until item needs to be swapped.
        // With the strict < sign, the pivot itself will always prevent running off the
        // entire array.
        ilo++;
        while (arr[ilo] < pivot) ilo++;
        // Repeat joke on other side.
        ihi--;
        while (arr[ihi] > pivot) ihi--;
        // If the array is in perfect order, ilo should be one higher than ihi.
        // If the indices are equal, it must be the pivot number, because it stopped
        // both of the progressors and it is fine to put partitioning index p that
        // that element.
        if (ilo >= ihi) {
            p = ihi; // The lower of the two.
            return;
        }
        // Index ihi is still higher than ilo. We can swap elements and progress further.
        // Swap elements.
        T swap = arr[ilo];
        arr[ilo] = arr[ihi];
        arr[ihi] = swap;
        if (indexarray != 0) {
            size_t idx_swap = indexarray[ilo];
            indexarray[ilo] = indexarray[ihi];
            indexarray[ihi] = idx_swap;
        }
    }

} // }}}
template<class T>
void Array<T>::bruteforce_sort( // {{{
    int lo,
    int hi,
    T *arr,
    size_t *indexarray
)
{
    for (int i=lo+1 ; i<=hi ; i++) {
        T cur = arr[i];
        size_t idx_cur;
        if (indexarray != 0) idx_cur = indexarray[i];
        // We have i sorted numbers on the left. Element i has to progress to the
        // left as far as it should go.
        int testel = i-1; // Last of the sorted array.
        while (true) {
            if (arr[testel] > cur) {
                // Move the original element a step to the right.
                arr[testel+1] = arr[testel];
                if (indexarray != 0) indexarray[testel+1] = indexarray[testel];
                testel--; // Check next element.
            } else {
                arr[testel+1] = cur;
                if (indexarray != 0) indexarray[testel+1] = idx_cur;
                break;
            }
            if (testel < lo) {
                arr[lo] = cur;
                if (indexarray != 0) indexarray[lo] = idx_cur;
                break;
            }
        }
    }
} // }}}
template<class T>
void Array<T>::execute_sort( // {{{
    int lo,
    int hi,
    T *arr,
    size_t *indexarray
)
{
    if (hi - lo + 1 <= SORT_BRUTEFORCE_THRESHOLD) {
        bruteforce_sort(lo,hi,arr,indexarray);
    } else {
        int p;
        partition(lo,hi,arr,p,indexarray);
        execute_sort(lo,p,arr,indexarray);
        execute_sort(p+1,hi,arr,indexarray);
    }

} // }}}
template<class T>
void Array<T>::sort( // {{{
    size_t sz,
    const T *arr,
    T *a_sorted,
    size_t *indexarray
)
{
    // Copy array to sorted.
    vector<T> sorted_new;
    T* sorted;
    if (a_sorted == 0) {
        // No sorted array. Use new one.
        sorted_new.resize(sz);
        sorted = sorted_new.data();
    } else {
        // Use argument.
        sorted = a_sorted;
    }
    for (size_t i=0 ; i<sz ; i++) sorted[i] = arr[i];
    // Initialize index array.
    if (indexarray != 0) {
        for (size_t i=0 ; i<sz ; i++) indexarray[i] = i;
    }
    execute_sort(0,sz-1,sorted,indexarray);
} // }}}
template<class T>
void Array<T>::execute_select( // {{{
    int lo,
    int hi,
    int tar_lo,
    int tar_hi,
    T *arr,
    T &res
)
{
    if (tar_lo == lo && tar_hi == hi) {
        // The entire array. No need for further sorting.
        res = 0.0;
        for (int tar = tar_lo ; tar <= tar_hi ; tar++) res += arr[tar];
    } else if (hi - lo + 1 <= SORT_BRUTEFORCE_THRESHOLD) {
        // Small array. Apply brute-force sorting.
        bruteforce_sort(lo,hi,arr);
        res = 0.0;
        for (int tar = tar_lo ; tar <= tar_hi ; tar++) res += arr[tar];
    } else {
        int p;
        partition(lo,hi,arr,p);
        if (tar_hi <= p) { // Everything on the left.
            execute_select(lo,p,tar_lo,tar_hi,arr,res);
        } else if (tar_lo > p) { // Everything on the right.
            execute_select(p+1,hi,tar_lo,tar_hi,arr,res);
        } else {
            // Split left and right.
            T res1;
            T res2;
            execute_select(lo,p,tar_lo,p,arr,res1);
            execute_select(p+1,hi,p+1,tar_hi,arr,res2);
            res = res1 + res2;
        }
    }
} // }}}
template<class T>
void Array<T>::select( // {{{
    int sz,
    int tar_lo,
    int tar_hi,
    T *arr,
    T &res
)
{
    execute_select(0,sz-1,tar_lo,tar_hi,arr,res);
    res /= (tar_hi-tar_lo+1);
} // }}}

} // namespace tango
