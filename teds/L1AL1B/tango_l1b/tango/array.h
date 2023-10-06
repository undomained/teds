// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#ifndef ARRAY_H
#define ARRAY_H

#include "header.h"

#define SORT_BRUTEFORCE_THRESHOLD 7

// This class has no objects. The constructor and destructor are private and even
// their own methods do not create object. All member functions are static, so it
// is as if there is no class. The only thing the class does is to wrap up the
// functions for the template.
template <typename T>
class Array {
    // Private constructor and destructor. No instances of this class are created,
    // it is a module with static functions. Use Array<TYPE>::FUNCTION to call
    // a function.
    private:
    Array() {}
    ~Array() {}
    public:
    // Sort interface. Turns array into sorted array and/or creates index array.
    // The idea is that sorted[i] = unsorted[indexarray[i]].
    static void sort(
        size_t sz,
        const T *arr,
        T *a_sorted,
        size_t *indexarray=0
    );
    static void select(
        int sz,
        int tar_lo,
        int tar_hi,
        T *arr,
        T &res
    );
    private:
    static void partition(
        int lo,
        int hi,
        T *arr,
        int &p,
        size_t *indexarray=0
    );
    static void bruteforce_sort(
        int lo,
        int hi,
        T *arr,
        size_t *indexarray=0
    );
    static void execute_sort(
        int lo,
        int hi,
        T *arr,
        size_t *indexarray=0
    );
    static void execute_select(
        int lo,
        int hi,
        int tar_lo,
        int tar_hi,
        T *arr,
        T &res
    );
};

// Instances that must be known at link time.
// It seems it works if these are in the C++ file or in the H-file.
template class Array<double>;
template class Array<uint64_t>;

#endif
