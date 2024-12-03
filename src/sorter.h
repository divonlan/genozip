// ------------------------------------------------------------------
//   sorter.h
//   Copyright (C) 2019-2024 Genozip Limited. Patent Pending.
//   Please see terms and conditions in the file LICENSE.txt
//
//   WARNING: Genozip is proprietary, not open source software. Modifying the source code is strictly prohibited
//   and subject to penalties specified in the license.

#pragma once

#include "buffer.h"

// used for qsort sort function - receives two integers of any type and returns -1/0/1 as required to sort in ascending order
#define SORTER(func) int func (const void *a, const void *b)
typedef SORTER ((*Sorter));
#define ASCENDING_RAW(a,b) (((a) > (b)) ? 1 : (a) < (b) ? -1 : 0)
#define DESCENDING_RAW(a,b) (-ASCENDING_RAW((a), (b)))
#define ASCENDING(struct_type,struct_field) ASCENDING_RAW (((struct_type *)a)->struct_field, ((struct_type *)b)->struct_field)
#define DESCENDING(struct_type,struct_field) (-ASCENDING(struct_type, struct_field))

#define ASCENDING_SORTER(func_name,struct_type,struct_field) \
    SORTER (func_name) { return ASCENDING (struct_type, struct_field); }

#define DESCENDING_SORTER(func_name,struct_type,struct_field) \
    SORTER (func_name) { return DESCENDING (struct_type, struct_field); }

// declaration of binary search recursive function - array of struct, must be sorted ascending by field
typedef enum { IfNotExact_ReturnLower, IfNotExact_ReturnHigher, IfNotExact_ReturnNULL, IfNotExact_Error } BinarySearchIfNotExact; 
#define BINARY_SEARCHER(func, struct_type, field_type, struct_field,                                                    \
                        is_unique,    /* if false, we return the first entry of requested value */                      \
                        if_not_exact) /* what to do if exact match is not found. Note: if IfNotExact_ReturnLower, and value is less than the first, will return pointer to an item before the first. Likewise with IfNotExact_ReturnHigher */ \
    struct_type *func (struct_type *first, struct_type *last, field_type value, int level)                              \
    {                                                                                                                   \
        ASSERT (level < 50, "recursion too deep (level=50) first=%p last=%p", first, last);                             \
        if (first > last) { /* not found */                                                                             \
            ASSERT (if_not_exact!=IfNotExact_Error, "requested %s=%"PRIu64" not found", #struct_field, (uint64_t)value);\
            return if_not_exact==IfNotExact_ReturnNULL?NULL : if_not_exact==IfNotExact_ReturnLower?last : first;        \
        }                                                                                                               \
        struct_type *mid = first + (last - first) / 2;                                                                  \
        if (mid->struct_field == value)                                                                                 \
            return (is_unique || mid == first || (mid->struct_field != (mid-1)->struct_field))                          \
                ? mid /* found and element is known to be unique, or is verified to be the first */                     \
                : func (first, mid-1, value, level+1); /* found, but it isn't first: continue searching for the first */\
        return (mid->struct_field > value) ? func (first, mid-1, value, level+1) : func (mid+1, last, value, level+1);  \
    }

// declaration of binary search recursive function - array of integral (i.e. comparable) values, must be sorted ascending
#define BINARY_SEARCHER_INTEGRAL(func, element_type,                                                                    \
                                 is_unique,    /* if false, we return the first entry of requested value*/              \
                                 if_not_exact) /* what to do if exact match is not found. Note: if IfNotExact_ReturnLower, and value is less than the first, will return pointer to an item before the first. Likewise with IfNotExact_ReturnHigher */ \
    element_type *func (element_type *first, element_type *last, element_type value, int level)                         \
    {                                                                                                                   \
        if (first > last) /* not found */                                                                               \
            return if_not_exact==IfNotExact_ReturnNULL?NULL : if_not_exact==IfNotExact_ReturnLower?last : first;        \
        element_type *mid = first + (last - first) / 2;                                                                 \
        if (*mid == value)                                                                                              \
            return (is_unique || mid == first || (*mid != *(mid-1)))                                                    \
                ? mid /* found and element is known to be unique, or is verified to be the first */                     \
                : func (first, mid-1, value, level+1); /* found, but it isn't first: continue searching for the first */\
        return (*mid > value) ? func (first, mid-1, value, level+1) : func (mid+1, last, value, level+1);               \
    }

// actually do a binary search. buf is required to be sorted in an ascending order of field
#define binary_search(func, type, buf, value) ((buf).len ? func (B1ST(type,(buf)), BLST(type,(buf)), value, 0) : NULL)

// declaration of binary search recursive function: sorted index array into an unsorted data array
// index_arr (an array of uint32_t) contains indices into data_arr (an array of struct_type).
// index_arr sorted-ascending by struct_field (a field of struct_type), while data_arr is not sorted.
#define BINARY_SEARCHER_WITH_INDEX(func, struct_type, field_type, struct_field,                                         \
                        is_unique,    /* true if data_arr is guaranteed to contain unique values of struct_field. if false, we return the first entry of requested value */                      \
                        if_not_exact) /* what to do if exact match is not found. Note: if IfNotExact_ReturnLower, and value is less than the first, will return pointer to an item before the first. Likewise with IfNotExact_ReturnHigher */ \
    int32_t func (int32_t first, int32_t last, field_type value, struct_type *data_arr, uint32_t *index_arr, int level) \
    {                                                                                                                   \
        if (first > last) /* not found */                                                                               \
            return if_not_exact==IfNotExact_ReturnNULL?-1 : if_not_exact==IfNotExact_ReturnLower?last : first;          \
        int32_t mid = first + (last - first) / 2;                                                                       \
        field_type mid_value = data_arr[index_arr[mid]].struct_field;                                                   \
        if (mid_value == value)                                                                                         \
            return (is_unique || mid == first || (mid_value != data_arr[index_arr[mid-1]].struct_field))                \
                ? mid /* found and element is known to be unique, or is verified to be the first */                     \
                : func (first, mid-1, value, data_arr, index_arr, level+1); /* found, but it isn't first: continue searching for the first */\
        return (mid_value > value) ? func (first, mid-1, value, data_arr, index_arr, level+1)                           \
                                   : func (mid+1, last,  value, data_arr, index_arr, level+1);                          \
    }

#define binary_search_with_index(func, struct_type, data_buf, index_buf, value, index_i) \
    ( (data_buf).len ? func (0, (data_buf).len-1, (value), B1ST(struct_type,(data_buf)), B1ST32(index_buf), 0) : -1)                                                                                                       \

