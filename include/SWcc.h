#ifndef SWCC_H
#define SWCC_H
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

typedef uint16_t bitmap_t;

static inline void clear_bit(bitmap_t *bitmap, int bit_index) {
    *bitmap &= ~((bitmap_t)1 << bit_index);
}

static inline void clear_all_bits(bitmap_t *bitmap) {
    *bitmap = 0;
}

static inline bool is_bit_set(const bitmap_t *bitmap, int bit_index) {
    return ((*bitmap) & ((bitmap_t)1 << bit_index)) != 0;
}

#endif