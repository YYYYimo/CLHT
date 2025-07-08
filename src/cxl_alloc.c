#include "cxl_alloc.h"

#include <emmintrin.h>
#include <immintrin.h>
#include <stdlib.h>
#include <malloc.h>
#include <xmmintrin.h>

void* cxl_alloc(size_t size, size_t align) {
    void* ptr = NULL;
    if (align <= 8) {
        ptr = malloc(size);
    } else {
        ptr = memalign(align, size);
    }
    return ptr;
}

void cxl_free(void* ptr) {
    free(ptr);
}

void force_read_from_mem(void* ptr) {
    _mm_clflush(ptr);
}

void force_write_to_mem(void* ptr) {
    _mm_clwb(ptr);
}