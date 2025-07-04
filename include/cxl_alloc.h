#ifndef CXL_ALLOC_H
#define CXL_ALLOC_H
#pragma once
#include <stddef.h>

void* cxl_alloc(size_t size, size_t align);

void cxl_free(void* ptr);

void force_read_from_mem(void* ptr);

void force_write_to_mem(void* ptr);

#endif