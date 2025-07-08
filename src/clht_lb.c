/*
 *   File: clht_lb.c
 *   Author: Vasileios Trigonakis <vasileios.trigonakis@epfl.ch>
 *   Description: lock-based cache-line hash table with no resizing
 *   clht_lb.c is part of ASCYLIB
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2014 Vasileios Trigonakis <vasileios.trigonakis@epfl.ch>
 *	      	      Distributed Programming Lab (LPD), EPFL
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 *of this software and associated documentation files (the "Software"), to deal
 *in the Software without restriction, including without limitation the rights
 *to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *copies of the Software, and to permit persons to whom the Software is
 *furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 *all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *SOFTWARE.
 *
 */
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <xmmintrin.h>
#include <stdatomic.h>

#include "clht_lb.h"
#include "cxl_alloc.h"
#include "SWcc.h"
#include "clht_thread_local.h"

__thread ssmem_allocator_t *clht_alloc;

#ifdef DEBUG
__thread uint32_t put_num_restarts = 0;
__thread uint32_t put_num_failed_expand = 0;
__thread uint32_t put_num_failed_on_new = 0;
#endif

#include "assert.h"
#include "stdlib.h"

const char *clht_type_desc() { return "CLHT-LB-NO-RESIZE"; }

inline int is_power_of_two(unsigned int x) {
    return ((x != 0) && !(x & (x - 1)));
}

static inline int is_odd(int x) { return x & 1; }

/** Jenkins' hash function for 64-bit integers. */
inline uint64_t __ac_Jenkins_hash_64(uint64_t key) {
    key += ~(key << 32);
    key ^= (key >> 22);
    key += ~(key << 13);
    key ^= (key >> 8);
    key += (key << 3);
    key ^= (key >> 15);
    key += ~(key << 27);
    key ^= (key >> 31);
    return key;
}

/* Create a new bucket. */
bucket_t *clht_bucket_create() {
    bucket_t *bucket = NULL;
    bucket = (bucket_t *)cxl_alloc(sizeof(bucket_t), CACHE_LINE_SIZE);
    /* bucket = malloc(sizeof(bucket_t)); */
    if (bucket == NULL) {
        return NULL;
    }

    bucket->lock = 0;

    uint32_t j;
    for (j = 0; j < ENTRIES_PER_BUCKET; j++) {
        bucket->key[j] = 0;
    }
    bucket->next = NULL;

    return bucket;
}

clht_t *clht_create(uint64_t num_buckets) {
    // clht_t* w = (clht_t*) memalign(CACHE_LINE_SIZE, sizeof(clht_t));
    clht_t *w = (clht_t *)cxl_alloc(sizeof(clht_t), CACHE_LINE_SIZE);
    if (w == NULL) {
        printf("** malloc @ hatshtalbe\n");
        return NULL;
    }

    w->ht = clht_hashtable_create(num_buckets);
    if (w->ht == NULL) {
        free(w);
        return NULL;
    }

    return w;
}

clht_hashtable_t *clht_hashtable_create(uint64_t num_buckets) {
    clht_hashtable_t *hashtable = NULL;

    if (num_buckets == 0) {
        return NULL;
    }

    /* Allocate the table itself. */
    // hashtable = (clht_hashtable_t*) memalign(CACHE_LINE_SIZE,
    // sizeof(clht_hashtable_t));
    hashtable = (clht_hashtable_t *)cxl_alloc(sizeof(clht_hashtable_t),
                                              CACHE_LINE_SIZE);
    if (hashtable == NULL) {
        printf("** malloc @ hatshtalbe\n");
        return NULL;
    }

    /* hashtable->table = calloc(num_buckets, (sizeof(bucket_t))); */
    // hashtable->table = (bucket_t*) memalign(CACHE_LINE_SIZE, num_buckets *
    // (sizeof(bucket_t)));
    hashtable->table =
        (bucket_t *)cxl_alloc(num_buckets * sizeof(bucket_t), CACHE_LINE_SIZE);
    hashtable->bitmap = 
        (bitmap_t *)cxl_alloc(num_buckets * sizeof(bitmap_t), CACHE_LINE_SIZE);
    if (hashtable->table == NULL) {
        printf("** alloc: hashtable->table\n");
        fflush(stdout);
        free(hashtable);
        return NULL;
    }

    memset(hashtable->table, 0, num_buckets * (sizeof(bucket_t)));
    memset(hashtable->bitmap, 0, num_buckets * sizeof(uint16_t));

    uint64_t i;
    for (i = 0; i < num_buckets; i++) {
        hashtable->table[i].lock = 0;
        hashtable->bitmap[i] = 0;
        uint32_t j;
        for (j = 0; j < ENTRIES_PER_BUCKET; j++) {
            hashtable->table[i].key[j] = 0;
        }
    }

    hashtable->num_buckets = num_buckets;

    return hashtable;
}

/* Hash a key for a particular hash table. */
uint64_t clht_hash(clht_hashtable_t *hashtable, clht_addr_t key) {
    /* uint64_t hashval; */
    /* hashval = __ac_Jenkins_hash_64(key); */
    /* return hashval % hashtable->num_buckets; */
    /* return key % hashtable->num_buckets; */
    return key & (hashtable->num_buckets - 1);
}

/* Retrieve a key-value entry from a hash table. */
clht_val_t clht_get(clht_hashtable_t *hashtable, clht_addr_t key) {
    size_t bin = clht_hash(hashtable, key);
    volatile bucket_t *bucket = hashtable->table + bin;

    force_read_from_mem((void*)&hashtable->bitmap[bin]);
    bool should_flush = !is_bit_set(&hashtable->bitmap[bin], thread_id);
    volatile bucket_t* cur = bucket;
    if (should_flush) {
        while (cur != NULL) {
            force_read_from_mem((void*)cur);
            cur = cur->next;
        }
        atomic_fetch_or((_Atomic bitmap_t*)&hashtable->bitmap[bin], ((bitmap_t)1 << thread_id));
    }

    uint32_t j;
    do {
        for (j = 0; j < ENTRIES_PER_BUCKET; j++) {
            clht_val_t val = bucket->val[j];
#ifdef __tile__
            _mm_lfence();
#endif
            if (bucket->key[j] == key) {
                return val;
            }
        }

        bucket = bucket->next;
    } while (bucket != NULL);
    return 0;
}

inline clht_addr_t bucket_exists(bucket_t *bucket, clht_addr_t key) {
    force_read_from_mem((void*)bucket);
    uint32_t j;
    do {
        for (j = 0; j < ENTRIES_PER_BUCKET; j++) {
            if (bucket->key[j] == key) {
                return true;
            }
        }

        bucket = bucket->next;
    } while (bucket != NULL);
    return false;
}

/* Insert a key-value entry into a hash table. */
int clht_put(clht_t *h, clht_addr_t key, clht_val_t val) {
    clht_hashtable_t *hashtable = h->ht;
    size_t bin = clht_hash(hashtable, key);
    bucket_t *bucket = hashtable->table + bin;

#if defined(READ_ONLY_FAIL)
    if (bucket_exists(bucket, key)) {
        return false;
    }
#endif
    clht_lock_t* lock = &bucket->lock;
    clht_addr_t *empty = NULL;
    clht_val_t *empty_v = NULL;

    uint32_t j;
    LOCK_ACQ(lock);
    do {
        
        for (j = 0; j < ENTRIES_PER_BUCKET; j++) {
            if (bucket->key[j] == key) {
                LOCK_RLS(lock);
                return false;
            } else if (empty == NULL && bucket->key[j] == 0) {
                empty = &bucket->key[j];
                empty_v = &bucket->val[j];
            }
        }

        if (bucket->next == NULL) {
            if (empty == NULL) {
                DPP(put_num_failed_expand);
                bucket->next = clht_bucket_create();
                bucket->next->key[0] = key;
#ifdef __tile__
                _mm_sfence();
#endif
                bucket->next->val[0] = val;
            } else {
                *empty_v = val;
#ifdef __tile__
                _mm_sfence();
#endif
                *empty = key;
            }
            force_write_to_mem((void*)bucket);
            _mm_sfence();
            atomic_store((_Atomic bitmap_t*)&hashtable->bitmap[bin], ((bitmap_t)1 << thread_id));
            LOCK_RLS(lock);
            return true;
        }

        bucket = bucket->next;
    } while (true);
}

/* Remove a key-value entry from a hash table. */
clht_val_t clht_remove(clht_t *h, clht_addr_t key) {
    clht_hashtable_t *hashtable = h->ht;
    size_t bin = clht_hash(hashtable, key);
    bucket_t *bucket = hashtable->table + bin;

#if defined(READ_ONLY_FAIL)
    if (!bucket_exists(bucket, key)) {
        return false;
    }
#endif /* READ_ONLY_FAIL */

    clht_lock_t *lock = &bucket->lock;
    uint32_t j;

    LOCK_ACQ(lock);
    do {
        for (j = 0; j < ENTRIES_PER_BUCKET; j++) {
            if (bucket->key[j] == key) {
                clht_val_t val = bucket->val[j];
                bucket->key[j] = 0;
                force_write_to_mem((void*)bucket);
                _mm_sfence();
                atomic_store((_Atomic bitmap_t*)&hashtable->bitmap[bin], ((bitmap_t)1 << thread_id));
                LOCK_RLS(lock);
                return val;
            }
        }
        bucket = bucket->next;
    } while (bucket != NULL);
    LOCK_RLS(lock);
    return false;
}

static uint32_t clht_put_seq(clht_hashtable_t *hashtable, clht_addr_t key,
                             clht_val_t val, uint64_t bin) {
    bucket_t *bucket = hashtable->table + bin;
    clht_addr_t *empty = NULL;
    clht_val_t *empty_v = NULL;

    uint32_t j;

    do {
        for (j = 0; j < ENTRIES_PER_BUCKET; j++) {
            if (bucket->key[j] == key) {
                return false;
            } else if (empty == NULL && bucket->key[j] == 0) {
                empty = &bucket->key[j];
                empty_v = &bucket->val[j];
            }
        }

        if (bucket->next == NULL) {
            if (empty == NULL) {
                DPP(put_num_failed_expand);
                bucket->next = clht_bucket_create();
                bucket->next->key[0] = key;
                bucket->next->val[0] = val;
            } else {
                *empty_v = val;
                *empty = key;
            }
            return true;
        }
        
		force_write_to_mem((void*)&bucket);
        _mm_sfence();
        atomic_store((_Atomic bitmap_t*)&hashtable->bitmap[bin], ((bitmap_t)1 << thread_id));
        bucket = bucket->next;
    } while (true);
}

static inline void bucket_cpy(bucket_t *bucket, clht_hashtable_t *ht_new) {
    LOCK_ACQ(&bucket->lock);
    uint32_t j;
    do {
        for (j = 0; j < ENTRIES_PER_BUCKET; j++) {
            clht_addr_t key = bucket->key[j];
            if (key != 0) {
                uint64_t bin = clht_hash(ht_new, key);
                clht_val_t val = bucket->key[j];
                clht_put_seq(ht_new, key, val, bin);
            }
        }
        bucket = bucket->next;
    } while (bucket != NULL);
    LOCK_RLS(&bucket->lock);
}

void clht_destroy(clht_hashtable_t *hashtable) {
    cxl_free(hashtable->table);
    cxl_free(hashtable);
}

size_t clht_size(clht_hashtable_t *hashtable) {
    uint64_t num_buckets = hashtable->num_buckets;
    bucket_t *bucket = NULL;
    size_t size = 0;

    uint64_t bin;
    for (bin = 0; bin < num_buckets; bin++) {
        bucket = hashtable->table + bin;
        force_read_from_mem((void*)&bucket);
        uint32_t j;
        do {
            for (j = 0; j < ENTRIES_PER_BUCKET; j++) {
                if (bucket->key[j] > 0) {
                    size++;
                }
            }

            bucket = bucket->next;
        } while (bucket != NULL);
    }
    return size;
}

void clht_print(clht_hashtable_t *hashtable) {
    uint64_t num_buckets = hashtable->num_buckets;
    bucket_t *bucket;

    printf("Number of buckets: %u\n", num_buckets);

    uint64_t bin;
    for (bin = 0; bin < num_buckets; bin++) {
        bucket = hashtable->table + bin;

        printf("[[%05d]] ", bin);

        uint32_t j;
        do {
            for (j = 0; j < ENTRIES_PER_BUCKET; j++) {
                if (bucket->key[j]) {
                    printf("(%-5llu)-> ",
                           (long long unsigned int)bucket->key[j]);
                }
            }

            bucket = bucket->next;
            printf(" ** -> ");
        } while (bucket != NULL);
        printf("\n");
    }
    fflush(stdout);
}
