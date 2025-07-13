#ifndef ZIPF_H_
#define ZIPF_H_

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

typedef struct {
    uint64_t n_;        // number of items
    double theta_;      // skewness (0=uniform, >0=zipf)
    double alpha_;      // only depends on theta
    double thres_;      // only depends on theta
    uint64_t last_n_;   // last n used to calculate zetan_
    double dbl_n_;
    double zetan_;
    double eta_;
    uint64_t seq_;      // for sequential number generation
    uint64_t rand_seed_;
} zipf_gen_t;

// Function declarations
zipf_gen_t* zipf_create(uint64_t n, double theta, uint64_t rand_seed);
void zipf_destroy(zipf_gen_t* zg);
uint64_t zipf_next(zipf_gen_t* zg);
void zipf_change_n(zipf_gen_t* zg, uint64_t n);
void zipf_test(double theta);

#endif // ZIPF_H_