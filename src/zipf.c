#include "zipf.h"

// Simple random number generator
static double zipf_rand_f64(uint64_t* seed) {
    *seed = (*seed * 1103515245ULL + 12345ULL) & 0x7fffffffULL;
    return (double)(*seed) / (double)0x7fffffffULL;
}

// Fast power approximation
static double pow_approx(double a, double b) {
    // From http://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/
    int e = (int)b;
    union {
        double d;
        int x[2];
    } u = {a};
    u.x[1] = (int)((b - (double)e) * (double)(u.x[1] - 1072632447) + 1072632447.);
    u.x[0] = 0;

    // Fast exponentiation
    double r = 1.;
    while (e) {
        if (e & 1) r *= a;
        a *= a;
        e >>= 1;
    }

    return r * u.d;
}

// Calculate zeta function
static double zeta(uint64_t last_n, double last_sum, uint64_t n, double theta) {
    if (last_n > n) {
        last_n = 0;
        last_sum = 0.;
    }
    while (last_n < n) {
        last_sum += 1. / pow_approx((double)last_n + 1., theta);
        last_n++;
    }
    return last_sum;
}

// Create Zipf generator
zipf_gen_t* zipf_create(uint64_t n, double theta, uint64_t rand_seed) {
    assert(n > 0);
    
    if (theta > 0.992 && theta < 1)
        fprintf(stderr, "warning: theta > 0.992 will be inaccurate due to approximation\n");
    
    if (theta >= 1. && theta < 40.) {
        fprintf(stderr, "error: theta in [1., 40.) is not supported\n");
        return NULL;
    }
    
    assert(theta == -1. || (theta >= 0. && theta < 1.) || theta >= 40.);
    
    zipf_gen_t* zg = (zipf_gen_t*)malloc(sizeof(zipf_gen_t));
    if (!zg) return NULL;
    
    zg->n_ = n;
    zg->theta_ = theta;
    zg->rand_seed_ = rand_seed;
    
    if (theta == -1.) {
        zg->seq_ = rand_seed % n;
        zg->alpha_ = 0;
        zg->thres_ = 0;
    } else if (theta > 0. && theta < 1.) {
        zg->seq_ = 0;
        zg->alpha_ = 1. / (1. - theta);
        zg->thres_ = 1. + pow_approx(0.5, theta);
    } else {
        zg->seq_ = 0;
        zg->alpha_ = 0.;
        zg->thres_ = 0.;
    }
    
    zg->last_n_ = 0;
    zg->zetan_ = 0.;
    zg->eta_ = 0;
    
    return zg;
}

// Destroy Zipf generator
void zipf_destroy(zipf_gen_t* zg) {
    if (zg) {
        free(zg);
    }
}

// Generate next Zipf distributed number
uint64_t zipf_next(zipf_gen_t* zg) {
    if (!zg) return 0;
    
    if (zg->last_n_ != zg->n_) {
        if (zg->theta_ > 0. && zg->theta_ < 1.) {
            zg->zetan_ = zeta(zg->last_n_, zg->zetan_, zg->n_, zg->theta_);
            zg->eta_ = (1. - pow_approx(2. / (double)zg->n_, 1. - zg->theta_)) /
                       (1. - zeta(0, 0., 2, zg->theta_) / zg->zetan_);
        }
        zg->last_n_ = zg->n_;
        zg->dbl_n_ = (double)zg->n_;
    }

    if (zg->theta_ == -1.) {
        uint64_t v = zg->seq_;
        if (++(zg->seq_) >= zg->n_) zg->seq_ = 0;
        return v;
    } else if (zg->theta_ == 0.) {
        double u = zipf_rand_f64(&zg->rand_seed_);
        return (uint64_t)(zg->dbl_n_ * u);
    } else if (zg->theta_ >= 40.) {
        return 0UL;
    } else {
        // Algorithm from J. Gray et al.
        double u = zipf_rand_f64(&zg->rand_seed_);
        double uz = u * zg->zetan_;
        if (uz < 1.)
            return 0UL;
        else if (uz < zg->thres_)
            return 1UL;
        else {
            uint64_t v = (uint64_t)(zg->dbl_n_ * pow_approx(zg->eta_ * (u - 1.) + 1., zg->alpha_));
            if (v >= zg->n_) v = zg->n_ - 1;
            return v;
        }
    }
}

// Change n value
void zipf_change_n(zipf_gen_t* zg, uint64_t n) {
    if (zg) {
        zg->n_ = n;
    }
}

// Test function
void zipf_test(double theta) {
    double zetan = 0.;
    const uint64_t n = 1000000UL;
    uint64_t i;

    for (i = 0; i < n; i++) 
        zetan += 1. / pow((double)i + 1., theta);

    if (theta < 1. || theta >= 40.) {
        zipf_gen_t* zg = zipf_create(n, theta, 0);
        if (!zg) return;

        uint64_t num_key0 = 0;
        const uint64_t num_samples = 10000000UL;
        for (i = 0; i < num_samples; i++)
            if (zipf_next(zg) == 0) num_key0++;

        printf("theta = %lf; using pow(): %.10lf", theta, 1. / zetan);
        printf(", using approx-pow(): %.10lf", (double)num_key0 / (double)num_samples);
        printf("\n");
        
        zipf_destroy(zg);
    }
}