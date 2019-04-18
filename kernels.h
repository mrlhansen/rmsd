#ifndef KERNELS_H
#define KERNELS_H

void delta0(mpfr_t, mpfr_t, mpfr_t);
void delta1(mpfr_t, mpfr_t, mpfr_t);
void delta2(mpfr_t, mpfr_t, mpfr_t);
void delta3(mpfr_t, mpfr_t, mpfr_t);

void int_delta0_sq(mpfr_t, mpfr_t, mpfr_t);
void int_delta1_sq(mpfr_t, mpfr_t, mpfr_t);
void int_delta2_sq(mpfr_t, mpfr_t, mpfr_t);
void int_delta3_sq(mpfr_t, mpfr_t, mpfr_t);

void int_delta0(mpfr_t, mpfr_t, mpfr_t, mpfr_t);
void int_delta1(mpfr_t, mpfr_t, mpfr_t, mpfr_t);
void int_delta2(mpfr_t, mpfr_t, mpfr_t, mpfr_t);
void int_delta3(mpfr_t, mpfr_t, mpfr_t, mpfr_t);

void init_kernel(double);

#endif
