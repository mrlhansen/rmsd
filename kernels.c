#include <mpfr.h>
#include "corefcts.h"

#include <arb.h>
#include <acb.h>
#include <acb_hypgeom.h>
#include <arb_hypgeom.h>

static mpfr_t pi;
static mpfr_t sqrt2;
static mpfr_t sigma;

void delta0(mpfr_t res, mpfr_t estar, mpfr_t omega)
{
	mpfr_t norm, arg;
	mpfr_inits(norm, arg, NULL);

	mpfr_mul_d(norm, pi, 2.0, ROUNDING);
	mpfr_sqrt(norm, norm, ROUNDING);
	mpfr_mul(norm, norm, sigma, ROUNDING);

	mpfr_sub(arg, estar, omega, ROUNDING);
	mpfr_div(arg, arg, sigma, ROUNDING);
	mpfr_mul(arg, arg, arg, ROUNDING);
	mpfr_mul_d(arg, arg, -0.5, ROUNDING);

	mpfr_exp(res, arg, ROUNDING);
	mpfr_div(res, res, norm, ROUNDING);

	mpfr_mul(norm, sqrt2, sigma, ROUNDING);
	mpfr_div(norm, estar, norm, ROUNDING);
	mpfr_erf(norm, norm, ROUNDING);
	mpfr_add_d(norm, norm, 1.0, ROUNDING);
	mpfr_mul_d(norm, norm, 0.5, ROUNDING);
	mpfr_div(res, res, norm, ROUNDING);

	mpfr_clears(norm, arg, NULL);
}

void delta1(mpfr_t res, mpfr_t estar, mpfr_t omega)
{
	mpfr_t norm, arg;
	mpfr_inits(norm, arg, NULL);

	if(mpfr_cmp(estar, omega) == 0)
	{
		mpfr_d_div(res, 1.0, sigma, ROUNDING);
		return;
	}

	mpfr_sub(norm, estar, omega, ROUNDING);
	mpfr_mul(norm, norm, pi, ROUNDING);
	mpfr_div(arg, norm, sigma, ROUNDING);
	mpfr_sin(res, arg, ROUNDING);
	mpfr_div(res, res, norm, ROUNDING);

	mpfr_clears(norm, arg, NULL);
}

void delta2(mpfr_t res, mpfr_t estar, mpfr_t omega)
{
	mpfr_t atmp, btmp;
	mpfr_inits(atmp, btmp, NULL);

	mpfr_sub(res, estar, omega, ROUNDING);
	mpfr_sqr(atmp, res, ROUNDING);
	mpfr_sqr(btmp, sigma, ROUNDING);
	mpfr_add(atmp, atmp, btmp, ROUNDING);
	mpfr_div(res, sigma, atmp, ROUNDING);

	mpfr_clears(atmp, btmp, NULL);
}

void delta3(mpfr_t res, mpfr_t estar, mpfr_t omega)
{
	mpfr_t atmp, btmp;
	mpfr_inits(atmp, btmp, NULL);

	mpfr_sub(res, estar, omega, ROUNDING);
	mpfr_sqr(atmp, res, ROUNDING);
	mpfr_sqr(btmp, sigma, ROUNDING);
	mpfr_add(atmp, atmp, btmp, ROUNDING);
	mpfr_div(res, res, atmp, ROUNDING);

	mpfr_clears(atmp, btmp, NULL);
}

void int_delta0_sq(mpfr_t res, mpfr_t estar, mpfr_t e0)
{
	mpfr_t norm, arg;
	mpfr_inits(norm, arg, NULL);

	mpfr_sqrt(norm, pi, ROUNDING);
	mpfr_mul(norm, norm, sigma, ROUNDING);
	mpfr_mul_d(norm, norm, 4.0, ROUNDING);

	mpfr_sub(arg, estar, e0, ROUNDING);
	mpfr_div(arg, arg, sigma, ROUNDING);
	mpfr_erf(res, arg, ROUNDING);
	mpfr_add_d(res, res, 1.0, ROUNDING);
	mpfr_div(res, res, norm, ROUNDING);

	mpfr_mul(norm, sqrt2, sigma, ROUNDING);
	mpfr_div(norm, estar, norm, ROUNDING);
	mpfr_erf(norm, norm, ROUNDING);
	mpfr_add_d(norm, norm, 1.0, ROUNDING);
	mpfr_mul_d(norm, norm, 0.5, ROUNDING);
	mpfr_sqr(norm, norm, ROUNDING);
	mpfr_div(res, res, norm, ROUNDING);

	mpfr_clears(norm, arg, NULL);
}

void int_delta1_sq(mpfr_t res, mpfr_t estar, mpfr_t e0)
{
	mpfr_t arg, etmp;
	arb_t x;
	arf_t y;

	mpfr_inits(arg, etmp, NULL);
	arb_init(x);
	arf_init(y);

	mpfr_sub(etmp, estar, e0, ROUNDING);
	mpfr_mul(arg, etmp, pi, ROUNDING);
	mpfr_div(arg, arg, sigma, ROUNDING);
	mpfr_sin(res, arg, ROUNDING);
	mpfr_div(res, res, pi, ROUNDING);
	mpfr_sqr(res, res, ROUNDING);
	mpfr_div(res, res, etmp, ROUNDING);

	mpfr_d_div(etmp, 0.5, sigma, ROUNDING);
	mpfr_sub(res, etmp, res, ROUNDING);
	mpfr_mul_d(arg, arg, 2.0, ROUNDING);

	arf_set_mpfr(y, arg);
	arb_set_arf(x, y);
	arb_hypgeom_si(x, x, 500);
	arf_get_mpfr(etmp, arb_midref(x), ROUNDING);

	mpfr_div(etmp, etmp, sigma, ROUNDING);
	mpfr_div(etmp, etmp, pi, ROUNDING);
	mpfr_add(res, res, etmp, ROUNDING);

	mpfr_clears(arg, etmp, NULL);
	arb_clear(x);
	arf_clear(y);
}

void int_delta2_sq(mpfr_t res, mpfr_t estar, mpfr_t e0)
{
	mpfr_t atmp, btmp, ctmp;
	mpfr_inits(atmp, btmp, ctmp, NULL);

	mpfr_sub(ctmp, estar, e0, ROUNDING);
	mpfr_div(atmp, ctmp, sigma, ROUNDING);
	mpfr_atan(res, atmp, ROUNDING);
	mpfr_div_d(atmp, pi, 2.0, ROUNDING);
	mpfr_add(res, res, atmp, ROUNDING);
	mpfr_div(res, res, sigma, ROUNDING);

	mpfr_sqr(atmp, ctmp, ROUNDING);
	mpfr_sqr(btmp, sigma, ROUNDING);
	mpfr_add(atmp, atmp, btmp, ROUNDING);
	mpfr_div(atmp, ctmp, atmp, ROUNDING);

	mpfr_add(res, res, atmp, ROUNDING);
	mpfr_div_d(res, res, 2.0, ROUNDING);

	mpfr_clears(atmp, btmp, ctmp, NULL);
}

void int_delta3_sq(mpfr_t res, mpfr_t estar, mpfr_t e0)
{
	mpfr_t atmp, btmp, ctmp;
	mpfr_inits(atmp, btmp, ctmp, NULL);

	mpfr_sub(ctmp, estar, e0, ROUNDING);
	mpfr_div(atmp, ctmp, sigma, ROUNDING);
	mpfr_atan(res, atmp, ROUNDING);
	mpfr_div_d(atmp, pi, 2.0, ROUNDING);
	mpfr_add(res, res, atmp, ROUNDING);
	mpfr_div(res, res, sigma, ROUNDING);

	mpfr_sqr(atmp, ctmp, ROUNDING);
	mpfr_sqr(btmp, sigma, ROUNDING);
	mpfr_add(atmp, atmp, btmp, ROUNDING);
	mpfr_div(atmp, ctmp, atmp, ROUNDING);

	mpfr_sub(res, res, atmp, ROUNDING);
	mpfr_div_d(res, res, 2.0, ROUNDING);

	mpfr_clears(atmp, btmp, ctmp, NULL);
}

void int_delta0(mpfr_t res, mpfr_t theta, mpfr_t estar, mpfr_t e0)
{
	mpfr_t norm, arg, tmp, rep;
	mpfr_inits(norm, arg, tmp, rep, NULL);

	mpfr_mul(rep, sigma, sigma, ROUNDING);
	mpfr_mul(rep, rep, theta, ROUNDING);

	mpfr_sub(tmp, rep, estar, ROUNDING);
	mpfr_sub(tmp, tmp, estar, ROUNDING);
	mpfr_mul(tmp, tmp, theta, ROUNDING);
	mpfr_mul_d(tmp, tmp, 0.5, ROUNDING);
	mpfr_exp(norm, tmp, ROUNDING);
	mpfr_mul_d(norm, norm, 0.5, ROUNDING);

	mpfr_sub(arg, estar, e0, ROUNDING);
	mpfr_sub(arg, arg, rep, ROUNDING);
	mpfr_mul(tmp, sqrt2, sigma, ROUNDING);
	mpfr_div(arg, arg, tmp, ROUNDING);
	mpfr_erf(res, arg, ROUNDING);
	mpfr_add_d(res, res, 1.0, ROUNDING);
	mpfr_mul(res, res, norm, ROUNDING);

	mpfr_mul(norm, sqrt2, sigma, ROUNDING);
	mpfr_div(norm, estar, norm, ROUNDING);
	mpfr_erf(norm, norm, ROUNDING);
	mpfr_add_d(norm, norm, 1.0, ROUNDING);
	mpfr_mul_d(norm, norm, 0.5, ROUNDING);
	mpfr_div(res, res, norm, ROUNDING);

	mpfr_clears(norm, arg, tmp, rep, NULL);
}

void int_delta1(mpfr_t res, mpfr_t theta, mpfr_t estar, mpfr_t e0)
{
	mpfr_t re, im;
	acb_t z, s;
	arb_t x, y;
	arf_t ctmp;

	mpfr_inits(re, im, NULL);
	arf_init(ctmp);
	arb_init(x);
	arb_init(y);
	acb_init(z);
	acb_init(s);

	mpfr_sub(re, e0, estar, ROUNDING);
	mpfr_mul(im, re, pi, ROUNDING);
	mpfr_div(im, im, sigma, ROUNDING);
	mpfr_mul(re, re, theta, ROUNDING);

	arf_set_mpfr(ctmp, re);
	arb_set_arf(x, ctmp);
	arf_set_mpfr(ctmp, im);
	arb_set_arf(y, ctmp);

	acb_set_arb_arb(z, x, y);
	acb_one(s);
	acb_hypgeom_expint(z, s, z, 500);
	acb_get_imag(x, z);
	arb_const_pi(y, 500);
	arb_div(x, x, y, 500);
	arb_set_d(y, 1.0);
	arb_sub(x, y, x, 500);

	arf_get_mpfr(res, arb_midref(x), ROUNDING);

	mpfr_mul(re, estar, theta, ROUNDING);
	mpfr_neg(re, re, ROUNDING);
	mpfr_exp(re, re, ROUNDING);
	mpfr_mul(res, res, re, ROUNDING);

	mpfr_clears(re, im, NULL);
	arf_clear(ctmp);
	arb_clear(x);
	arb_clear(y);
	acb_clear(z);
	acb_clear(s);
}

void int_delta2(mpfr_t res, mpfr_t theta, mpfr_t estar, mpfr_t e0)
{
	mpfr_t re0, re1, im;
	arf_t ctmp;
	arb_t x0, x1, y;
	acb_t u, v, s;

	mpfr_inits(re0, re1, im, NULL);
	arf_init(ctmp);
	arb_init(x0);
	arb_init(x1);
	arb_init(y);
	acb_init(u);
	acb_init(v);
	acb_init(s);

	mpfr_mul(re0, estar, theta, ROUNDING);
	mpfr_neg(re0, re0, ROUNDING);
	mpfr_sub(re1, e0, estar, ROUNDING);
	mpfr_mul(re1, re1, theta, ROUNDING);
	mpfr_mul(im, sigma, theta, ROUNDING);

	arf_set_mpfr(ctmp, re0);
	arb_set_arf(x0, ctmp);
	arf_set_mpfr(ctmp, re1);
	arb_set_arf(x1, ctmp);
	arf_set_mpfr(ctmp, im);
	arb_set_arf(y, ctmp);

	acb_set_arb_arb(u, x0, y);
	acb_set_arb_arb(v, x1, y);

	acb_exp(u, u, 500);
	acb_one(s);
	acb_hypgeom_expint(v, s, v, 500);
	acb_mul(v, u, v, 500);

	acb_get_imag(y, v);
	arb_neg(y, y);
	arf_get_mpfr(res, arb_midref(y), ROUNDING);

	mpfr_clears(re0, re1, im, NULL);
	arf_init(ctmp);
	arb_clear(x0);
	arb_clear(x1);
	arb_clear(y);
	acb_clear(u);
	acb_clear(v);
	acb_clear(s);
}

void int_delta3(mpfr_t res, mpfr_t theta, mpfr_t estar, mpfr_t e0)
{
	mpfr_t re0, re1, im;
	arf_t ctmp;
	arb_t x0, x1, y;
	acb_t u, v, s;

	mpfr_inits(re0, re1, im, NULL);
	arf_init(ctmp);
	arb_init(x0);
	arb_init(x1);
	arb_init(y);
	acb_init(u);
	acb_init(v);
	acb_init(s);

	mpfr_mul(re0, estar, theta, ROUNDING);
	mpfr_neg(re0, re0, ROUNDING);
	mpfr_sub(re1, e0, estar, ROUNDING);
	mpfr_mul(re1, re1, theta, ROUNDING);
	mpfr_mul(im, sigma, theta, ROUNDING);

	arf_set_mpfr(ctmp, re0);
	arb_set_arf(x0, ctmp);
	arf_set_mpfr(ctmp, re1);
	arb_set_arf(x1, ctmp);
	arf_set_mpfr(ctmp, im);
	arb_set_arf(y, ctmp);

	acb_set_arb_arb(u, x0, y);
	acb_set_arb_arb(v, x1, y);

	acb_exp(u, u, 500);
	acb_one(s);
	acb_hypgeom_expint(v, s, v, 500);
	acb_mul(v, u, v, 500);

	acb_get_real(y, v);
	arb_neg(y, y);
	arf_get_mpfr(res, arb_midref(y), ROUNDING);

	mpfr_clears(re0, re1, im, NULL);
	arf_init(ctmp);
	arb_clear(x0);
	arb_clear(x1);
	arb_clear(y);
	acb_clear(u);
	acb_clear(v);
	acb_clear(s);
}

void init_kernel(double s)
{
	static int init = 0;

	if(init == 0)
	{
		mpfr_inits(pi, sqrt2, sigma, NULL);

		mpfr_set_d(pi, 1.0, ROUNDING);
		mpfr_asin(pi, pi, ROUNDING);
		mpfr_mul_d(pi, pi, 2.0, ROUNDING);
		mpfr_set_d(sqrt2, 2.0, ROUNDING);
		mpfr_sqrt(sqrt2, sqrt2, ROUNDING);
	}

	mpfr_set_d(sigma, s, ROUNDING);
}
