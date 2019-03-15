#include <stdlib.h>
#include <stdio.h>
#include <mpfr.h>
#include "corefcts.h"
#include "utils.h"

#define CORR_EXP  0x01
#define CORR_COSH 0x02

static int type;
static int inorm = 1;

static mpfr_t pi;
static mpfr_t sqrt2;
static mpfr_t lambda, clambda;
static mpfr_t alpha, sigma;

static int tmax = 0;
static int mmax = 0;

static mpfr_t *w = 0;
static mpfr_t *g = 0;
static mpfr_t *r = 0;
static mpfr_t *f = 0;
static mpfr_t *wr = 0;
static mpfr_t *wf = 0;
static mpfr_t *ws = 0;

static void ludcmp(int n, mpfr_t *a, int *mutate)
{
	mpfr_t big, tmp;
	int row, ndx;
	mpfr_inits(big, tmp, NULL);

	for(int i = 0; i < n; i++)
	{
		mutate[i] = i;
	}

	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j <= i; j++)
		{
			for(int k = 0; k < j; k++)
			{
				mpfr_mul(tmp, a[j*n+k], a[k*n+i], ROUNDING);
				mpfr_sub(a[j*n+i], a[j*n+i], tmp, ROUNDING);
			}
		}

		mpfr_abs(big, a[i*n+i], ROUNDING);
		row = i;

		for(int j = i+1; j < n; j++)
		{
			for(int k = 0; k < i; k++)
			{
				mpfr_mul(tmp, a[j*n+k], a[k*n+i], ROUNDING);
				mpfr_sub(a[j*n+i], a[j*n+i], tmp, ROUNDING);
			}

			mpfr_abs(tmp, a[j*n+i], ROUNDING);
			if(mpfr_greater_p(tmp, big))
			{
				mpfr_set(big, tmp, ROUNDING);
				row = j;
			}
		}

		if(mpfr_cmp_d(big, 0.0) <= 0)
		{
			error("LU decomposition failed - matrix is singular");
		}

		if(row != i)
		{
			for(int k = 0; k < n; k++)
			{
				mpfr_set(tmp, a[row*n+k], ROUNDING);
				mpfr_set(a[row*n+k], a[i*n+k], ROUNDING);
				mpfr_set(a[i*n+k], tmp, ROUNDING);
			}

			ndx = mutate[row];
			mutate[row] = mutate[i];
			mutate[i] = ndx;
		}

		for(int k = i+1; k < n; k++)
		{
			mpfr_div(a[k*n+i], a[k*n+i], a[i*n+i], ROUNDING);
		}
	}

	mpfr_clears(big, tmp, NULL);
}

static void lubksb(int n, mpfr_t *a, mpfr_t *x, mpfr_t *b, int *mutate)
{
	mpfr_t tmp;
	mpfr_init(tmp);

	for(int i = 0; i < n; i++)
	{
		mpfr_set(x[i], b[mutate[i]], ROUNDING);
		for(int k = 0; k < i; k++)
		{
			mpfr_mul(tmp, a[i*n+k], x[k], ROUNDING);
			mpfr_sub(x[i], x[i], tmp, ROUNDING);
		}
	}

	for(int i = n-1; i >= 0; i--)
	{
		for(int k = i+1; k < n; k++)
		{
			mpfr_mul(tmp, a[i*n+k], x[k], ROUNDING);
			mpfr_sub(x[i], x[i], tmp, ROUNDING);
		}
		mpfr_div(x[i], x[i], a[i*n+i], ROUNDING);
	}

	mpfr_clear(tmp);
}

static void int_delta(mpfr_t res, mpfr_t theta, mpfr_t e0, mpfr_t estar)
{
	mpfr_t arg, norm, tmp, rep;
	mpfr_inits(arg, norm, tmp, rep, NULL);

	mpfr_mul(tmp, sigma, sigma, ROUNDING);
	mpfr_mul(rep, theta, tmp, ROUNDING);

	mpfr_add(tmp, rep, estar, ROUNDING);
	mpfr_add(tmp, tmp, estar, ROUNDING);
	mpfr_mul(tmp, tmp, theta, ROUNDING);
	mpfr_mul_d(tmp, tmp, 0.5, ROUNDING);
	mpfr_exp(norm, tmp, ROUNDING);
	mpfr_mul_d(norm, norm, 0.5, ROUNDING);

	mpfr_sub(arg, rep, e0, ROUNDING);
	mpfr_add(arg, arg, estar, ROUNDING);
	mpfr_mul(tmp, sqrt2, sigma, ROUNDING);
	mpfr_div(arg, arg, tmp, ROUNDING);
	mpfr_erf(tmp, arg, ROUNDING);
	mpfr_add_d(res, tmp, 1.0, ROUNDING);
	mpfr_mul(res, res, norm, ROUNDING);

	mpfr_mul(norm, sqrt2, sigma, ROUNDING);
	mpfr_div(norm, estar, norm, ROUNDING);
	mpfr_erf(norm, norm, ROUNDING);
	mpfr_add_d(norm, norm, 1.0, ROUNDING);
	mpfr_mul_d(norm, norm, 0.5, ROUNDING);
	mpfr_div(res, res, norm, ROUNDING);

	mpfr_clears(arg, norm, tmp, rep, NULL);
}

static void int_delta_sq(mpfr_t res, mpfr_t theta, mpfr_t e0, mpfr_t estar)
{
	mpfr_t arg, norm, tmp, rep;
	mpfr_inits(arg, norm, tmp, rep, NULL);

	mpfr_mul(tmp, sigma, sigma, ROUNDING);
	mpfr_mul(rep, theta, tmp, ROUNDING);

	mpfr_mul_d(tmp, estar, 4.0, ROUNDING);
	mpfr_add(tmp, tmp, rep, ROUNDING);
	mpfr_mul(tmp, tmp, theta, ROUNDING);
	mpfr_mul_d(tmp, tmp, 0.25, ROUNDING);
	mpfr_exp(norm, tmp, ROUNDING);
	mpfr_sqrt(tmp, pi, ROUNDING);
	mpfr_mul(tmp, tmp, sigma, ROUNDING);
	mpfr_mul_d(tmp, tmp, 4.0, ROUNDING);
	mpfr_div(norm, norm, tmp, ROUNDING);

	mpfr_sub(arg, estar, e0, ROUNDING);
	mpfr_mul_d(arg, arg, 2.0, ROUNDING);
	mpfr_add(arg, arg, rep, ROUNDING);
	mpfr_mul_d(tmp, sigma, 2.0, ROUNDING);
	mpfr_div(arg, arg, tmp, ROUNDING);
	mpfr_erf(tmp, arg, ROUNDING);
	mpfr_add_d(res, tmp, 1.0, ROUNDING);
	mpfr_mul(res, res, norm, ROUNDING);

	mpfr_mul(norm, sqrt2, sigma, ROUNDING);
	mpfr_div(norm, estar, norm, ROUNDING);
	mpfr_erf(norm, norm, ROUNDING);
	mpfr_add_d(norm, norm, 1.0, ROUNDING);
	mpfr_mul_d(norm, norm, 0.5, ROUNDING);
	mpfr_sqr(norm, norm, ROUNDING);
	mpfr_div(res, res, norm, ROUNDING);

	mpfr_clears(arg, norm, tmp, rep, NULL);
}

void delta(mpfr_t res, mpfr_t estar, mpfr_t e)
{
	mpfr_t arg, norm;
	mpfr_inits(arg, norm, NULL);

	mpfr_mul_d(norm, pi, 2.0, ROUNDING);
	mpfr_sqrt(norm, norm, ROUNDING);
	mpfr_mul(norm, norm, sigma, ROUNDING);

	mpfr_sub(arg, e, estar, ROUNDING);
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

	mpfr_clears(arg, norm, NULL);
}

void bg_method(mpfr_t estar, mpfr_t *cov)
{
	int *mutate, n;
	mpfr_t inm, tmp;

	type = CORR_EXP;
	n = tmax-1;

	mpfr_inits(inm, tmp, NULL);
	mutate = malloc(n*sizeof(int));
	if(mutate == NULL)
	{
		error("Unable to allocate auxiliary array");
	}

	for(int i = 0; i < n; i++)
	{
		mpfr_set_d(tmp, (double)(i+1), ROUNDING);
		mpfr_d_div(r[i], 1.0, tmp, ROUNDING);

		for(int j = i; j < n; j++)
		{
			mpfr_set_d(inm, (double)(i+j+2), ROUNDING);
			mpfr_mul(tmp, estar, estar, ROUNDING);
			mpfr_div(w[j*n+i], tmp, inm, ROUNDING);

			mpfr_mul(tmp, inm, inm, ROUNDING);
			mpfr_div(tmp, estar, tmp, ROUNDING);
			mpfr_mul_d(tmp, tmp, -2.0, ROUNDING);
			mpfr_add(w[j*n+i], w[j*n+i], tmp, ROUNDING);

			mpfr_mul(tmp, inm, inm, ROUNDING);
			mpfr_mul(tmp, tmp, inm, ROUNDING);
			mpfr_d_div(tmp, 2.0, tmp, ROUNDING);
			mpfr_add(w[j*n+i], w[j*n+i], tmp, ROUNDING);

			mpfr_mul(w[j*n+i], w[j*n+i], clambda, ROUNDING);
			mpfr_set(w[i*n+j], w[j*n+i], ROUNDING);
		}

		mpfr_mul(tmp, cov[i+1], lambda, ROUNDING);
		mpfr_add(w[i*n+i], w[i*n+i], tmp, ROUNDING);
	}

	ludcmp(n, w, mutate);
	lubksb(n, w, wr, r, mutate);

	mpfr_set_zero(inm, 1);
	for(int i = 0; i < n; i++)
	{
		mpfr_mul(tmp, r[i], wr[i], ROUNDING);
		mpfr_add(inm, inm, tmp, ROUNDING);
	}

	for(int i = 0; i < n; i++)
	{
		mpfr_div(g[i], wr[i], inm, ROUNDING);
	}

	mpfr_clears(inm, tmp, NULL);
	free(mutate);
}

void rm_method(mpfr_t e0, mpfr_t estar, mpfr_t *cov)
{
	int *mutate, n;
	mpfr_t tmp, x0, x1;

	type = CORR_EXP;
	n = tmax-1;

	mpfr_inits(tmp, x0, x1, NULL);
	mutate = malloc(n*sizeof(int));
	if(mutate == NULL)
	{
		error("Unable to allocate auxiliary array");
	}

	for(int i = 0; i < n; i++)
	{
		mpfr_sub_d(tmp, alpha, (double)(i+1), ROUNDING);
		int_delta(f[i], tmp, e0, estar);
		mpfr_mul(f[i], f[i], clambda, ROUNDING);

		if(inorm)
		{
			mpfr_set_d(r[i], 1.0, ROUNDING);
			mpfr_div_d(r[i], r[i], (double)(i+1), ROUNDING);
		}

		for(int j = i; j < n; j++)
		{
			mpfr_d_sub(tmp, (double)(i+j+2), alpha, ROUNDING);
			mpfr_mul(x0, e0, tmp, ROUNDING);
			mpfr_neg(x0, x0, ROUNDING);
			mpfr_exp(x0, x0, ROUNDING);
			mpfr_div(x0, x0, tmp, ROUNDING);

			mpfr_mul(x0, x0, clambda, ROUNDING);
			mpfr_set(w[j*n+i], x0, ROUNDING);
			mpfr_set(w[i*n+j], x0, ROUNDING);
			mpfr_set(ws[j*n+i], x0, ROUNDING);
			mpfr_set(ws[i*n+j], x0, ROUNDING);
		}

		mpfr_mul(tmp, cov[i+1], lambda, ROUNDING);
		mpfr_add(w[i*n+i], w[i*n+i], tmp, ROUNDING);

	}

	ludcmp(n, w, mutate);
	lubksb(n, w, wf, f, mutate);

	if(inorm)
	{
		lubksb(n, w, wr, r, mutate);

		mpfr_set_d(x0, 1.0, ROUNDING);
		mpfr_set_d(x1, 0.0, ROUNDING);

		for(int i = 0; i < n; i++)
		{
			mpfr_mul(tmp, r[i], wf[i], ROUNDING);
			mpfr_sub(x0, x0, tmp, ROUNDING);
			mpfr_mul(tmp, r[i], wr[i], ROUNDING);
			mpfr_add(x1, x1, tmp, ROUNDING);
		}

		mpfr_div(x0, x0, x1, ROUNDING);
	}

	for(int i = 0; i < n; i++)
	{
		if(inorm)
		{
			mpfr_mul(g[i], wr[i], x0, ROUNDING);
			mpfr_add(g[i], g[i], wf[i], ROUNDING);
		}
		else
		{
			mpfr_set(g[i], wf[i], ROUNDING);
		}
	}

	mpfr_clears(tmp, x0, x1, NULL);
	free(mutate);
}

void rm_method_cosh(mpfr_t e0, mpfr_t estar, mpfr_t *cov)
{
	int *mutate, n;
	mpfr_t tmp, x0, x1;

	type = CORR_COSH;
	n = tmax/2;

	mpfr_inits(tmp, x0, x1, NULL);
	mutate = malloc(n*sizeof(int));
	if(mutate == NULL)
	{
		error("Unable to allocate auxiliary array");
	}

	for(int i = 0; i < n; i++)
	{
		mpfr_sub_d(tmp, alpha, (double)(i+1), ROUNDING);
		int_delta(x0, tmp, e0, estar);
		mpfr_sub_d(tmp, alpha, (double)(tmax-i-1), ROUNDING);
		int_delta(x1, tmp, e0, estar);
		mpfr_add(f[i], x0, x1, ROUNDING);
		mpfr_mul(f[i], f[i], clambda, ROUNDING);

		if(inorm)
		{
			mpfr_set_d(x0, (double)(i+1), ROUNDING);
			mpfr_d_div(x0, 1.0, x0, ROUNDING);
			mpfr_set_d(x1, (double)(tmax-i-1), ROUNDING);
			mpfr_d_div(x1, 1.0, x1, ROUNDING);
			mpfr_add(r[i], x0, x1, ROUNDING);
		}

		for(int j = i; j < n; j++)
		{
			mpfr_d_sub(tmp, (double)(i+j+2), alpha, ROUNDING);
			mpfr_mul(x0, e0, tmp, ROUNDING);
			mpfr_neg(x0, x0, ROUNDING);
			mpfr_exp(x0, x0, ROUNDING);
			mpfr_div(x0, x0, tmp, ROUNDING);
			mpfr_set(x1, x0, ROUNDING);

			mpfr_d_sub(tmp, (double)(tmax+i-j), alpha, ROUNDING);
			mpfr_mul(x0, e0, tmp, ROUNDING);
			mpfr_neg(x0, x0, ROUNDING);
			mpfr_exp(x0, x0, ROUNDING);
			mpfr_div(x0, x0, tmp, ROUNDING);
			mpfr_add(x1, x1, x0, ROUNDING);

			mpfr_d_sub(tmp, (double)(tmax-i+j), alpha, ROUNDING);
			mpfr_mul(x0, e0, tmp, ROUNDING);
			mpfr_neg(x0, x0, ROUNDING);
			mpfr_exp(x0, x0, ROUNDING);
			mpfr_div(x0, x0, tmp, ROUNDING);
			mpfr_add(x1, x1, x0, ROUNDING);

			mpfr_d_sub(tmp, (double)(2*tmax-i-j-2), alpha, ROUNDING);
			mpfr_mul(x0, e0, tmp, ROUNDING);
			mpfr_neg(x0, x0, ROUNDING);
			mpfr_exp(x0, x0, ROUNDING);
			mpfr_div(x0, x0, tmp, ROUNDING);
			mpfr_add(x1, x1, x0, ROUNDING);

			mpfr_mul(x1, x1, clambda, ROUNDING);
			mpfr_set(w[j*n+i], x1, ROUNDING);
			mpfr_set(w[i*n+j], x1, ROUNDING);
			mpfr_set(ws[j*n+i], x1, ROUNDING);
			mpfr_set(ws[i*n+j], x1, ROUNDING);
		}

		mpfr_mul(tmp, cov[i+1], lambda, ROUNDING);
		mpfr_add(w[i*n+i], w[i*n+i], tmp, ROUNDING);
	}

	ludcmp(n, w, mutate);
	lubksb(n, w, wf, f, mutate);

	if(inorm)
	{
		lubksb(n, w, wr, r, mutate);

		mpfr_set_d(x0, 1.0, ROUNDING);
		mpfr_set_d(x1, 0.0, ROUNDING);

		for(int i = 0; i < n; i++)
		{
			mpfr_mul(tmp, r[i], wf[i], ROUNDING);
			mpfr_sub(x0, x0, tmp, ROUNDING);
			mpfr_mul(tmp, r[i], wr[i], ROUNDING);
			mpfr_add(x1, x1, tmp, ROUNDING);
		}

		mpfr_div(x0, x0, x1, ROUNDING);
	}

	for(int i = 0; i < n; i++)
	{
		if(inorm)
		{
			mpfr_mul(g[i], wr[i], x0, ROUNDING);
			mpfr_add(g[i], g[i], wf[i], ROUNDING);
		}
		else
		{
 			mpfr_set(g[i], wf[i], ROUNDING);
		}
	}

	mpfr_clears(tmp, x0, x1, NULL);
	free(mutate);
}

void transform(mpfr_t e0, mpfr_t estar, mpfr_t *corr, mpfr_t *cov, double *rho, double *ag, double *bg)
{
	mpfr_t term, sum, bsum, asum;
	int n;

	mpfr_inits(term, sum, bsum, asum, NULL);
	mpfr_set_zero(sum, 1);
	mpfr_set_zero(bsum, 1);
	mpfr_set_zero(asum, 1);

	if(type == CORR_EXP)
	{
		n = tmax-1;
	}
	else
	{
		n = tmax/2;
	}

	for(int i = 0; i < n; i++)
	{
		mpfr_mul(term, g[i], corr[i+1], ROUNDING);
		mpfr_add(sum, sum, term, ROUNDING);

		mpfr_mul(term, g[i], f[i], ROUNDING);
		mpfr_mul_d(term, term, -2.0, ROUNDING);
		mpfr_add(asum, asum, term, ROUNDING);

		for(int j = 0; j < n; j++)
		{
			mpfr_mul(term, g[i], ws[i*n+j], ROUNDING);
			mpfr_mul(term, term, g[j], ROUNDING);
			mpfr_add(asum, asum, term, ROUNDING);
		}

		mpfr_mul(term, g[i], cov[i+1], ROUNDING);
		mpfr_mul(term, term, g[i], ROUNDING);
		mpfr_mul(term, term, lambda, ROUNDING);
		mpfr_add(bsum, bsum, term, ROUNDING);
	}

	int_delta_sq(term, alpha, e0, estar);
	mpfr_mul(term, term, clambda, ROUNDING);
	mpfr_add(asum, asum, term, ROUNDING);

	*rho = mpfr_get_d(sum, ROUNDING);
	*ag = mpfr_get_d(asum, ROUNDING);
	*bg = mpfr_get_d(bsum, ROUNDING);

	mpfr_clears(term, sum, asum, bsum, NULL);
}

void deltabar(mpfr_t sum, mpfr_t e)
{
	mpfr_t term;
	int n;

	mpfr_init(term);
	mpfr_set_zero(sum, 1);

	if(type == CORR_EXP)
	{
		n = tmax-1;
	}
	else
	{
		n = tmax/2;
	}

	for(int i = 0; i < n; i++)
	{
		mpfr_mul_d(term, e, (double)(-i-1), ROUNDING);
		mpfr_exp(term, term, ROUNDING);
		mpfr_mul(term, term, g[i], ROUNDING);
		mpfr_add(sum, sum, term, ROUNDING);

		if(type == CORR_COSH)
		{
			mpfr_mul_d(term, e, (double)(-tmax+i+1), ROUNDING);
			mpfr_exp(term, term, ROUNDING);
			mpfr_mul(term, term, g[i], ROUNDING);
			mpfr_add(sum, sum, term, ROUNDING);
		}
	}

	mpfr_clear(term);
}

void set_params(double s, double a, double l)
{
	static int init = 0;

	if(init == 0)
	{
		mpfr_inits(pi, sqrt2, sigma, alpha, lambda, clambda, NULL);

		mpfr_set_d(pi, 1.0, ROUNDING);
		mpfr_asin(pi, pi, ROUNDING);
		mpfr_mul_d(pi, pi, 2.0, ROUNDING);
		mpfr_set_d(sqrt2, 2.0, ROUNDING);
		mpfr_sqrt(sqrt2, sqrt2, ROUNDING);

		type = CORR_EXP;
		init = 1;
	}

	mpfr_set_d(sigma, s, ROUNDING);
	mpfr_set_d(alpha, a, ROUNDING);
	mpfr_set_d(lambda, l, ROUNDING);
	mpfr_set_d(clambda, 1.0-l, ROUNDING);
}

void set_tmax(int tmx)
{
	int nmx;

	if(tmx < mmax)
	{
		tmax = tmx;
		return;
	}

	if(g)
	{
		nmx = (2*mmax+5)*mmax;
		for(int i = 0; i < nmx; i++)
		{
			mpfr_clear(g[i]);
		}
		free(g);
	}

	nmx = (2*tmx+5)*tmx;
	g = malloc(nmx*sizeof(mpfr_t));
	if(g == NULL)
	{
		error("Failed to allocate MPFR variables");
	}

	r = g+tmx;
	wr = r+tmx;
	f = wr+tmx;
	wf = f+tmx;
	w = wf+tmx;
	ws = w+tmx*tmx;

	for(int i = 0; i < nmx; i++)
	{
		mpfr_init(g[i]);
	}

	tmax = tmx;
	mmax = tmx;
}
