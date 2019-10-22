#include <stdlib.h>
#include <stdio.h>
#include <mpfr.h>
#include "corefcts.h"
#include "kernels.h"
#include "utils.h"

#define CORR_EXP  0x01
#define CORR_COSH 0x02

static int type;
static int kernel;
static int inorm = 1;

static mpfr_t lambda;
static mpfr_t clambda;

static int tmin = 1;
static int tmax = 0;
static int tlen = 0;
static int mlen = 0;

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

void delta(mpfr_t res, mpfr_t estar, mpfr_t omega)
{
	switch(kernel)
	{
		case 0:
			delta0(res, estar, omega);
			break;
		case 1:
			delta1(res, estar, omega);
			break;
		case 2:
			delta2(res, estar, omega);
			break;
		case 3:
			delta3(res, estar, omega);
			break;
	}
}

void int_delta_sq(mpfr_t res, mpfr_t estar, mpfr_t e0)
{
	switch(kernel)
	{
		case 0:
			int_delta0_sq(res, estar, e0);
			break;
		case 1:
			int_delta1_sq(res, estar, e0);
			break;
		case 2:
			int_delta2_sq(res, estar, e0);
			break;
		case 3:
			int_delta3_sq(res, estar, e0);
			break;
	}
}

void int_delta(mpfr_t res, mpfr_t theta, mpfr_t estar, mpfr_t e0)
{
	switch(kernel)
	{
		case 0:
			int_delta0(res, theta, estar, e0);
			break;
		case 1:
			int_delta1(res, theta, estar, e0);
			break;
		case 2:
			int_delta2(res, theta, estar, e0);
			break;
		case 3:
			int_delta3(res, theta, estar, e0);
			break;
	}
}

void bg_method(mpfr_t estar, mpfr_t *cov)
{
	int *mutate, n;
	mpfr_t inm, tmp;

	type = CORR_EXP;
	n = tlen;

	mpfr_inits(inm, tmp, NULL);
	mutate = malloc(n*sizeof(int));
	if(mutate == NULL)
	{
		error("Unable to allocate auxiliary array");
	}

	for(int i = 0; i < n; i++)
	{
		mpfr_set_d(tmp, (double)(i+tmin), ROUNDING);
		mpfr_d_div(r[i], 1.0, tmp, ROUNDING);

		for(int j = i; j < n; j++)
		{
			mpfr_set_d(inm, (double)(i+j+2*tmin), ROUNDING);
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

		mpfr_mul(tmp, cov[i+tmin], lambda, ROUNDING);
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
	n = tlen;

	mpfr_inits(tmp, x0, x1, NULL);
	mutate = malloc(n*sizeof(int));
	if(mutate == NULL)
	{
		error("Unable to allocate auxiliary array");
	}

	for(int i = 0; i < n; i++)
	{
		mpfr_set_d(tmp, (double)(i+tmin), ROUNDING);
		int_delta(f[i], tmp, estar, e0);
		mpfr_mul(f[i], f[i], clambda, ROUNDING);

		if(inorm)
		{
			mpfr_set_d(r[i], 1.0, ROUNDING);
			mpfr_div_d(r[i], r[i], (double)(i+tmin), ROUNDING);
		}

		for(int j = i; j < n; j++)
		{
			mpfr_set_d(tmp, (double)(i+j+2*tmin), ROUNDING);
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

		mpfr_mul(tmp, cov[i+tmin], lambda, ROUNDING);
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
	n = tlen;

	mpfr_inits(tmp, x0, x1, NULL);
	mutate = malloc(n*sizeof(int));
	if(mutate == NULL)
	{
		error("Unable to allocate auxiliary array");
	}

	for(int i = 0; i < n; i++)
	{
		mpfr_set_d(tmp, (double)(i+tmin), ROUNDING);
		int_delta(x0, tmp, estar, e0);
		mpfr_set_d(tmp, (double)(tmax-i-tmin), ROUNDING);
		int_delta(x1, tmp, estar, e0);
		mpfr_add(f[i], x0, x1, ROUNDING);
		mpfr_mul(f[i], f[i], clambda, ROUNDING);

		if(inorm)
		{
			mpfr_set_d(x0, (double)(i+tmin), ROUNDING);
			mpfr_d_div(x0, 1.0, x0, ROUNDING);
			mpfr_set_d(x1, (double)(tmax-i-tmin), ROUNDING);
			mpfr_d_div(x1, 1.0, x1, ROUNDING);
			mpfr_add(r[i], x0, x1, ROUNDING);
		}

		for(int j = i; j < n; j++)
		{
			mpfr_set_d(tmp, (double)(i+j+2*tmin), ROUNDING);
			mpfr_mul(x0, e0, tmp, ROUNDING);
			mpfr_neg(x0, x0, ROUNDING);
			mpfr_exp(x0, x0, ROUNDING);
			mpfr_div(x0, x0, tmp, ROUNDING);
			mpfr_set(x1, x0, ROUNDING);

			mpfr_set_d(tmp, (double)(tmax+i-j), ROUNDING);
			mpfr_mul(x0, e0, tmp, ROUNDING);
			mpfr_neg(x0, x0, ROUNDING);
			mpfr_exp(x0, x0, ROUNDING);
			mpfr_div(x0, x0, tmp, ROUNDING);
			mpfr_add(x1, x1, x0, ROUNDING);

			mpfr_set_d(tmp, (double)(tmax-i+j), ROUNDING);
			mpfr_mul(x0, e0, tmp, ROUNDING);
			mpfr_neg(x0, x0, ROUNDING);
			mpfr_exp(x0, x0, ROUNDING);
			mpfr_div(x0, x0, tmp, ROUNDING);
			mpfr_add(x1, x1, x0, ROUNDING);

			mpfr_set_d(tmp, (double)(2*tmax-i-j-2*tmin), ROUNDING);
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

		mpfr_mul(tmp, cov[i+tmin], lambda, ROUNDING);
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
	n = tlen;

	for(int i = 0; i < n; i++)
	{
		mpfr_mul(term, g[i], corr[i+tmin], ROUNDING);
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

		mpfr_mul(term, g[i], cov[i+tmin], ROUNDING);
		mpfr_mul(term, term, g[i], ROUNDING);
		mpfr_mul(term, term, lambda, ROUNDING);
		mpfr_add(bsum, bsum, term, ROUNDING);
	}

	int_delta_sq(term, estar, e0);
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
	n = tlen;

	for(int i = 0; i < n; i++)
	{
		mpfr_mul_d(term, e, (double)(-i-tmin), ROUNDING);
		mpfr_exp(term, term, ROUNDING);
		mpfr_mul(term, term, g[i], ROUNDING);
		mpfr_add(sum, sum, term, ROUNDING);

		if(type == CORR_COSH)
		{
			mpfr_mul_d(term, e, (double)(-tmax+i+tmin), ROUNDING);
			mpfr_exp(term, term, ROUNDING);
			mpfr_mul(term, term, g[i], ROUNDING);
			mpfr_add(sum, sum, term, ROUNDING);
		}
	}

	mpfr_clear(term);
}

void set_params(double s, double l, int k)
{
	static int init = 0;

	if(init == 0)
	{
		mpfr_inits(lambda, clambda, NULL);
		type = CORR_EXP;
		init = 1;
	}

	if(k < 0 || k > 3)
	{
		error("Invalid smearing kernel");
	}
	else
	{
		kernel = k;
		inorm = (k == 0);
	}

	mpfr_set_d(lambda, l, ROUNDING);
	mpfr_set_d(clambda, 1.0-l, ROUNDING);

	init_kernel(s);
}

void set_time_parms(int ta, int tb, int tc)
{
	int nmx;
	int len;

	if((ta < 1) || (ta >= tb) || (tc <= tb))
	{
		error("Invalid parameters for Euclidean time range");
	}

	len = (tb - ta + 1);

	if(len < mlen)
	{
		tmin = ta;
		tmax = tc;
		tlen = len;
		return;
	}

	if(g)
	{
		nmx = (2*mlen+5)*mlen;
		for(int i = 0; i < nmx; i++)
		{
			mpfr_clear(g[i]);
		}
		free(g);
	}

	nmx = (2*len+5)*len;
	g = malloc(nmx*sizeof(mpfr_t));
	if(g == NULL)
	{
		error("Failed to allocate MPFR variables");
	}

	r = g+len;
	wr = r+len;
	f = wr+len;
	wf = f+len;
	w = wf+len;
	ws = w+len*len;

	for(int i = 0; i < nmx; i++)
	{
		mpfr_init(g[i]);
	}

	tmin = ta;
	tmax = tc;
	tlen = len;
	mlen = len;
}
