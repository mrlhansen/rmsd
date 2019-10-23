#include <sys/stat.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpfr.h>
#include <math.h>
#include "corefcts.h"
#include "utils.h"

double sigma = 0.1;
double lambda = 0.05;
int nsteps = 100;
double ei = 0.0;
double ef = 1.0;
int tmx = 0;
int nms = 1;
int prec = 128;
double emin = 0.0;
double escan = 0.0;
int sym = 0;
int L = 32;
int kid = -1;

mpfr_t *corr, *cov;
mpfr_t *expected;
mpfr_t *mcorr;

char path[128];
double de;

FILE *fopen_path(char *name)
{
	char fn[128];
	FILE *fp;

	sprintf(fn, "%s/%s", path, name);
	fp = fopen(fn, "w");
	if(fp == NULL)
	{
		error("Failed to open file for writing");
	}

	return fp;
}

void write_fakecorr()
{
	FILE *fp;

	fp = fopen_path("exact.txt");
	fprintf(fp, "%d %d\n", tmx, nsteps);

	for(int i = 0; i < tmx; i++)
	{
		mpfr_out_str(fp, 10, 0, corr[i], ROUNDING);
		fprintf(fp, "\n");
	}

	for(int i = 0; i < nsteps; i++)
	{
		mpfr_out_str(fp, 10, 0, expected[i], ROUNDING);
		fprintf(fp, "\n");
	}

	fclose(fp);
}

void read_fakecorr(const char *fn)
{
	FILE *fp;
	char tmp[16];
	int a, b;

	fp = fopen(fn, "r");

	if(fp == NULL)
	{
		error("Unable to open file for reading");
	}

	fgets(tmp, 16, fp);
	sscanf(tmp, "%d %d", &a, &b);

	if(a != tmx || b != nsteps)
	{
		error("Temporal extent or nsteps does not match");
	}

	for(int i = 0; i < tmx; i++)
	{
		mpfr_inp_str(corr[i], fp, 10, ROUNDING);
		fseek(fp, 1, SEEK_CUR);
	}

	for(int i = 0; i < nsteps; i++)
	{
		mpfr_inp_str(expected[i], fp, 10, ROUNDING);
		fseek(fp, 1, SEEK_CUR);
	}

	fclose(fp);
}

void add_peak(mpfr_t mass, mpfr_t w)
{
	mpfr_t estar, res;
	mpfr_inits(estar, res, NULL);

	for(int n = 0; n < nsteps; n++)
	{
		mpfr_set_d(estar, ei+n*de, ROUNDING);
		delta(res, estar, mass);
		mpfr_mul(res, res, w, ROUNDING);
		mpfr_add(expected[n], expected[n], res, ROUNDING);
	}

	mpfr_clears(estar, res, NULL);
}

void fakecorr()
{
	double dmpi, dmk, dmphi, dcutoff;
	mpfr_t mpi, mk, mphi, cutoff;
	mpfr_t E1, E2, E3, Esum;
	mpfr_t mf, mom, w, expE;
	mpfr_t f1, f2, cg, cl;
	int psq, Lmax;

	mpfr_inits(E1, E2, E3, Esum, mf, mom, w, expE, f1, f2, NULL);
	mpfr_inits(mpi, mk, mphi, cutoff, cg, cl, NULL);

	dmpi = 0.066;
	dmk = 3.55*dmpi;
	dmphi = 7.30*dmpi;
	dcutoff = 1.25*ef;
	Lmax = 1.0 + (L/(2.0*M_PI))*sqrt(dcutoff*dcutoff-dmpi*dmpi);
	printf("Generating correlator with Lmax = %d\n", Lmax);

	mpfr_set_d(cg, 1.0, ROUNDING);
	mpfr_set_d(cl, 10.0*sqrt(8.0), ROUNDING);
	mpfr_set_d(mpi, dmpi, ROUNDING);
	mpfr_set_d(mk, dmk, ROUNDING);
	mpfr_set_d(mphi, dmphi, ROUNDING);
	mpfr_set_d(cutoff, dcutoff, ROUNDING);

	for(int n = 0; n < nsteps; n++)
	{
		mpfr_set_zero(expected[n], 1);
	}

	mpfr_set_d(mf, 1.0, ROUNDING);
	mpfr_asin(mf, mf, ROUNDING);
	mpfr_mul_d(mf, mf, 4.0, ROUNDING);
	mpfr_div_d(mf, mf, (double)L, ROUNDING);
	mpfr_mul(mf, mf, mf, ROUNDING);

	mpfr_mul(f1, cg, cg, ROUNDING);
	mpfr_mul(f1, f1, mphi, ROUNDING);
	mpfr_mul(f1, f1, mphi, ROUNDING);
	mpfr_div_d(f1, f1, 2.0, ROUNDING);
	mpfr_div(f1, f1, mpi, ROUNDING);
	mpfr_div(f1, f1, mpi, ROUNDING);
	mpfr_div(f1, f1, mpi, ROUNDING);
	mpfr_div_d(f1, f1, (double)L, ROUNDING);
	mpfr_div_d(f1, f1, (double)L, ROUNDING);
	mpfr_div_d(f1, f1, (double)L, ROUNDING);

	mpfr_mul(f2, cl, cl, ROUNDING);
	mpfr_div_d(f2, f2, 48.0, ROUNDING);
	mpfr_div(f2, f2, mpi, ROUNDING);
	mpfr_div(f2, f2, mpi, ROUNDING);
	mpfr_div(f2, f2, mpi, ROUNDING);
	mpfr_div_d(f2, f2, (double)L, ROUNDING);
	mpfr_div_d(f2, f2, (double)L, ROUNDING);
	mpfr_div_d(f2, f2, (double)L, ROUNDING);
	mpfr_div_d(f2, f2, (double)L, ROUNDING);
	mpfr_div_d(f2, f2, (double)L, ROUNDING);
	mpfr_div_d(f2, f2, (double)L, ROUNDING);

	for(int i = 0; i < tmx; i++)
	{
		mpfr_set_zero(corr[i], 1);
	}

	for(int p1 = 1-Lmax; p1 < Lmax; p1++)
	for(int p2 = 1-Lmax; p2 < Lmax; p2++)
	for(int p3 = 1-Lmax; p3 < Lmax; p3++)
	{
		psq = p1*p1 + p2*p2 + p3*p3;

		mpfr_mul_d(mom, mf, (double)psq, ROUNDING);
		mpfr_mul(E1, mk, mk, ROUNDING);
		mpfr_add(mom, mom, E1, ROUNDING);
		mpfr_sqrt(E1, mom, ROUNDING);
		mpfr_mul_d(E1, E1, 2.0, ROUNDING);

		if(mpfr_less_p(E1, cutoff))
		{
			mpfr_mul(mom, E1, E1, ROUNDING);
			mpfr_div(w, f1, mom, ROUNDING);
			add_peak(E1, w);

			for(int i = 0; i < tmx; i++)
			{
				mpfr_mul_d(expE, E1, (double)(-i), ROUNDING);
				mpfr_exp(expE, expE, ROUNDING);
				mpfr_mul(expE, expE, w, ROUNDING);
				mpfr_add(corr[i], corr[i], expE, ROUNDING);
				if(sym)
				{
					mpfr_mul_d(expE, E1, (double)(-tmx+i), ROUNDING);
					mpfr_exp(expE, expE, ROUNDING);
					mpfr_mul(expE, expE, w, ROUNDING);
					mpfr_add(corr[i], corr[i], expE, ROUNDING);
				}
			}
		}

		for(int q1 = 1-Lmax; q1 < Lmax; q1++)
		for(int q2 = 1-Lmax; q2 < Lmax; q2++)
		for(int q3 = 1-Lmax; q3 < Lmax; q3++)
		{
			psq = p1*p1 + p2*p2 + p3*p3;
			mpfr_mul_d(mom, mf, (double)psq, ROUNDING);
			mpfr_mul(E1, mpi, mpi, ROUNDING);
			mpfr_add(mom, mom, E1, ROUNDING);
			mpfr_sqrt(E1, mom, ROUNDING);

			psq = q1*q1 + q2*q2 + q3*q3;
			mpfr_mul_d(mom, mf, (double)psq, ROUNDING);
			mpfr_mul(E2, mpi, mpi, ROUNDING);
			mpfr_add(mom, mom, E2, ROUNDING);
			mpfr_sqrt(E2, mom, ROUNDING);

			psq = (p1+q1)*(p1+q1) + (p2+q2)*(p2+q2) + (p3+q3)*(p3+q3);
			mpfr_mul_d(mom, mf, (double)psq, ROUNDING);
			mpfr_mul(E3, mpi, mpi, ROUNDING);
			mpfr_add(mom, mom, E3, ROUNDING);
			mpfr_sqrt(E3, mom, ROUNDING);

			mpfr_add(Esum, E1, E2, ROUNDING);
			mpfr_add(Esum, Esum, E3, ROUNDING);

			if(mpfr_greater_p(Esum, cutoff))
			{
				continue;
			}

			mpfr_mul(mom, E1, E2, ROUNDING);
			mpfr_mul(mom, mom, E3, ROUNDING);
			mpfr_div(w, f2, mom, ROUNDING);
			add_peak(Esum, w);

			for(int i = 0; i < tmx; i++)
			{
				mpfr_mul_d(expE, Esum, (double)(-i), ROUNDING);
				mpfr_exp(expE, expE, ROUNDING);
				mpfr_mul(expE, expE, w, ROUNDING);
				mpfr_add(corr[i], corr[i], expE, ROUNDING);
				if(sym)
				{
					mpfr_mul_d(expE, Esum, (double)(-tmx+i), ROUNDING);
					mpfr_exp(expE, expE, ROUNDING);
					mpfr_mul(expE, expE, w, ROUNDING);
					mpfr_add(corr[i], corr[i], expE, ROUNDING);
				}
			}
		}
	}

	mpfr_clears(E1, E2, E3, Esum, mf, mom, w, expE, f1, f2, NULL);
	mpfr_clears(mpi, mk, mphi, cutoff, cg, cl, NULL);
}

void fakecov()
{
	double rsigma, *rd, res, dres;
	mpfr_t tmp;
	FILE *fp;

	fp = fopen_path("correlator.txt");
	rd = malloc(nms*sizeof(double));
	if(rd == NULL)
	{
		error("Failed to allocate auxiliary array");
	}

	rsigma = 0.25/sqrt((double)nms);
	mpfr_init(tmp);

	for(int i = 0; i < tmx; i++)
	{
		mpfr_set_zero(cov[i], 1);

		if(nms > 1)
		{
			randg(rd, nms);
		}
		else
		{
			continue;
		}

		for(int j = 0; j < nms; j++)
		{
			mpfr_set_d(tmp, rsigma*rd[j], ROUNDING);
			mpfr_mul(tmp, tmp, corr[i], ROUNDING);
			mpfr_add(mcorr[j*tmx+i], corr[i], tmp, ROUNDING);
			mpfr_mul(tmp, tmp, tmp, ROUNDING);
			mpfr_add(cov[i], cov[i], tmp, ROUNDING);
		}

		mpfr_div_d(cov[i], cov[i], (double)(nms-1), ROUNDING);
	}

	for(int i = 0; i < tmx; i++)
	{
		res = mpfr_get_d(corr[i], ROUNDING);
		dres = mpfr_get_d(cov[i], ROUNDING);
		dres = sqrt(dres);
		fprintf(fp, "%3d %1.8e %1.8e\n", i, res, dres);
	}

	free(rd);
	mpfr_clear(tmp);
	fclose(fp);
}

void full_sample()
{
	mpfr_t tmp;
	mpfr_init(tmp);

	if(nms <= 1)
	{
		return;
	}

	for(int i = 0; i < tmx; i++)
	{
		mpfr_set_zero(corr[i], 1);
		for(int j = 0; j < nms; j++)
		{
			mpfr_add(corr[i], corr[i], mcorr[j*tmx+i], ROUNDING);
		}
		mpfr_div_d(corr[i], corr[i], (double)nms, ROUNDING);

		mpfr_set_zero(cov[i], 1);
		for(int j = 0; j < nms; j++)
		{
			mpfr_sub(tmp, mcorr[j*tmx+i], corr[i], ROUNDING);
			mpfr_mul(tmp, tmp, tmp, ROUNDING);
			mpfr_add(cov[i], cov[i], tmp, ROUNDING);
		}
		mpfr_div_d(cov[i], cov[i], (double)(nms-1), ROUNDING);
		mpfr_mul(tmp, corr[0], corr[0], ROUNDING);
		mpfr_div(cov[i], cov[i], tmp, ROUNDING);
	}

	mpfr_clear(tmp);
}

void bootstrap_sample()
{
	mpfr_t tmp;
	int *rd;

	if(nms <= 1)
	{
		return;
	}

	mpfr_init(tmp);
	rd = malloc(sizeof(int)*nms);
	randi(rd, nms, 0, nms);

	for(int i = 0; i < tmx; i++)
	{
		mpfr_set_zero(corr[i], 1);
		for(int j = 0; j < nms; j++)
		{
			mpfr_add(corr[i], corr[i], mcorr[rd[j]*tmx+i], ROUNDING);
		}
		mpfr_div_d(corr[i], corr[i], (double)nms, ROUNDING);

		mpfr_set_zero(cov[i], 1);
		for(int j = 0; j < nms; j++)
		{
			mpfr_sub(tmp, mcorr[rd[j]*tmx+i], corr[i], ROUNDING);
			mpfr_mul(tmp, tmp, tmp, ROUNDING);
			mpfr_add(cov[i], cov[i], tmp, ROUNDING);
		}
		mpfr_div_d(cov[i], cov[i], (double)(nms-1), ROUNDING);
		mpfr_mul(tmp, corr[0], corr[0], ROUNDING);
		mpfr_div(cov[i], cov[i], tmp, ROUNDING);
	}

	mpfr_clear(tmp);
	free(rd);
}

void get_spectral_density()
{
	double val, ag, mag, bg, mbg;
	double *res, rho, stat, sys;
	mpfr_t estar, ec, tmp0, tmp1, e0;
	FILE *sfp, *rfp;

	mpfr_inits(estar, ec, tmp0, tmp1, e0, NULL);
	mpfr_set_d(e0, emin, ROUNDING);
	res = malloc(nms*sizeof(double));

	if(res == 0)
	{
		error("Failed to allocate auxiliary array");
	}

	rfp = fopen_path("results.txt");
	sfp = fopen_path("systematic.txt");

	for(int n = 0; n < nsteps; n++)
	{
		mpfr_set_d(estar, ei+n*de, ROUNDING);
		rho = 0;
		mag = 0;
		mbg = 0;
		sys = 0;

		printf("\rStarting step %d/%d", n+1, nsteps);
		fflush(stdout);

		for(int k = 0; k < nms; k++)
		{
			bootstrap_sample();

			if(sym)
			{
				rm_method_cosh(e0, estar, cov);
			}
			else
			{
				rm_method(e0, estar, cov);
			}

			transform(e0, estar, corr, cov, &val, &ag, &bg);
			res[k] = val;
			rho += val;
			mag += ag;
			mbg += bg;

			mpfr_set(ec, estar, ROUNDING);
			if(kid == 3)
			{
				mpfr_add_d(ec, ec, sigma, ROUNDING);
			}

			deltabar(tmp0, ec);
			delta(tmp1, estar, ec);
			mpfr_sub(tmp0, tmp0, tmp1, ROUNDING);
			mpfr_div(tmp0, tmp0, tmp1, ROUNDING);

			val = mpfr_get_d(tmp0, ROUNDING);
			val = fabs(val * res[k]);
			sys += val;
		}

		rho /= nms;
		mag /= nms;
		mbg /= nms;
		sys *= 0.6827;
		sys /= nms;
		stat = 0;

		if(nms > 1)
		{
			for(int k = 0; k < nms; k++)
			{
				stat += (res[k]-rho)*(res[k]-rho);
			}
			stat /= (nms-1);
			stat = sqrt(stat);
		}

		fprintf(rfp, "%1.8e %1.8e %1.8e %1.8e %1.8e %1.8e %1.8e\n",
					ei+n*de, mpfr_get_d(expected[n], ROUNDING), rho,
					stat, sys, sqrt(stat*stat+sys*sys), mag+mbg);

		full_sample();

		if(sym)
		{
			rm_method_cosh(e0, estar, cov);
		}
		else
		{
			rm_method(e0, estar, cov);
		}

		for(int k = 0; k < nsteps; k++)
		{
			mpfr_set_d(ec, ei+k*de, ROUNDING);
			deltabar(tmp0, ec);
			delta(tmp1, estar, ec);
			fprintf(sfp, "%1.8e %1.8e %1.8e %1.8e\n",
						ei+n*de, ei+k*de, mpfr_get_d(tmp1, ROUNDING),
						mpfr_get_d(tmp0, ROUNDING));
		}

		fflush(sfp);
		fflush(rfp);
	}

	fclose(rfp);
	fclose(sfp);

	printf("\n");
	mpfr_clears(estar, ec, tmp0, tmp1, e0, NULL);
	free(res);
}

void scan_lambda()
{
	double val, ag, bg, dl;
	mpfr_t estar, e0;
	FILE *fp;

	fp = fopen_path("lambda.txt");

	mpfr_inits(e0, estar, NULL);
	mpfr_set_d(e0, emin, ROUNDING);
	mpfr_set_d(estar, escan, ROUNDING);

	full_sample();
	dl = 0.001;

	for(double l = dl; l < 0.5+dl; l += dl)
	{
		set_params(sigma, l, kid);

		if(sym)
		{
			rm_method_cosh(e0, estar, cov);
		}
		else
		{
			rm_method(e0, estar, cov);
		}

		transform(e0, estar, corr, cov, &val, &ag, &bg);
		fprintf(fp, "%1.8e %1.8e %1.8e %1.8e\n", l, ag, bg, ag+bg);
	}

	mpfr_clears(e0, estar, NULL);
	fclose(fp);
}

void prepare_path()
{
	FILE *fp;

	sprintf(path, "./%ld", time(NULL));
	if(mkdir(path, 0755))
	{
		error("Failed to create directory");
	}

	fp = fopen_path("params.txt");
	fprintf(fp, "T      = %d\n", tmx);
	fprintf(fp, "nms    = %d\n", nms);
	fprintf(fp, "kernel = %d\n", kid);
	fprintf(fp, "prec   = %d\n", prec);
	fprintf(fp, "lambda = %lg\n", lambda);
	fprintf(fp, "sigma  = %lg\n", sigma);
	fprintf(fp, "ei     = %lg\n", ei);
	fprintf(fp, "ef     = %lg\n", ef);
	fprintf(fp, "nsteps = %d\n", nsteps);
	fprintf(fp, "L      = %d\n", L);
	fprintf(fp, "sym    = %d\n", sym);
	fprintf(fp, "emin   = %lg\n", emin);
	fprintf(fp, "scan   = %lg\n", escan);
	fclose(fp);
}

void allocate_data()
{
	int count;

	count = 2*tmx+nsteps+tmx*nms;
	corr = malloc(count*sizeof(mpfr_t));
	cov = corr+tmx;
	expected = cov+tmx;
	mcorr = expected+nsteps;

	if(corr == NULL)
	{
		error("Failed to allocate MPFR variables");
	}

	for(int i = 0; i < count; i++)
	{
		mpfr_init(corr[i]);
	}
}

int main(int argc, char *argv[])
{
	find_int(argc, argv, "-T", &tmx);
	find_int(argc, argv, "-nms", &nms);
	find_int(argc, argv, "-kernel", &kid);

	if(tmx <= 0 || nms <= 0 || kid == -1)
	{
		printf("Usage: %s -T <int> -nms <int> -kernel <int> [ options ]\n", argv[0]);
		printf("Options:\n");
		printf("  -T       <int>     temporal extent of correlator\n");
		printf("  -nms     <int>     number of measurements\n");
		printf("  -lambda  <float>   trade-off parameter (default %1.6f)\n", lambda);
		printf("  -sigma   <float>   width for smearing function (default %1.6f)\n", sigma);
		printf("  -prec    <int>     working precision (default %d decimal places)\n", prec);
		printf("  -ei      <float>   start of energy range for spectral density (default %1.6f)\n", ei);
		printf("  -ef      <float>   end of energy range for spectral density (default %1.6f)\n", ef);
		printf("  -nsteps  <int>     number of steps used in the energy range (default %d)\n", nsteps);
		printf("  -L       <int>     spatial extent used to generate the synthetic correlator (default %d)\n", L);
		printf("  -sym               generatic symmetric (periodic) correlator (not default)\n");
		printf("  -emin    <float>   lower limit for integration over energy (default %1.6f)\n", emin);
		printf("  -scan    <float>   perform a scan in lambda at the specified energy (not default)\n");
		printf("\n");
		printf("Kernels:\n");
		printf("  0        Gaussian\n");
		printf("  1        Sinc\n");
		printf("  2        i-epsilon (real part)\n");
		printf("  3        i-epsilon (imaginary part)\n");
		printf("  4        Bump\n");
		return 0;
	}

	find_dbl(argc, argv, "-sigma", &sigma);
	find_dbl(argc, argv, "-lambda", &lambda);

	find_dbl(argc, argv, "-ei", &ei);
	find_dbl(argc, argv, "-ef", &ef);
	find_int(argc, argv, "-nsteps", &nsteps);
	find_int(argc, argv, "-prec", &prec);
	find_int(argc, argv, "-L", &L);
	find_dbl(argc, argv, "-emin", &emin);
	find_opt(argc, argv, "-sym", &sym);
	find_dbl(argc, argv, "-scan", &escan);

	srand48(1337);
	mpfr_set_default_prec(3.322*prec);

	de = (ef-ei)/nsteps;
	set_params(sigma, lambda, kid);

	if(sym)
	{
		set_time_parms(1, tmx/2, tmx);
	}
	else
	{
		set_time_parms(1, tmx-1, tmx);
	}

	prepare_path();
	allocate_data();

	fakecorr();
	write_fakecorr();
	fakecov();

	if(escan)
	{
		scan_lambda();
	}
	else
	{
		get_spectral_density();
	}

	return 0;
}
