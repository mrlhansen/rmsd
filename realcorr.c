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
int kid = -1;

mpfr_t *corr, *cov;
mpfr_t *expected;
mpfr_t *mcorr;

char path[128];
char file[128];
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

void read_corr()
{
	double res, dres, norm;
	char *token;
	char *buf;
	int offset;
	FILE *fp;

	fp = fopen(file, "r");
	if(fp == NULL)
	{
		error("Unable to open file for reading");
	}

	buf = malloc(2048*sizeof(char));
	offset = 0;

	while(fgets(buf, 2048, fp))
	{
		token = strtok(buf, " ");

		if(strlen(buf) == 0)
		{
			break;
		}

		while(token)
		{
			mpfr_set_d(mcorr[offset], atof(token), ROUNDING);
			token = strtok(NULL, " ");
			offset++;
		}
	}

	if(offset != tmx*nms)
	{
		error("Error reading correlator data, number of elements is incorrect");
	}

	fclose(fp);
	fp = fopen_path("correlator.txt");
	full_sample();
	norm = mpfr_get_d(corr[0], ROUNDING);

	for(int i = 0; i < tmx; i++)
	{
		res = mpfr_get_d(corr[i], ROUNDING);
		dres = mpfr_get_d(cov[i], ROUNDING);
		dres = norm*sqrt(dres);
		fprintf(fp, "%3d %1.8e %1.8e\n", i, res, dres);
	}

	fclose(fp);
	free(buf);
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
			rm_method_cosh(e0, estar, cov);

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

		fprintf(rfp, "%1.8e %1.8e %1.8e %1.8e %1.8e %1.8e\n",
					ei+n*de, rho, stat, sys, sqrt(stat*stat+sys*sys), mag+mbg);

		full_sample();
		rm_method_cosh(e0, estar, cov);

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
		rm_method_cosh(e0, estar, cov);
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
	fprintf(fp, "file   = %s\n", file);
	fprintf(fp, "prec   = %d\n", prec);
	fprintf(fp, "lambda = %lg\n", lambda);
	fprintf(fp, "sigma  = %lg\n", sigma);
	fprintf(fp, "ei     = %lg\n", ei);
	fprintf(fp, "ef     = %lg\n", ef);
	fprintf(fp, "nsteps = %d\n", nsteps);
	fprintf(fp, "emin   = %lg\n", emin);
	fprintf(fp, "scan   = %lg\n", escan);
	fclose(fp);
}

void allocate_data()
{
	int count;

	count = 2*tmx+tmx*nms;
	corr = malloc(count*sizeof(mpfr_t));
	cov = corr+tmx;
	mcorr = cov+tmx;

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
	memset(file, 0, sizeof(file));
	find_int(argc, argv, "-T", &tmx);
	find_int(argc, argv, "-nms", &nms);
	find_int(argc, argv, "-kernel", &kid);
	find_str(argc, argv, "-file", file);

	if(tmx <= 0 || nms <= 0 || kid == -1 || strlen(file) == 0)
	{
		printf("Usage: %s -T <int> -nms <int> -kernel <int> -file <string> [ options ]\n", argv[0]);
		printf("Options:\n");
		printf("  -T       <int>     temporal extent of correlator\n");
		printf("  -nms     <int>     number of measurements\n");
		printf("  -file    <string>  path for reading the correlator data\n");
		printf("  -lambda  <float>   trade-off parameter (default %1.6f)\n", lambda);
		printf("  -sigma   <float>   width for smearing function (default %1.6f)\n", sigma);
		printf("  -prec    <int>     working precision (default %d decimal places)\n", prec);
		printf("  -ei      <float>   start of energy range for spectral density (default %1.6f)\n", ei);
		printf("  -ef      <float>   end of energy range for spectral density (default %1.6f)\n", ef);
		printf("  -nsteps  <int>     number of steps used in the energy range (default %d)\n", nsteps);
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
	find_dbl(argc, argv, "-emin", &emin);
	find_dbl(argc, argv, "-scan", &escan);

	srand48(1337);
	mpfr_set_default_prec(3.322*prec);

	de = (ef-ei)/nsteps;
	set_params(sigma, lambda, kid);
	set_time_parms(1, tmx/2, tmx);

	prepare_path();
	allocate_data();
	read_corr();

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
