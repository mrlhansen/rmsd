#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"

void find_dbl(int argc, char **argv, char *name, double *value)
{
	for(int k = 1; k < argc; k++)
	{
		if(strcmp(argv[k], name) == 0)
		{
			if(k+1 < argc)
			{
				sscanf(argv[k+1], "%lg", value);
				return;
			}
			else
			{
				error("Missing value for command line option");
			}
		}
	}
}

void find_int(int argc, char **argv, char *name, int *value)
{
	for(int k = 1; k < argc; k++)
	{
		if(strcmp(argv[k], name) == 0)
		{
			if(k+1 < argc)
			{
				sscanf(argv[k+1], "%d", value);
				return;
			}
			else
			{
				error("Missing value for command line option");
			}
		}
	}
}

void find_str(int argc, char **argv, char *name, char *value)
{
	for(int k = 1; k < argc; k++)
	{
		if(strcmp(argv[k], name) == 0)
		{
			if(k+1 < argc)
			{
				sscanf(argv[k+1], "%s", value);
				return;
			}
			else
			{
				error("Missing value for command line option");
			}
		}
	}
}

void find_opt(int argc, char **argv, char *name, int *value)
{
	for(int k = 1; k < argc; k++)
	{
		if(strcmp(argv[k], name) == 0)
		{
			*value = 1;
			return;
		}
	}
	*value = 0;
}

void randu(double *rd, int n)
{
	for(int i = 0; i < n; i++)
	{
		rd[i] = drand48();
	}
}

void randi(int *rd, int n, int a, int b)
{
	for(int i = 0; i < n; i++)
	{
		rd[i] = a + drand48()*(b-a);
	}
}

void randg(double *rd, int n)
{
	double u1, u2;
	double r1, r2;

	for(int i = 0; i < n;)
	{
		u1 = drand48();
		u1 = -2.0*log(u1);
		u1 = sqrt(u1);
		u2 = 2.0*M_PI*drand48();
		r1 = u1*cos(u2);
		r2 = u1*sin(u2);
		rd[i++] = r1;
		if(i < n)
		{
			rd[i++] = r2;
		}
	}
}
