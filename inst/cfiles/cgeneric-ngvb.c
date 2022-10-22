
#include <assert.h>
#if !defined(__FreeBSD__)
#include <malloc.h>
#endif
#include <math.h>
#include <stdlib.h>
#include <strings.h>

#include "cgeneric.h"

#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#define SQR(x) ((x)*(x))

double logprior_ar1(double tau, double phi, double U1, double alpha1, double U2, double alpha2){
	double lambda = -log(alpha1)/U1;
	double log_tau_prior = log(lambda/2) -1.5*log(tau) - lambda/sqrt(tau);

	double theta = -log(alpha2)/sqrt(-log(1-SQR(U2)));
	double log_phi_prior = -log(2.0) + log(theta) - theta*sqrt(-log(1-SQR(phi))) + log(fabs(phi)) - log(1-SQR(phi)) -0.5*log(-log(1-SQR(phi)));
	return (log_tau_prior + log_phi_prior);
}
//prior for tau (driving noise precision) is Gamma(alpha,beta)
//Pior for autocorrelation parameter is given by 3.1 inhttps://repository.kaust.edu.sa/bitstream/handle/10754/625107/revision-sorbye-rue-jtsa.pdf?sequence=1&isAllowed=y
//P(|phi|>U)=a
//cmodel <- inla.cgeneric.define(model = "inla_cgeneric_ar1_model",
//                               shlib = "cgeneric-VINIG.so", n = n, 
//                               V = rep(1,n),
//                               alpha = alpha,
//                               beta  = beta,
//                               gamma = -log(a)/sqrt(-log(1-U^2)))

double *inla_cgeneric_ar1_model(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data)
{
	// this reimplement `inla.rgeneric.ar1.model` using cgeneric

	double *ret = NULL, prec, lprec, rho, rho_intern;

	if (theta) {
		lprec = theta[0];
		prec = exp(lprec);
		rho_intern = theta[1];
		rho = 2.0 * exp(rho_intern) / (1.0 + exp(rho_intern)) - 1.0;
	} else {
		prec = lprec = rho = rho_intern = NAN;
	}

	assert(!strcasecmp(data->ints[0]->name, "n"));	       // this will always be the case
	int N = data->ints[0]->ints[0];			       // this will always be the case
	assert(N > 0);

	assert(data->n_doubles > 0);
	assert(!strcasecmp(data->doubles[0]->name, "V"));
	inla_cgeneric_vec_tp *V = data->doubles[0];

	assert(!strcasecmp(data->doubles[1]->name, "U1"));
	inla_cgeneric_vec_tp *U1 = data->doubles[1];
	assert(!strcasecmp(data->doubles[2]->name, "alpha1"));
	inla_cgeneric_vec_tp *alpha1  = data->doubles[2];
	assert(!strcasecmp(data->doubles[3]->name, "U2"));
    inla_cgeneric_vec_tp *U2 = data->doubles[3];
	assert(!strcasecmp(data->doubles[4]->name, "alpha2"));
	inla_cgeneric_vec_tp *alpha2 = data->doubles[4];

	switch (cmd) {
	case INLA_CGENERIC_VOID:
	{
		assert(!(cmd == INLA_CGENERIC_VOID));
		break;
	}

	case INLA_CGENERIC_GRAPH:
	{
		// return a vector of indices with format
		// c(N, M, ii, jj)
		// where ii<=jj, ii is non-decreasing and jj is non-decreasing for the same ii
		// so like the loop
		// for i=0, ...
		// for j=i, ...
		// G_ij = 
		// and M is the length of ii

		int M = N + N - 1, offset, i, k;
		ret = Calloc(2 + 2 * M, double);
		assert(ret);
		offset = 2;
		ret[0] = N;				       /* dimension */
		ret[1] = M;				       /* number of (i <= j) */
		for (k = i = 0; i < N; i++) {
			ret[offset + k] = i;		       /* i */
			ret[offset + M + k++] = i;	       /* j */
			if (i < N - 1) {
				ret[offset + k] = i;	       /* i */
				ret[offset + M + k++] = i + 1; /* j */
			}
		}
		break;
	}

	case INLA_CGENERIC_Q:
	{
		// optimized format
		// return c(-1, M, Qij) in the same order as defined in INLA_CGENERIC_GRAPH
		// where M is the length of Qij

		double param = prec / (1.0 - SQR(rho));
		int M = N + N - 1;
		int offset, i, k;
		ret = Calloc(2 + M, double);
		assert(ret);
		offset = 2;
		ret[0] = -1;				       /* REQUIRED */
		ret[1] = M;
		for (i = k = 0; i < N; i++) {
			ret[offset + k++] = param * (i == 0 || i == N - 1 ? (V->doubles[i]) : ((V->doubles[i-1]) + (V->doubles[i])*SQR(rho)));
			if (i < N - 1) {
				ret[offset + k++] = -param * rho * (V->doubles[i]);
			}
		}
		break;
	}

	case INLA_CGENERIC_MU:
	{
		// return (N, mu)
		// if N==0 then mu is not needed as its taken to be mu[]==0

		ret = Calloc(1, double);
		assert(ret);
		ret[0] = 0;
		break;
	}

	case INLA_CGENERIC_INITIAL:
	{
		// return c(M, initials)
		// where M is the number of hyperparameters

		ret = Calloc(3, double);
		assert(ret);
		ret[0] = 2.0;
		ret[1] = 1.0;
		ret[2] = 1.0;
		break;
	}

	case INLA_CGENERIC_LOG_NORM_CONST:
	{
		// return c(NORM_CONST) or a NULL-pointer if INLA should compute it by itself

		double prec_innovation = prec / (1.0 - SQR(rho));
		ret = Calloc(1, double);
		assert(ret);
		ret[0] = N * (-0.5 * log(2.0 * M_PI) + 0.5 * log(prec_innovation)) + 0.5 * log(1.0 - SQR(rho));
		break;
	}

	case INLA_CGENERIC_LOG_PRIOR:
	{
		// return c(LOG_PRIOR)


		ret = Calloc(1, double);
		assert(ret);
		ret[0] = logprior_ar1(prec, rho, U1->doubles[0], alpha1->doubles[0], U2->doubles[0], alpha2->doubles[0]) + lprec +  log(2.0) + rho_intern - 2.0*log(1.0 + exp(rho_intern));
		//ret[0] = -prec + lprec - 0.5 * log(2.0 * M_PI) - 0.5  * SQR(rho_intern);
		break;
	}

	case INLA_CGENERIC_QUIT:
	default:
		break;
	}

	return (ret);
}