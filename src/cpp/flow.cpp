//
// Created by mafu on 24-3-26.
//
#include "../headers/flow.h"
#include "gsl/gsl_sf_bessel.h"
#include <iostream>

int Coal::equation_multi(const gsl_vector *x, void *p, gsl_vector *f) {
    const params *params = static_cast<struct params *>(p);
    const double a     = params->a;
    const int k        = params->k;
    const double x_val = gsl_vector_get(x, 0);
    const double arg   = x_val * x_val / 4;

    const double I_k_minus_1_over_2 = gsl_sf_bessel_Inu((k - 1) / 2.0, arg);
    const double I_k_plus_1_over_2  = gsl_sf_bessel_Inu((k + 1) / 2.0, arg);

    const double y = M_PI / 2 / sqrt(2) * x_val * exp(-arg) * (I_k_minus_1_over_2 + I_k_plus_1_over_2) - a;

    gsl_vector_set(f, 0, y);
    return GSL_SUCCESS;
}
double Coal::solveEquation(const double a,const int k) {
    const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver *s            = gsl_multiroot_fsolver_alloc(T, 1);

    params par                  = {a, k};
    gsl_multiroot_function func = {&equation_multi, 1, &par};

     constexpr double x_init[1] = {0.5};
    gsl_vector *x    = gsl_vector_alloc(1);
    gsl_vector_set(x, 0, x_init[0]);

    gsl_multiroot_fsolver_set(s, &func, x);

    int status, iter = 0;
    do {
        iter++;
        status = gsl_multiroot_fsolver_iterate(s);
        if (status)
            break;

        status = gsl_multiroot_test_residual(s->f, 1e-7);
    } while (status == GSL_CONTINUE && iter < 1000);

    const double root = gsl_vector_get(s->x, 0);
    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);

    return root;
}
double Coal::getRef(const double x, const int k) {
    const double arg = x * x/ 4;
    const double I_k_minus_1_over_2 = gsl_sf_bessel_Inu((k - 1) / 2.0, arg);
    const double I_k_plus_1_over_2  = gsl_sf_bessel_Inu((k + 1) / 2.0, arg);
    return M_PI / 2 / sqrt(2) * x * std::exp(-arg) * (I_k_minus_1_over_2 + I_k_plus_1_over_2);
}
