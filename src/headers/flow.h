//
// Created by mafu on 24-3-26.
//
#pragma once
#include "gsl/gsl_multiroots.h"
namespace Coal {

    struct params {
        double a;
        int k;
    };

    double equation(double x, void *p);
    int equation_multi(const gsl_vector *x, void *p, gsl_vector *f);

    double solveEquation(double a, int k);

    double getRef(double x, int k);

}// namespace Coal
