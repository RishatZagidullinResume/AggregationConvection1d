#pragma once
#include <stdio.h>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>

namespace solver
{
    class BaikalMod
    {
    protected:
        double alpha = 1.0;
        double beta = 3.0;
        double gamma = 5.0;
        double q_max = 5.0;
        double delta = 3.0;
        double * ozone, *I;
    public:
        double * q;
        BaikalMod(int TIME, int max_x);
        ~BaikalMod();
        void iteration(std::ofstream & check, int max_x, int size,
                       double dt, int MOD,
                       double * init_layer, int t);
    };

    class AdvDifCoag1d
    {
    protected:
        double dt, dx;
        int MOD;

        double J = 1.0;
        int XMOD = 25;

        double *n_k, *smoluch_op;
        double *vel_coefs, *dif_coefs;
        double *coef_a, *coef_b, *coef_c;

        BaikalMod * baikal;
    public:
        int TIME, size, max_x;
        double *init_layer;
        AdvDifCoag1d();
        ~AdvDifCoag1d();
        double L1(const int &i, double q);
        double L2(const int &i, double q);
        void solveMatrix(double *a, double *c, double *b, 
                         double *f, double *x);
        void iteration(std::ofstream & output, std::ofstream & check,
                       int t);
        void fillMatrix();
        double K(const int & u, const int &v);
    };


}

