#include "solver.h"
#include <cmath>


namespace solver
{
    BaikalMod::BaikalMod(int TIME, int max_x)
    {
        ozone = new double [max_x];
        q = new double [max_x];
        for (int i = 0; i < max_x; i++)
        {
            ozone[i] = 0.0;
            q[i] = 0.0;
        }
        I = new double [TIME];
        for (int i = 0; i < TIME; i++)
        {
            I[i] = 0.5*exp(-0.0000005*(i-TIME/2)*(i-TIME/2));
        }
    }
    
    BaikalMod::~BaikalMod()
    {
        delete [] q;
        delete [] ozone;
        delete [] I;
    }

    void BaikalMod::iteration(std::ofstream & check, int max_x,
                              int size, double dt, int MOD,
                              double * init_layer, int t)
    {
        //modelling Baikal data
        for (int x = 0; x < max_x; x++)
        {
            double sums = 0.0;
            for (int m = 0; m < size; m++)
            {
                sums += init_layer[m+x*size] * pow(m, 2./3.);
            }
            ozone[x] = ozone[x] + dt * (I[t] - alpha * ozone[x] 
                                        -beta * sums * ozone[x]);
            if (t%MOD==0) check << ozone[x] << " ";
            if (ozone[x] < 0.0)
                ozone[x] = 0.0;
            q[x] = q[x] + dt * (gamma * ozone[x] * (q_max - q[x]) * 
                       (q[x] >= q_max ? 0. : q_max) - delta * q[x]);
            if (q[x] < 0.0)
                q[x] = 0.0;
        }
        if (t%MOD==0) check << "\n";
    }

    AdvDifCoag1d::AdvDifCoag1d()
    {
        ///*
        //res_baikal.txt

        dt = 0.002;
        TIME = 20001;
        MOD = 20;
        size = 64;
        max_x = 500;
        dx = 0.1;
        //*/

        /*
        //res_full.txt

        dt = 0.01;
        TIME = 20001;
        MOD = 2000;
        size = 256;
        max_x = 5000;
        dx = 0.02;
        */

        /*
        //res_ballistic.txt

        dt = 0.001;
        TIME = 50001;
        MOD = 5000;
        size = 256;
        max_x = 5000;
        dx = 0.02;
        */

        /*
        //res_adv.txt
    
        dt = 0.1;
        TIME = 5001; 
        MOD = 500;
        size = 128;
        max_x = 1000;
        dx = 0.2;
        */

        /*
        //res_diff.txt

        dt = 0.2;
        TIME = 2001; 
        MOD = 200;
        size = 256;
        max_x = 500;
        dx = 0.4;
        */

        baikal = new BaikalMod(TIME, max_x);

        n_k = new double [size];
        init_layer = new double [size * max_x];
        smoluch_op = new double [size * max_x];
        vel_coefs = new double [size];
        dif_coefs = new double [size];
        //filling velocity and diffusion coefficients arrays
        for (int m = 0; m < size; m++)
        {
            /*
            //res_diff.txt
            vel_coefs[m] = 0.0;
            dif_coefs[m] = 1.0;
            */

            /*
            //res_adv.txt
            vel_coefs[m] = 1.0;
            dif_coefs[m] = 1.0;
            */

            ///*
            //res_full.txt, res_ballistic.txt or res_baikal.txt
            vel_coefs[m] = 1.0*pow(m+1, 2./3.);
            dif_coefs[m] = 1.0*pow(m+1, -1./3.);
            //*/
        }

        coef_a = new double [max_x * size];
        coef_b = new double [max_x * size];
        coef_c = new double [max_x * size];
        fillMatrix();
    }

    AdvDifCoag1d::~AdvDifCoag1d()
    {
        delete baikal;
        delete [] smoluch_op;
        delete [] init_layer;
        delete [] n_k;
        delete [] vel_coefs;
        delete [] dif_coefs;
        delete [] coef_a;
        delete [] coef_b;
        delete [] coef_c;
    }

    //solve first Smoluchowski integral
    double AdvDifCoag1d::L1(const int &i, double q)
    {
        double l1 = 0;
        for (int i1 = 0; i1 < i; i1++)
        {
            l1 += n_k[i1] * n_k[i - i1 - 1] * K((i - i1 - 1), i1) * 
                        (1.0 + q);
        }
        return l1;
    }

    //solve second Smoluchowski integral
    double AdvDifCoag1d::L2(const int &i, double q)
    {
        double l2 = 0;
        for (int i1 = 0; i1 < size; i1++)
        {
            l2 += n_k[i1] * K(i, i1)*(1.0+q);
        }
        return l2;
    }

    //solve tridiagonal linear system
    void AdvDifCoag1d::solveMatrix(double *a, double *c, double *b,
                                   double *f, double *x)
    {
        double m;
        for (int i = 1; i < max_x; i++)
        {
            m = a[i] / c[i-1];
            c[i] = c[i] - m * b[i-1];
            f[i] = f[i] - m * f[i-1];
        }
        x[max_x-1] = f[max_x-1]/c[max_x-1];

        for (int i = max_x - 2; i >= 0; i--)
        {
            x[i] = ( f[i] - b[i] * x[i+1] ) / c[i];
        }
        return;
    }

    void AdvDifCoag1d::fillMatrix()
    {
        //filling tridiagonal system matrix
        #pragma omp parallel for
        for (int m = 0; m < size; m++)
        {
            for (int x = 0; x < max_x; x++)
            {
                int ind = x + m * max_x;
                coef_b[ind] = - dt * dif_coefs[m]
                            * exp(-vel_coefs[m]/dif_coefs[m]*dx/2.0)
                            / (dx*dx);
                coef_c[ind] = 1.0 + dt * 
                              (dif_coefs[m] * 
                               exp(vel_coefs[m]/dif_coefs[m]*dx/2.0) +
                               dif_coefs[m] * 
                               exp(-vel_coefs[m]/dif_coefs[m]*dx/2.0)) 
                              / (dx*dx);
                coef_a[ind] = - dt * dif_coefs[m] * 
                            exp(vel_coefs[m]/dif_coefs[m]*dx/2.0)
                            / (dx*dx);
            }

            //boundary conditions can be specified here
            coef_b[0 + m*max_x] = 1.0;

            coef_b[(max_x-1) + m*max_x] = 0.0;
            coef_a[0 + m*max_x] = 0.0;
            coef_a[(max_x-1) + m*max_x] = 0.0;

            coef_c[0 + m*max_x] = -1.0;

            coef_c[(max_x-1) + m*max_x] = 1.0;
        }
    }

    void AdvDifCoag1d::iteration(std::ofstream & output,
                                 std::ofstream & check, int t)
    {
        //this is calculating coagulation
        for (int x = 0; x < max_x; x++)
        {
            for (int m = 0; m < size; m++)
            {
                int ind = m+x*size;
                n_k[m] = init_layer[ind];
                if (t%MOD==0 && x%XMOD==0) output 
                                           << init_layer[ind] << " ";
            }
            if (t%MOD==0 && x%XMOD==0) output << std::endl;

            #pragma omp parallel for
            for (int m = 0; m < size; m++)
            {
            	int ind = m+x*size;
                smoluch_op[ind] = (L1(m, baikal->q[x]) * 0.5
                     -n_k[m] * L2(m, baikal->q[x])) * dt + n_k[m];
                if (smoluch_op[ind] < 0.0) smoluch_op[ind] = 0.0;
            }
        }
        //this is calculating advection
        #pragma omp parallel for
        for (int m = 0; m < size; m++)
        {
            double * rhs = new double [max_x];
            for (int x = 0; x < max_x; x++)
            {
                //boundary conditions can be specified here
                rhs[x] = smoluch_op[m+x*size];
                if (m==0 && x==0) rhs[x] = -J*dx*0.5;
                
                //uncomment this line if res_diff.txt else comment it
                //else if (x==0) rhs[x] = 0.0;
            }
            double * output = new double [max_x];
            solveMatrix(&coef_a[m*max_x], &coef_c[m*max_x], 
                        &coef_b[m*max_x], rhs, output);
            
            for (int x = 0; x < max_x; x++)
            {
                init_layer[m+x*size] = output[x];
            }
            delete [] output;
            delete [] rhs;
        }
        baikal->iteration(check, max_x, size, dt, MOD, init_layer, t);
        fillMatrix();
    }

    double AdvDifCoag1d::K(const int & u, const int &v)
    {
        //simplified kernel for res_full.txt
        //double u1=pow( (u + 1.0) , 2.0/3.0);
        //double v1=pow( (v + 1.0) , 2.0/3.0);
        //double result = u1*v1;

        //ballistic kernel for res_ballistic.txt and res_baikal.txt
        double u1=pow( (u + 1.0) , 1.0/3.0);
        double v1=pow( (v + 1.0) , 1.0/3.0);
        double result = 0.0;
        if (u==v)
            result = (u1+v1)*(1./u1+1./v1);
        else
            result = pow((u1+v1), 2)*fabs(u1*u1-v1*v1);
        
        //kernel for res_diff.txt and res_adv.txt
        //double result = 2.0;
        return result;
    }
}
