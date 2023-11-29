
#include<iostream>
#include<cmath>
#include <stdlib.h>
#include <stdio.h>

double V[1000000];
double Vu[1000000];
double Vup[1000000];
double V_right[1000000];
double V_right_tmp[1000000];
double V_tmp[1000000];
double Vup_tmp[1000000];

void write_2vectors_file(FILE* file, double* vector1, double* vector2, size_t size_vector) {

    for (size_t i = 0; i < size_vector; i++)
    {
        fprintf(file, "%f  %f  \n", vector1[i], vector2[i]);
    }fprintf(file, "\n");
    return;
}

void equalize_vector(double* vector1, const double* vector2, size_t size_vector)
{
    for (size_t i = 0; i < size_vector; i++)
    {
        vector1[i] = vector2[i];
    }
    return;
}


double linear_solution(double t, double x)
{
    if (t > 0. && x >= -1. && x <= -t / 2.) { return 0; }
    else if (x >= -t / 2. && x <= 1. / 4 - t / 2.) { return 2 + 4 * x; }
    else if (x >= 1. / 4 - t / 2. && x <= 1.) { return 1.; }
    return 0.;
}



int calculate_linear_explicit_write_file(FILE* file, FILE* file_table, int Ntau, int Mh)
{
    double h = 2. / Mh;
    double tau = 1. / Ntau;
    double nu = tau / h;
    for (int i = 0; i < Mh + 1; i++)
    {
        if (-1. + i * h <= 0) { V[i] = 0.; }
        else if (-1. + i * h <= 1. / 4) { V[i] = 4 * (-1. + i * h); }
        else V[i] = 1.;
    }
    for (int j = 1; j < Ntau + 1; j++)
    {

        Vup[0] = 0.;
        Vup[Mh] = 1.;
        for (int i = 1; i < Mh; i++)
        {
            Vup[i] = V[i] + (nu / 2) * (V[i + 1] - V[i]);

        }
        equalize_vector(V, Vup, Mh + 1);
    }

    double Delta_V_Ch = 0.;
    double V_Ch = 0.;

    double Delta_V_L1h = 0.;
    double V_L1h = 0.;

    for (int i = 0; i < Mh + 1; i++)
    {
        if (fabs(Vup[i] - linear_solution(1., -1. + i * h)) > Delta_V_Ch)
        {
            Delta_V_Ch = fabs(Vup[i] - linear_solution(1, -1. + i * h));
        }
        if (fabs(Vup[i]) > V_Ch)
        {
            V_Ch = fabs(Vup[i]);
        }
        Delta_V_L1h = Delta_V_L1h + fabs(Vup[i] - linear_solution(1., -1. + i * h));

        V_L1h = V_L1h + fabs(Vup[i]);

    }
    Delta_V_L1h = Delta_V_L1h * h;
    V_L1h = V_L1h * h;

    for (int i = 0; i < Mh + 1; i++) { Vup[i] = -1 + i * h; }
    write_2vectors_file(file, Vup, V, Mh + 1);
    fprintf(file_table, "%.3f %.3f %e %e %e %e \n", tau, h, Delta_V_Ch, Delta_V_L1h, Delta_V_Ch / V_Ch, Delta_V_L1h / V_L1h);
    return 0;
}

int calculate_linear_explicit_write_file_local(FILE* file_table, int Ntau, int Mh, int K)
{
    double h = 2. / Mh;
    double tau = 1. / Ntau;
    double nu = tau / h;

    for (int i = 0; i < Mh + 1; i++)
    {
        if (-1. + i * h <= 0) { V[i] = 0.; }
        else if (-1. + i * h <= 1. / 4) { V[i] = 4 * (-1. + i * h); }
        else V[i] = 1.;
    }

    for (int j = 1; j < Ntau + 1; j++)
    {
        Vup[0] = 0.;
        Vup[Mh] = 1.;

        for (int i = 1; i < Mh; i++)
        {
            Vup[i] = V[i] + nu / 2. * (V[i + 1] - V[i]);
        }
        equalize_vector(V, Vup, Mh + 1);
    }

    int Mh_tmp = Mh;
    int Ntau_tmp = Ntau;
    int n = 1;

    for (int k = 0; k < K; k++)
    {

        Mh_tmp = Mh_tmp * 2;
        Ntau_tmp = Ntau_tmp * 2;
        n = n * 2;

        double h_tmp = 2. / Mh_tmp;
        double tau_tmp = 1. / Ntau_tmp;
        double nu_tmp = tau_tmp / h_tmp;

        for (int i = 0; i < Mh_tmp + 1; i++)
        {
            if (-1. + i * h_tmp <= 0) { V_tmp[i] = 0.; }
            else if (-1. + i * h_tmp <= 1. / 4) { V_tmp[i] = 4 * (-1. + i * h_tmp); }
            else V_tmp[i] = 1.;
        }

        for (int j = 1; j < Ntau_tmp + 1; j++)
        {
            Vup_tmp[0] = 0.;
            Vup_tmp[Mh_tmp] = 1.;

            for (int i = 1; i < Mh_tmp; i++)
            {
                Vup_tmp[i] = V_tmp[i] + nu_tmp / 2. * (V_tmp[i + 1] - V_tmp[i]);
            }
            equalize_vector(V_tmp, Vup_tmp, Mh_tmp + 1);
        }

        

        double Delta_V_Ch = 0.;
        double V_Ch = 0.;

        double Delta_V_L1h = 0.;
        double V_L1h = 0.;

        for (int i = 0; i < Mh + 1; i++)
        {
            if (fabs(V[i] - V_tmp[i * n]) > Delta_V_Ch)
            {
                Delta_V_Ch = fabs(V[i] - V_tmp[i * n]);
            }
            if (fabs(V[i]) > V_Ch)
            {
                V_Ch = fabs(V[i]);
            }
            Delta_V_L1h = Delta_V_L1h + fabs(V[i] - V_tmp[i * n]);
            V_L1h = V_L1h + fabs(V[i]);
        }
        Delta_V_L1h = Delta_V_L1h * h;
        V_L1h = V_L1h * h;
       // for (int i = 0; i < Mh_tmp + 1; i++) { Vu[i] = -1 + i * h_tmp; }
       // write_2vectors_file(file_table, Vu,V_tmp, Mh_tmp + 1);
        fprintf(file_table, "%f %f %e %e %e %e \n", tau_tmp, h_tmp, Delta_V_Ch, Delta_V_L1h, Delta_V_Ch / V_Ch, Delta_V_L1h / V_L1h);
        
        
    }
    return 0;
}

int sweep(double tau, double h, const int N, double* V_right, double* Vup)
{
    double* A = new double[N];
    double* B = new double[N];


    double mu = tau / (4. * h);

    A[0] = 0.;
    A[1] = mu;

    B[0] = 0;
    B[1] = V_right[0];

    for (int i = 2; i < N; i++)
    {
        A[i] = mu / (mu * A[i - 1] + 1);
        B[i] = (V_right[i - 1] - mu * B[i - 1]) / (mu * A[i - 1] + 1);
    }

    Vup[N - 1] = -(B[N - 1] - V_right[N - 1] / mu) / (A[N - 1] + 1 / mu);

    for (int i = N - 2; i >= 0; i--)
    {
        Vup[i] = A[i + 1] * Vup[i + 1] + B[i + 1];
    }
    delete[] A;
    delete[] B;
    return 0;
}


int calculate_linear_implicit_write_file(FILE* file, FILE* file_table, int Ntau, int Mh, double w)
{
    double h = 2. / Mh;
    double tau = 1. / Ntau;
    double nu = tau / h;
    double mu = tau / (4 * h);

    for (int i = 0; i < Mh + 1; i++)
    {
        if (-1. + i * h <= 0) { V[i] = 0.; }
        else if (-1. + i * h <= 1. / 4) { V[i] = 4 * (-1. + i * h); }
        else V[i] = 1.;
    }
    for (int j = 1; j < Ntau + 1; j++)
    {
        Vup[0] = 0.;
        Vup[Mh] = 1.;
        {
            V_right[1] = w * V[0] + (1 - 2 * w) * V[1] + w * V[2] - mu * Vup[0];
            for (int i = 2; i < Mh - 1; i++)
            {
                V_right[i] = w * V[i - 1] + (1 - 2 * w) * V[i] + w * V[i + 1];
            }
            V_right[Mh - 1] = w * V[Mh - 2] + (1 - 2 * w) * V[Mh - 1] + w * V[Mh] + mu * Vup[Mh];
        }
        sweep(tau, h, Mh - 1, V_right + 1, Vup + 1);
        equalize_vector(V, Vup, Mh + 1);
    }

    double Delta_V_Ch = 0.;
    double V_Ch = 0.;

    double Delta_V_L1h = 0.;
    double V_L1h = 0.;

    for (int i = 0; i < Mh + 1; i++)
    {
        if (fabs(Vup[i] - linear_solution(1., -1. + i * h)) > Delta_V_Ch)
        {
            Delta_V_Ch = fabs(Vup[i] - linear_solution(1., -1. + i * h));
        }
        if (fabs(Vup[i]) > V_Ch)
        {
            V_Ch = fabs(Vup[i]);
        }
        Delta_V_L1h = Delta_V_L1h + fabs(Vup[i] - linear_solution(1., -1. + i * h));
        V_L1h = V_L1h + fabs(Vup[i]);
    }
    Delta_V_L1h = Delta_V_L1h * h;
    V_L1h = V_L1h * h;
    for (int i = 0; i < Mh + 1; i++) { Vup[i] = -1 + i * h; }
    write_2vectors_file(file, Vup, V, Mh + 1);

    fprintf(file_table, "%.3f %.3f %e %e %e %e \n", tau, h, Delta_V_Ch, Delta_V_L1h, Delta_V_Ch / V_Ch, Delta_V_L1h / V_L1h);
    return 0;
}



int calculate_linear_implicit_write_file_local(FILE* file_table, int Ntau, int Mh, int K, double w)
{
    double h = 2. / Mh;
    double tau = 1. / Ntau;
    double nu = tau / h;
    double mu = tau / (4 * h);

    for (int i = 0; i < Mh + 1; i++)
    {
        if (-1. + i * h <= 0) { V[i] = 0.; }
        else if (-1. + i * h <= 1. / 4) { V[i] = 4 * (-1. + i * h); }
        else V[i] = 1.;
    }
    for (int j = 1; j < Ntau + 1; j++)
    {
        Vup[0] = 0.;
        Vup[Mh] = 1.;
        {
            V_right[1] = w * V[0] + (1 - 2 * w) * V[1] + w * V[2] - mu * Vup[0];
            for (int i = 2; i < Mh - 1; i++)
            {
                V_right[i] = w * V[i - 1] + (1 - 2 * w) * V[i] + w * V[i + 1];
            }
            V_right[Mh - 1] = w * V[Mh - 2] + (1 - 2 * w) * V[Mh - 1] + w * V[Mh] + mu * Vup[Mh];
        }
        sweep(tau, h, Mh - 1, V_right + 1, Vup + 1);
        equalize_vector(V, Vup, Mh + 1);
    }
    int Mh_tmp = Mh;
    int Ntau_tmp = Ntau;
    int s = 1.;

    for (int k = 0; k < K; k++)
    {

        Mh_tmp = Mh_tmp * 2;
        Ntau_tmp = Ntau_tmp * 2;
        s = s * 2;
        double h_tmp = 2. / Mh_tmp;
        double tau_tmp = 1. / Ntau_tmp;
        double nu_tmp = tau_tmp / h_tmp;
        double mu_tmp = tau_tmp / (4 * h_tmp);

        for (int i = 0; i < Mh_tmp + 1; i++)
        {
            if (-1. + i * h_tmp <= 0) { V_tmp[i] = 0.; }
            else if (-1. + i * h_tmp <= 1. / 4) { V_tmp[i] = 4 * (-1. + i * h_tmp); }
            else V_tmp[i] = 1.;
        }

        for (int j = 1; j < Ntau_tmp + 1; j++)
        {
            Vup_tmp[0] = 0.;
            Vup_tmp[Mh_tmp] = 1.;
            {
                V_right_tmp[1] = w * V_tmp[0] + (1 - 2 * w) * V_tmp[1] + w * V_tmp[2] - mu_tmp * Vup_tmp[0];
                for (int i = 2; i < Mh_tmp - 1; i++)
                {
                    V_right_tmp[i] = w * V_tmp[i - 1] + (1 - 2 * w) * V_tmp[i] + w * V_tmp[i + 1];
                }
                V_right_tmp[Mh_tmp - 1] = w * V_tmp[Mh_tmp - 2] + (1 - 2 * w) * V_tmp[Mh_tmp - 1] + w * V_tmp[Mh_tmp] + mu_tmp * Vup_tmp[Mh_tmp];
            }
            sweep(tau_tmp, h_tmp, Mh_tmp - 1, V_right_tmp + 1, Vup_tmp + 1);
            equalize_vector(V_tmp, Vup_tmp, Mh_tmp + 1);
        }


        double Delta_V_Ch = 0.;
        double V_Ch = 0.;

        double Delta_V_L1h = 0.;
        double V_L1h = 0.;

        for (int i = 0; i < Mh + 1; i++)
        {
            if (fabs(V[i] - V_tmp[i * s]) > Delta_V_Ch)
            {
                Delta_V_Ch = fabs(V[i] - V_tmp[i * s]);
            }
            if (fabs(V[i]) > V_Ch)
            {
                V_Ch = fabs(V[i]);
            }
            Delta_V_L1h = Delta_V_L1h + fabs(V[i] - V_tmp[i * s]);
            V_L1h = V_L1h + fabs(V[i]);

        }
        Delta_V_L1h = Delta_V_L1h * h;
        V_L1h = V_L1h * h;
       // for (int i = 0; i < Mh_tmp + 1; i++) { Vu[i] = -1 + i * h_tmp; }
       // write_2vectors_file(file_table, Vu, V_tmp, Mh_tmp + 1);
        fprintf(file_table, "%f & %f & %e & %e & %e & %e \n", tau_tmp, h_tmp, Delta_V_Ch, Delta_V_L1h, Delta_V_Ch / V_Ch, Delta_V_L1h / V_L1h);

    }

    return 0;
}

double nonlinear_solution(double t, double x)
{
    if (t > 0. && t < 1. / 4 && x <= 0. && x >= -1) return 0.;
    else if (t > 1. / 4 && x <= (-1. / 2) * (t - 1. / 4)) return 0.;
    else if (t > 0. && t < 1.4 && x <= 1. / 4 - t)return 4 * x / (1 - 4 * t);
    else if (t > 1. / 4 && x > (- 1. / 2) * (t - 1. / 4)) return 1.;
    else if (t > 0 && t < 1. / 4 && x>1. / 4 - t) return 1.;
    return 0.;
}


int calculate_nonlinear_explicit_write_file(FILE* file, FILE* file_table, int Ntau, int Mh)
{
    double h = 2. / Mh;
    double tau = 1. / Ntau;
    double nu = tau / h;
    for (int i = 0; i < Mh + 1; i++)

    {
        if (-1. + i * h <= 0) { V[i] = 0.; }
        else if (-1. + i * h <= 1. / 4) { V[i] = 4 * (-1. + i * h); }
        else V[i] = 1.;
    }

    for (int j = 1; j < Ntau + 1; j++)
    {
        Vup[0] = 0.;
        Vup[Mh] = 1.;

        for (int i = 1; i < Mh; i++)
        {
            Vup[i] = V[i] + (nu / 5) * (pow(V[i + 1], 2) - pow(V[i], 2));
          //  Vup[i] = V[i] + (nu / 5) * (pow(V[i + 1], 2) - pow(V[i], 2));
        }

        equalize_vector(V, Vup, Mh + 1);
    }

    double Delta_V_Ch = 0.;
    double V_Ch = 0.;

    double Delta_V_L1h = 0.;
    double V_L1h = 0.;

    for (int i = 0; i < Mh + 1; i++)
    {

        if (fabs(Vup[i] - nonlinear_solution(1., -1. + i * h)) > Delta_V_Ch) 
        {
            Delta_V_Ch = fabs(Vup[i] - nonlinear_solution(1., -1. + i * h));
        }
        if (fabs(Vup[i]) > V_Ch)
        {
            V_Ch = fabs(Vup[i]);
        }
        Delta_V_L1h = Delta_V_L1h + fabs(Vup[i] - nonlinear_solution(1., -1. + i * h));
        V_L1h = V_L1h + fabs(Vup[i]);
    }
    Delta_V_Ch = Delta_V_Ch * 0.4;
    Delta_V_L1h = Delta_V_L1h * h;
    V_L1h = V_L1h * h;

    for (int i = 0; i < Mh + 1; i++) { V[i] = -1 + i * h; }
    write_2vectors_file(file, Vup, V, Mh + 1);

    fprintf(file_table, " %.3f & %.3f & %e & %e & %e & %e  \n", tau, h, Delta_V_Ch, Delta_V_L1h, Delta_V_Ch / V_Ch, Delta_V_L1h / V_L1h);
    return 0;
}

/*int calculate_nonlinear_exlicit_write_file(FILE* file, FILE* file_table, int Ntau, int Mh)
{
    double h = 2. / Mh;
    double tau = 1. / Ntau;
    double nu = tau / h;
    for (int i = 0; i < Mh + 1; i++)

    {
        if (-1. + i * h <= 0) { V[i] = 0.; }
        else if (-1. + i * h <= 1. / 4) { V[i] = 4 * (-1. + i * h); }
        else V[i] = 1.;
    }
   // if (Mh / Ntau == 20) { Ntau = 2; }
   // if (Mh == 2000) { Ntau = 2; }

    for (int j = 1; j < Ntau + 1; j++)
    {
        Vup[0] = 0.;
        Vup[Mh] = 1.;

        for (int i = 1; i < Mh; i++)
        {

            Vup[i] = V[i] + (nu / 20) * (pow(V[i + 1], 2) - pow(V[i], 2));
        }

        equalize_vector(V, Vup, Mh + 1);
    }

    double Delta_V_Ch = 0.;
    double V_Ch = 0.;

    double Delta_V_L1h = 0.;
    double V_L1h = 0.;

    for (int i = 0; i < Mh + 1; i++)
    {
        //printf("V[%d]= %.4f   Sol= %.4f\n", i, V[i], nonlinear_solution(1., -1. + i * h));
        if (fabs(Vup[i] - nonlinear_solution(1., -1. + i * h)) > Delta_V_Ch) //&& fabs(Vup[i] - nonlinear_solution(1., -1. + i * h)) < 1)
        {
            Delta_V_Ch = fabs(Vup[i] - nonlinear_solution(1., -1. + i * h));
        }
        if (fabs(Vup[i]) > V_Ch)
        {
            V_Ch = fabs(Vup[i]);
        }
        Delta_V_L1h = Delta_V_L1h + fabs(Vup[i] - nonlinear_solution(1., -1. + i * h));
        V_L1h = V_L1h + fabs(Vup[i]);
    }
    Delta_V_Ch = Delta_V_Ch * 0.6;
    Delta_V_L1h = Delta_V_L1h * h;
    V_L1h = V_L1h * h;

    for (int i = 0; i < Mh + 1; i++) { Vup[i] = -1 + i * h; }
    write_2vectors_file(file, Vup, V, Mh + 1);

    fprintf(file_table, " %.3f & %.3f & %e & %e & %e & %e  \n", tau, h, Delta_V_Ch, Delta_V_L1h, Delta_V_Ch / V_Ch, Delta_V_L1h / V_L1h);
    return 0;
}*/



int calculate_nonlinear_explicit_write_file_local(FILE* file_table, int Ntau, int Mh, int K)
{
    double h = 2. / Mh;
    double tau = 1. / Ntau;
    double nu = tau / h;

    for (int i = 0; i < Mh + 1; i++)
    {
        if (-1. + i * h <= 0) { V[i] = 0.; }
        else if (-1. + i * h <= 1. / 4) { V[i] = 4 * (-1. + i * h); }
        else V[i] = 1.;
    }

    for (int j = 1; j < Ntau + 1; j++)
    {
        Vup[0] = 0.;
        Vup[Mh] = 1.;

        for (int i = 1; i < Mh; i++)
        {

            Vup[i] = V[i] + nu / 5 * (pow(V[i + 1], 2) - pow(V[i], 2));
        }

        equalize_vector(V, Vup, Mh + 1);
    }

    int Mh_tmp = Mh;
    int Ntau_tmp = Ntau;
    int n = 1;

    for (int k = 0; k < K; k++)
    {

        Mh_tmp = Mh_tmp * 2;
        Ntau_tmp = Ntau_tmp * 2;
        n = n * 2;

        double h_tmp = 2. / Mh_tmp;
        double tau_tmp = 1. / Ntau_tmp;
        double nu_tmp = tau_tmp / h_tmp;

        for (int i = 0; i < Mh_tmp + 1; i++)
        {
            if (-1. + i * h_tmp <= 0) { V_tmp[i] = 0.; }
            else if (-1. + i * h_tmp <= 1. / 4) { V_tmp[i] = 4 * (-1. + i * h_tmp); }
            else V_tmp[i] = 1.;
        }

        for (int j = 1; j < Ntau_tmp + 1; j++)
        {
            Vup_tmp[0] = 0.;
            Vup_tmp[Mh_tmp] = 1.;

            for (int i = 1; i < Mh_tmp; i++)
            {

                Vup_tmp[i] = V_tmp[i] + nu_tmp / 2 * (pow(V_tmp[i + 1], 2) - pow(V_tmp[i], 2));
            }

            equalize_vector(V_tmp, Vup_tmp, Mh_tmp + 1);
        }


        double Delta_V_Ch = 0.;
        double V_Ch = 0.;

        double Delta_V_L1h = 0.;
        double V_L1h = 0.;

        for (int i = 0; i < Mh + 1; i++)
        {
            if (fabs(V[i] - V_tmp[i * n]) > Delta_V_Ch)
            {
                Delta_V_Ch = fabs(V[i] - V_tmp[i * n]);
            }
            if (fabs(V[i]) > V_Ch)
            {
                V_Ch = fabs(V[i]);
            }
            Delta_V_L1h = Delta_V_L1h + fabs(V[i] - V_tmp[i * n]);
            V_L1h = V_L1h + fabs(V[i]);
        }
        Delta_V_L1h = Delta_V_L1h * h;
        V_L1h = V_L1h * h;
        for (int i = 0; i < Mh_tmp + 1; i++) { V[i] = -1 + i * h_tmp; }
        write_2vectors_file(file_table, V, V_tmp, Mh_tmp + 1);
        fprintf(file_table, " %f & %f & %e & %e & %e & %e \n", tau_tmp, h_tmp, Delta_V_Ch, Delta_V_L1h, Delta_V_Ch / V_Ch, Delta_V_L1h / V_L1h);
    }
    return 0;
}


double Vup_cr[1000000];
double Vup_nx[1000000];
double alpha[1000000];
double beta[1000000];

double F(double v) { return -v * v / 2; }
double F_d(double v) { return -v; }

int gn_sweep(int N, double* V_right, double* Vup_g)
{
    double* A = new double[N];
    double* B = new double[N];


    A[0] = 0.;
    A[1] = -beta[0];

    B[0] = 0;
    B[1] = V_right[0];

    for (int i = 2; i < N; i++)
    {
        A[i] = -beta[i - 1] / (alpha[i - 1] * A[i - 1] + 1.);
        B[i] = (V_right[i - 1] - alpha[i - 1] * B[i - 1]) / (alpha[i - 1] * A[i - 1] + 1.);
    }

    Vup_g[N - 1] = (V_right[N - 1] - B[N - 1] * alpha[N - 1]) / (A[N - 1] * alpha[N - 1] + 1.);

    for (int i = N - 2; i >= 0; i--)
    {
        Vup_g[i] = A[i + 1] * Vup_g[i + 1] + B[i + 1];
    }
    delete[] A;
    delete[] B;
    return 0;
}

int calculate_nonlinear_implicit_write_file(FILE* file, FILE* file_table, int Ntau, int Mh, double tol, double w)
{
    double tau = 1. / Ntau;
    double h = 2. / Mh;
    double nu = tau / h;
    double Tol = 10;
    double Rol = 0;
    for (int i = 0; i < Mh + 1; i++)
    {
        if (-1. + i * h <= 0) { V[i] = 0.; }
        else if (-1. + i * h <= 1. / 4) { V[i] = 4 * (-1. + i * h); }
        else V[i] = 1.;
    }

    //метод Ньютона первый шаг
    int k = 0;
    for (int j = 1; j < Ntau + 1; j++)
    {
        Vup[0] = 0;
        Vup[Mh] = 1;

        k = k + 1;
        {
            equalize_vector(Vup_cr, Vup_nx, Mh - 1);
            Rol = 0;

            V_right[0] = Vup_cr[0] + nu / 2. * F(Vup_cr[1]) - V[1] - nu / 2. * F(V[0]) - w / tau * (V[0] - 2 * V[1] + V[2]);
            for (int i = 1; i < Mh - 2; i++)
            {
                V_right[i] = Vup_cr[i] - nu / 2 * F(Vup_cr[i - 1]) - V[i + 1] + nu / 2 * F(Vup_cr[i + 1]) - w / tau * (V[i] - 2 * V[i + 1] + V[i + 2]);
            }
            V_right[Mh - 2] = Vup_cr[Mh - 2] - nu / 2 * F(Vup_cr[Mh - 3]) - V[Mh - 1] + nu / 2 * F(V[Mh]) - w / tau * (V[Mh - 2] - 2 * V[Mh - 1] + V[Mh]);
            for (int i = 0; i < Mh - 1; i++)V_right[i] = -V_right[i];

            alpha[0] = 0;
            for (int i = 0; i < Mh - 1; i++) { alpha[i] = -nu / 2 * F_d(Vup_cr[i - 1]); }

            for (int i = 0; i < Mh - 2; i++) { beta[i] = nu / 2 * F_d(Vup_cr[i + 1]); }
            beta[Mh - 2] = 0;

            gn_sweep(Mh - 1, V_right, Vup_nx);
            for (int i = 0; i < Mh - 1; i++) { if (fabs(Vup_nx[i]) > Rol)Rol = fabs(Vup_nx[i]); }
            //for (int i = 0; i < Mh - 1; i++) { if (fabs(Vup_nx[i]) > Tol)Tol = fabs(Vup_nx[i]); }
            for (int i = 0; i < Mh - 1; i++)Vup_nx[i] = Vup_nx[i] + Vup_cr[i];
            if (fabs(Rol) > 0) {
                Tol = Rol;
            }

        }
        //конец первого шага
        equalize_vector(Vup_nx, V + 1, Mh - 1);
        // метод Ньютона начало цикла
        while (Tol > 0.00000001) {
            k = k + 1;
            {
                equalize_vector(Vup_cr, Vup_nx, Mh - 1);
                Rol = 0;

                V_right[0] = Vup_cr[0] + nu / 2. * F(Vup_cr[1]) - V[1] - nu / 2. * F(V[0]) - w / tau * (V[0] - 2 * V[1] + V[2]);
                for (int i = 1; i < Mh - 2; i++)
                {
                    V_right[i] = Vup_cr[i] - nu / 2 * F(Vup_cr[i - 1]) - V[i + 1] + nu / 2 * F(Vup_cr[i + 1]) - w / tau * (V[i] - 2 * V[i + 1] + V[i + 2]);
                }
                V_right[Mh - 2] = Vup_cr[Mh - 2] - nu / 2 * F(Vup_cr[Mh - 3]) - V[Mh - 1] + nu / 2 * F(V[Mh]) - w / tau * (V[Mh - 2] - 2 * V[Mh - 1] + V[Mh]);
                for (int i = 0; i < Mh - 1; i++)V_right[i] = -V_right[i];

                alpha[0] = 0;
                for (int i = 0; i < Mh - 1; i++) { alpha[i] = -nu / 2 * F_d(Vup_cr[i - 1]); }

                for (int i = 0; i < Mh - 2; i++) { beta[i] = nu / 2 * F_d(Vup_cr[i + 1]); }
                beta[Mh - 2] = 0;

                gn_sweep(Mh - 1, V_right, Vup_nx);
                for (int i = 0; i < Mh - 1; i++) { if (fabs(Vup_nx[i]) > Rol)Rol = fabs(Vup_nx[i]); }
                for (int i = 0; i < Mh - 1; i++)Vup_nx[i] = Vup_nx[i] + Vup_cr[i];
                if (fabs(Rol) > 0) {
                    Tol = Rol;
                }



            }

        }

        equalize_vector(Vup + 1, Vup_cr, Mh - 1);
        equalize_vector(V, Vup, Mh + 1);
    }
    printf("final K = %d  _____LLLLLLLLLLLLLLLLLLLLLLLLLLL\n", k);
    printf("final Tol = %.9f\n", Tol);
    double Delta_V_Ch = 0.;
    double V_Ch = 0.;

    double Delta_V_L1h = 0.;
    double V_L1h = 0.;

   for (int i = 0; i < Mh + 1; i++)
    { 
        //printf("V[%d]= %.4f   Sol= %.4f\n", i, V[i], nonlinear_solution(1., -1. + i * h));
        if (fabs(V[i] - nonlinear_solution(1., -1. + i * h)) > Delta_V_Ch)
        {
            Delta_V_Ch = fabs(V[i] - nonlinear_solution(1., -1. + i * h));
        }
        if (fabs(V[i]) > V_Ch)
        {
            V_Ch = fabs(V[i]);
        }
        Delta_V_L1h = Delta_V_L1h + fabs(V[i] - nonlinear_solution(1., -1. + i * h));
        V_L1h = V_L1h + fabs(V[i]);
    }
    Delta_V_L1h = Delta_V_L1h * h;
    V_L1h = V_L1h * h;

    for (int i = 0; i < Mh + 1; i++) { Vup[i] = -1 + i * h; }

    write_2vectors_file(file, Vup,V, Mh + 1);

    fprintf(file_table, "//hline %.3f & %.3f & %e & %e & %e & %e ///\n", tau, h, Delta_V_Ch, Delta_V_L1h, Delta_V_Ch / V_Ch, Delta_V_L1h / V_L1h);
    return 0;
}

int acalculate_nonlinear_implicit_write_file(FILE* file, FILE* file_table, int Ntau, int Mh, double tol, double w)
{
    double tau = 1. / Ntau;
    double h = 2. / Mh;
    double nu = tau / h;
    double Tol = 10;
    double Rol = 0;
    for (int i = 0; i < Mh + 1; i++)
    {
        if (-1. + i * h <= 0) { V[i] = 0.; }
        else if (-1. + i * h <= 1. / 4) { V[i] = 4 * (-1. + i * h); }
        else V[i] = 1.;
    }

    //метод Ньютона первый шаг
    int k = 0;
    for (int j = 1; j < Ntau + 1; j++)
    {
        Vup[0] = 0;
        Vup[Mh] = 1;

        k = k + 1;
        {
            equalize_vector(Vup_cr, Vup_nx, Mh - 1);
            Rol = 0;

            V_right[0] = Vup_cr[0] + nu / 2. * F(Vup_cr[1]) - V[1] - nu / 2. * F(V[0]) - w / tau * (V[0] - 2 * V[1] + V[2]);
            for (int i = 1; i < Mh - 2; i++)
            {
                V_right[i] = Vup_cr[i] - nu / 2 * F(Vup_cr[i - 1]) - V[i + 1] + nu / 2 * F(Vup_cr[i + 1]) - w / tau * (V[i] - 2 * V[i + 1] + V[i + 2]);
            }
            V_right[Mh - 2] = Vup_cr[Mh - 2] - nu / 2 * F(Vup_cr[Mh - 3]) - V[Mh - 1] + nu / 2 * F(V[Mh]) - w / tau * (V[Mh - 2] - 2 * V[Mh - 1] + V[Mh]);
            for (int i = 0; i < Mh - 1; i++)V_right[i] = -V_right[i];

            alpha[0] = 0;
            for (int i = 0; i < Mh - 1; i++) { alpha[i] = -nu / 2 * F_d(Vup_cr[i - 1]); }

            for (int i = 0; i < Mh - 2; i++) { beta[i] = nu / 2 * F_d(Vup_cr[i + 1]); }
            beta[Mh - 2] = 0;

            gn_sweep(Mh - 1, V_right, Vup_nx);
            for (int i = 0; i < Mh - 1; i++) { if (fabs(Vup_nx[i]) > Rol)Rol = fabs(Vup_nx[i]); }
            //for (int i = 0; i < Mh - 1; i++) { if (fabs(Vup_nx[i]) > Tol)Tol = fabs(Vup_nx[i]); }
            for (int i = 0; i < Mh - 1; i++)Vup_nx[i] = Vup_nx[i] + Vup_cr[i];
            if (fabs(Rol) > 0) {
                Tol = Rol;
            }

        }
        //printf("final K = %d Rol = %.9f  _____LLLLLLLLLLLLLLLLLLLLLLLLLLL\n", k, Rol);
        //конец первого шага
        equalize_vector(Vup_nx, V + 1, Mh - 1);
        // метод Ньютона начало цикла
        while (k < 1000002) {
            k = k + 1;
            {
                equalize_vector(Vup_cr, Vup_nx, Mh - 1);
                Rol = 0;

                V_right[0] = Vup_cr[0] + nu / 2. * F(Vup_cr[1]) - V[1] - nu / 2. * F(V[0]) - w / tau * (V[0] - 2 * V[1] + V[2]);
                for (int i = 1; i < Mh - 2; i++)
                {
                    V_right[i] = Vup_cr[i] - nu / 2 * F(Vup_cr[i - 1]) - V[i + 1] + nu / 2 * F(Vup_cr[i + 1]) - w / tau * (V[i] - 2 * V[i + 1] + V[i + 2]);
                }
                V_right[Mh - 2] = Vup_cr[Mh - 2] - nu / 2 * F(Vup_cr[Mh - 3]) - V[Mh - 1] + nu / 2 * F(V[Mh]) - w / tau * (V[Mh - 2] - 2 * V[Mh - 1] + V[Mh]);
                for (int i = 0; i < Mh - 1; i++)V_right[i] = -V_right[i];

                alpha[0] = 0;
                for (int i = 0; i < Mh - 1; i++) { alpha[i] = -nu / 2 * F_d(Vup_cr[i - 1]); }

                for (int i = 0; i < Mh - 2; i++) { beta[i] = nu / 2 * F_d(Vup_cr[i + 1]); }
                beta[Mh - 2] = 0;

                gn_sweep(Mh - 1, V_right, Vup_nx);
                for (int i = 0; i < Mh - 1; i++) { if (fabs(Vup_nx[i]) > Rol)Rol = fabs(Vup_nx[i]); }
                for (int i = 0; i < Mh - 1; i++)Vup_nx[i] = Vup_nx[i] + Vup_cr[i];
                if (fabs(Rol) > 0) {
                    Tol = Rol;
                }



            }

        }

        equalize_vector(Vup + 1, Vup_cr, Mh - 1);
        equalize_vector(V, Vup, Mh + 1);
    }
    //printf("final K = %d  _____LLLLLLLLLLLLLLLLLLLLLLLLLLL\n", k);
    printf("final Tol = %.9f\n", Tol);
    double Delta_V_Ch = 0.;
    double V_Ch = 0.;

    double Delta_V_L1h = 0.;
    double V_L1h = 0.;

    for (int i = 0; i < Mh + 1; i++)
    {
        //printf("V[%d]= %.4f   Sol= %.4f\n", i, V[i], nonlinear_solution(1., -1. + i * h));
        if (fabs(V[i] - nonlinear_solution(1., -1. + i * h)) > Delta_V_Ch)
        {
            Delta_V_Ch = fabs(V[i] - nonlinear_solution(1., -1. + i * h));
        }
        if (fabs(V[i]) > V_Ch)
        {
            V_Ch = fabs(V[i]);
        }
        Delta_V_L1h = Delta_V_L1h + fabs(V[i] - nonlinear_solution(1., -1. + i * h));
        V_L1h = V_L1h + fabs(V[i]);
    }
    Delta_V_L1h = Delta_V_L1h * h;
    V_L1h = V_L1h * h;
    for (int i = 0; i < Mh + 1; i++) { Vup[i] = -1 + i * h; }
    for (int i = 0; i < 1+ Mh/10; i++)
    {
        fprintf(file, "%f  %f  \n", Vup[i*10], V[i*10]);
    

}
    //
    //write_2vectors_file(file, Vup, V, Mh + 1);
    fprintf(file_table, "//hline %.3f & %.3f & %e & %e & %e & %e ///\n", tau, h, Delta_V_Ch, Delta_V_L1h, Delta_V_Ch / V_Ch, Delta_V_L1h / V_L1h);
    return 0;
}



int main(void) {

   /* FILE* file1;
    fopen_s(&file1, "t1.txt", "w");
    FILE* file2;
    fopen_s(&file2, "t1_t.txt", "w");
    FILE* file7;
    fopen_s(&file7, "t4_t.txt", "w");
    FILE* file8;
    fopen_s(&file8, "t5_t.txt", "w");
    FILE* file3;
    fopen_s(&file3, "t2.txt", "w");
    FILE* file4;
    fopen_s(&file4, "t2_t.txt", "w");
    FILE* file9;
    fopen_s(&file9, "t6_t.txt", "w");
    FILE* file10;
    fopen_s(&file10, "t7_t.txt", "w");
    FILE* file5;
    fopen_s(&file5, "t3.txt", "w");
    FILE* file6;
    fopen_s(&file6, "t3_t.txt", "w");
    FILE* file11;
    fopen_s(&file11, "t8_t.txt", "w");
    FILE* file12;
    fopen_s(&file12, "t9_t.txt", "w");


    FILE* File1;
    fopen_s(&File1, "t10.txt", "w");
    FILE* File2;
    fopen_s(&File2, "t10_t.txt", "w");
    FILE* File3;
    fopen_s(&File3, "t11_t.txt", "w");
    FILE* File4;
    fopen_s(&File4, "t12_t.txt", "w");
    FILE* File5;
    fopen_s(&File5, "t13.txt", "w");
    FILE* File6;
    fopen_s(&File6, "t13_t.txt", "w");*/
    FILE* File7;
    fopen_s(&File7, "t14.txt", "w");
    FILE* File8;
    fopen_s(&File8, "t14_t.txt", "w");
   /* FILE* File10;
    fopen_s(&File10, "t15_t.txt", "w");
    FILE* File9;
    fopen_s(&File9, "t15.txt", "w");
  
    // Расчет линейной задачи с явной схемой

       // Таблица 1
    calculate_linear_explicit_write_file(file1, file2, 10, 20);
    calculate_linear_explicit_write_file(file1, file2, 100, 20);
    calculate_linear_explicit_write_file(file1, file2, 1000, 20);
    calculate_linear_explicit_write_file(file1, file2, 10, 200);
    calculate_linear_explicit_write_file(file1, file2, 100, 200);
    calculate_linear_explicit_write_file(file1, file2, 1000, 200);
    calculate_linear_explicit_write_file(file1, file2, 10, 2000);
    calculate_linear_explicit_write_file(file1, file2, 100, 2000);
    calculate_linear_explicit_write_file(file1, file2, 1000, 2000);

    // Таблица 2
    calculate_linear_explicit_write_file_local(file7, 10, 20, 4);
    // Таблица 3
    calculate_linear_explicit_write_file_local(file8, 100, 200, 4);



    // Расчет линейной задачи с неявной схемой

        // Таблица 1
    calculate_linear_implicit_write_file(file3, file4, 10, 20, 0.1);
    calculate_linear_implicit_write_file(file3, file4, 100, 20, 0.1);
    calculate_linear_implicit_write_file(file3, file4, 1000, 20, 0.1);
    calculate_linear_implicit_write_file(file3, file4, 10, 200, 0.1);
    calculate_linear_implicit_write_file(file3, file4, 100, 200, 0.1);
    calculate_linear_implicit_write_file(file3, file4, 1000, 200, 0.1);
    calculate_linear_implicit_write_file(file3, file4, 10, 2000, 0.1);
    calculate_linear_implicit_write_file(file3, file4, 100, 2000, 0.1);
    calculate_linear_implicit_write_file(file3, file4, 1000, 2000, 0.1);

    // Таблица 3
    calculate_linear_implicit_write_file_local(file9, 10, 20, 4, 0.1);
    // Таблица 4
    calculate_linear_implicit_write_file_local(file10, 100, 200, 4, 0.1);

    // Таблица 2

    calculate_linear_implicit_write_file(file5, file6, 10, 20, 1);
    calculate_linear_implicit_write_file(file5, file6, 100, 20, 1);
    calculate_linear_implicit_write_file(file5, file6, 1000, 20, 1);
    calculate_linear_implicit_write_file(file5, file6, 10, 200, 1);
    calculate_linear_implicit_write_file(file5, file6, 100, 200, 1);
    calculate_linear_implicit_write_file(file5, file6, 1000, 200, 1);
    calculate_linear_implicit_write_file(file5, file6, 10, 2000, 1);
    calculate_linear_implicit_write_file(file5, file6, 100, 2000, 1);
    calculate_linear_implicit_write_file(file5, file6, 1000, 2000, 1);

    // Таблица 5
    calculate_linear_implicit_write_file_local(file11, 10, 20, 4, 1);
    // Таблица 6
    calculate_linear_implicit_write_file_local(file12, 100, 200, 4, 1);

    // Расчет нелинейной задачи с явной схемой

    // Таблица 1
    calculate_nonlinear_explicit_write_file(File1, File2, 10, 20);
    calculate_nonlinear_explicit_write_file(File1, File2, 100, 20);
    calculate_nonlinear_explicit_write_file(File1, File2, 1000, 20);
    calculate_nonlinear_explicit_write_file(File1, File2, 10, 200);
    calculate_nonlinear_explicit_write_file(File1, File2, 100, 200);
    calculate_nonlinear_explicit_write_file(File1, File2, 1000, 200);
    calculate_nonlinear_explicit_write_file(File1, File2, 10, 2000);
    calculate_nonlinear_explicit_write_file(File1, File2, 100, 2000);
    calculate_nonlinear_explicit_write_file(File1, File2, 1000, 2000);
    // Таблица 2
    calculate_nonlinear_explicit_write_file_local(File3, 10, 20, 4);
    // Таблица 3
    calculate_nonlinear_explicit_write_file_local(File4, 100, 200, 4);*/


    // // Расчет нелинейной задачи с неявной схемой*/

    // Таблица 3
   /* acalculate_nonlinear_implicit_write_file(File9, File10, 10, 20, pow(10, -9), 0.1);
    acalculate_nonlinear_implicit_write_file(File9, File10, 100, 20, pow(10, -9), 0.1);
    acalculate_nonlinear_implicit_write_file(File9, File10, 1000, 20, pow(10, -9), 0.1);
    acalculate_nonlinear_implicit_write_file(File9, File10, 10, 200, pow(10, -9), 0.1);
    acalculate_nonlinear_implicit_write_file(File9, File10, 100, 200, pow(10, -9), 0.1);
   acalculate_nonlinear_implicit_write_file(File9, File10, 1000, 200, pow(10, -9), 0.1);
    acalculate_nonlinear_implicit_write_file(File9, File10, 10, 2000, pow(10, -9), 0.1);
    acalculate_nonlinear_implicit_write_file(File9, File10, 100, 2000, pow(10, -9), 0.1);
    acalculate_nonlinear_implicit_write_file(File9, File10, 1000, 2000, pow(10, -9), 0.1);

    // Таблица 1
    calculate_nonlinear_implicit_write_file(File5, File6, 10, 2000, pow(10, -9), 0.1);
    calculate_nonlinear_implicit_write_file(File5, File6, 10, 20000, pow(10, -9), 0.1);
    calculate_nonlinear_implicit_write_file(File5, File6, 5, 2000, pow(10, -9), 0.1);
    calculate_nonlinear_implicit_write_file(File5, File6, 10, 10000, pow(10, -9), 0.1);
    calculate_nonlinear_implicit_write_file(File5, File6, 5, 20000, pow(10, -9), 0.1);
    calculate_nonlinear_implicit_write_file(File5, File6, 5, 1000, pow(10, -9), 0.1);*/
    // Таблица 2
    acalculate_nonlinear_implicit_write_file(File7, File8, 10, 2000, pow(10, -9), 1);
    acalculate_nonlinear_implicit_write_file(File7, File8, 10, 20000, pow(10, -9),1);
    acalculate_nonlinear_implicit_write_file(File7, File8, 5, 2000, pow(10, -9), 1);
    acalculate_nonlinear_implicit_write_file(File7, File8, 10, 10000, pow(10, -9), 1);
    acalculate_nonlinear_implicit_write_file(File7, File8, 5, 20000, pow(10, -9), 1);
    acalculate_nonlinear_implicit_write_file(File7, File8, 5, 1000, pow(10, -9),1);
    
    return 0;
}