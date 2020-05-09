//
// Created by Khaled Maksoud on 08/04/2019.
// Example for a Velocity Verlet integration scheme on a single 1D harmonic oscillator

#include <cmath>
#include <vector>
#include <iostream>
#define MDSteps 10000000

using namespace std;

// Set initial functions for calculating energies + forces



double calc_force(double x, const int k)
{
    return -k*x;
}

double calc_potential(double x, const int k)
{
    return 0.5*k*pow(x, 2);
}

double calc_kinetic(double v, const int mass)
{
    return 0.5*mass*pow(v, 2);
}

//Set functions for Velocity Verlet update of positions and velocities

double position_update(double x, double v, const float dt, const double stepfraction=1.0)
{
    return x + v*dt*stepfraction;
}

double velocity_update(double v, double f, const int mass, const float dt, const double stepfraction=1.0)
{
    return v + (0.5*dt*stepfraction*f)/mass;
}

struct MDarray{
    double * q_traj;
    double * p_traj;
    double * e_pot;
    double * e_kin;
    double * e_tot;
};



MDarray VelocityVerletIntegrator(double x, double v, const float dt=0.1, const int mass=1, const int k=10)
{

    const int length = MDSteps;

    static double  x_trj[length];
    static double  v_trj[length];
    static double  Epot[length];
    static double  Ekin[length];
    static double  Etot[length];

    double f;

    for (int step = 0; step <= MDSteps; step++)
    {
        f = calc_force(x, k);
        v = velocity_update(v, f, mass, dt, 0.5);
        x = position_update(x, v, dt, 1.0);
        f = calc_force(x, k);
        v = velocity_update(v, f, mass, dt, 0.5);

        x_trj[step] = x;
        v_trj[step] = v;

        double U = calc_potential(x, k);
        Epot[step] = U;

        double K = calc_kinetic(v, mass);
        Ekin[step] = K;

        double Tot = U + K;
        Etot[step] = Tot;

    }

    return MDarray{x_trj, v_trj, Epot, Ekin, Etot};
}

int main(int argc, const char **argv)
{

    //Set initial velocity and position of particle

    double v_0 = 1.0;
    double x_0 = 1.0;

    int numsteps = MDSteps;

    MDarray T;

    T = VelocityVerletIntegrator(x_0, v_0, 0.1, 1, 10);

    for (int i = 0; i <= numsteps; i++)
    {
        cout << "Step " << i << "/" << numsteps << " - " << numsteps-i << " steps remaining. " << endl;
        cout << T.q_traj[i] << "\t" << T.p_traj[i] << "\t" << T.e_pot[i] << "\t" << T.e_kin[i] << endl;
        cout << "Total Energy = " << T.e_tot[i] << "\n";
    }


    return 0;
}
