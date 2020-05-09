//
// Created by Khaled Maksoud on 2019-04-08.
// Example for a Velocity Verlet integration scheme on a single 3D harmonic oscillator

#include <cmath>
#include <vector>
#include <iostream>
#include <random>
#define MDSteps 5000

using namespace std;

// Set initial functions for calculating energies + forces

// Simulation temperature
const double temperature = 300;   // kelvin

const double k_boltz = 1.987206504191549E-003;  // kcal mol-1 K-1

double box_size[3] = {30.0, 30.0, 30.0};

// Subroutine to print a PDB of the coordinates
void print_pdb(double *coords, const int n_atoms, const int step)
{
    char filename[64];


    snprintf(filename, 64, "harmonic_%00000d.pdb", step);

    FILE *f = fopen(filename, "w");

    fprintf(f, "CRYST1 %8.3f %8.3f %8.3f  90.00  90.00  90.00\n",  box_size[0], box_size[1], box_size[2]);



    fprintf(f, "ATOM  %5d  Kr   Kr     1    %8.3f%8.3f%8.3f  1.00  0.00          Kr\n",
                1, coords[0], coords[1], coords[2]);
    fprintf(f, "TER\n");


    fclose(f);
}


// Function to assign a random set of velocities drawn from the kinetic energy distribution
double *assign_velocities(const double T, const double Kb)
{
    static double vels[3];

    default_random_engine gen_random;
    normal_distribution<double> P_v(0,Kb*T);


    for (int i = 0; i < 3; i++)
    {
        vels[i] = P_v(gen_random);
    }

    return vels;

}



// function to return a random number between 'start' to 'end'
double rand(const double start, const double end)
{
    return (end-start) * (double(rand()) / RAND_MAX) + start;
}

//double get_radial(double *pos, double *vel, bool choose) //Will be useful only in LJ potential example
//{
//    double r_x = 0.0;
//    double r_v = 0.0;
//
//    for (int i = 0; i < 3; i++){
//        r_x += pow(pos[i], 2);
//        r_v += pow(vel[i], 2);
//    }
//
//    r_x = sqrt(r_x);
//    r_v = sqrt(r_v);
//
//    if (choose)
//    {
//        return r_x;
//    }
//    else if (choose == false){
//        return r_v;
//    }
//}


double * calc_3Dforce(double *pos, const int k)
{
    static double f[3];

    f[0] = -k*pos[0];
    f[1] = -k*pos[1];
    f[2] = -k*pos[2];

    return f;
}

double calc_potential(double *pos, const int k)
{
    double V = 0.0;

    for (int i = 0; i < 3; i++){
        V += 0.5*k*pow(pos[i], 2);
    }

    return V;
}

double calc_kinetic(double *pos, const int mass)
{
    double K = 0.0;

    for (int i = 0; i < 3; i++){
        K += 0.5*mass*pow(pos[i], 2);
    }

    return K;
}

//Set functions for Velocity Verlet update of positions and velocities

double * position_3Dupdate(double *pos, double *vel, const float dt, const double stepfraction=1.0)
{
    static double q[3];

    for (int i = 0; i < 3; i++){
        q[i] = pos[i] + vel[i]*dt*stepfraction;
    }

    return q;
}

double * velocity_3Dupdate(double *vel, double *F, const int mass, const float dt, const double stepfraction=1.0)
{
    static double p[3];

    for (int i=0; i < 3; i++){
        p[i] = vel[i] + (0.5*dt*stepfraction/mass)*F[i];
    }

    return p;
}

struct MDarray{
    double ** q_traj;
    double ** p_traj;
    double * e_pot;
    double * e_kin;
    double * e_tot;
};



MDarray VelocityVerletIntegrator(double *x, double *v, double *F, const float dt=0.1, const int mass=1, const int k=10)
{

    const int length = MDSteps;

    auto **x_trj = new double*[length];
    auto **v_trj = new double*[length];
    static double  Epot[length];
    static double  Ekin[length];
    static double  Etot[length];



    for (int step = 0; step <= length; step++)
    {

        double U = calc_potential(x, k);
        Epot[step] = U;

        double K = calc_kinetic(v, mass);
        Ekin[step] = K;

        double Tot = U + K;
        Etot[step] = Tot;

        F = calc_3Dforce(x, k);
        v = velocity_3Dupdate(v, F, mass, dt, 0.5);
        x = position_3Dupdate(x, v, dt, 1.0);

        if (step % 100 == 0){
            print_pdb(x, 3, step);
        }

        F = calc_3Dforce(x, k);
        v = velocity_3Dupdate(v, F, mass, dt, 0.5);

        x_trj[step] = new double[3];
        v_trj[step] = new double[3];

        for (int i = 0; i < 3; i++){

            x_trj[step][i] = x[i];
            v_trj[step][i] = v[i];

        }


    }

    MDarray Traj;

    Traj.q_traj = x_trj;
    Traj.p_traj = v_trj;
    Traj.e_pot = Epot;
    Traj.e_kin = Ekin;
    Traj.e_tot = Etot;

    return Traj;
}

int main(int argc, const char **argv)
{

    // The total number of accepted moves
    int naccept = 0;

    // The total number of rejected moves
    int nreject = 0;


    //Set initial velocity and position of particle

    double *v_0 = assign_velocities(temperature, k_boltz);
    double x_0[3] = {rand(0,2), rand(0,2), rand(0,2)};
    double f_0[3] = {0.0, 0.0, 0.0};

    const double kT = temperature * k_boltz;

    int numsteps = MDSteps;

    MDarray T;

    T = VelocityVerletIntegrator(x_0, v_0, f_0, 0.1, 1, 5);

    //Metropolis-Hastings acceptance test for Hamiltonian

    const double H_inital = T.e_tot[0];
    const double H_final = T.e_tot[-1];

    bool accept = false;

    if (H_final <= H_inital)
    {

        accept = true;

    }
    else
    {

        double q = exp( H_final - H_inital / kT);

        if (q >= rand(0.0, 1.0))
        {
            accept = true;
        }
        else
        {
            accept = false;
        }


    }

    if (accept)
    {
        //Accept the MD trajectory
        naccept += 1;
        //double *v_next = T.q_traj[-1];
    }
    else
    {

        //Reject the MD trajectory and reset the coordinates
        nreject += 1;
        //double *v_next = T.q_traj[0];

    }

    for (int t = 0; t <= numsteps; t++)
    {
        cout << "Step " << t << "/" << numsteps << " - " << numsteps-t << " steps remaining. " << endl;

        for (int i = 0; i < 3; i++)
        {
            cout << T.q_traj[t][i] << "\t" << T.p_traj[t][i] << endl;
        }


        cout << "                                                                  " << endl;
        cout << "Potential = " << T.e_pot[t] << "\t" << ", Kinetic = " << T.e_kin[t] << endl;
        cout << "Total Energy = " << T.e_tot[t] << endl;
    }


    return 0;
}