//
// Created by Khaled Maksoud on 2019-04-08.
// Example for a Velocity Verlet integration scheme on a single 3D harmonic oscillator

#include <cmath>
#include <vector>
#include <iostream>
#include <random>
#define MDSteps 1000000

using namespace std;

// Set initial functions for calculating energies + forces

const int n_atoms = 2;

// Simulation temperature
const double temperature = 300;   // kelvin

const double k_boltz = 1.987206504191549E-003;  // kcal mol-1 K-1

double box_size[3] = {2.0, 2.0, 2.0};

const double r_eq  = 1.54; //Angstroms

// Subroutine to apply periodic boundaries
double make_periodic(double x, const double box)
{
    while (x < -0.5*box)
    {
        x = x + box;
    }

    while (x > 0.5*box)
    {
        x = x - box;
    }

    return x;
}

// Subroutine to wrap the coordinates into a box
double wrap_into_box(double x, double box)
{
    while (x > box)
    {
        x = x - box;
    }

    while (x < 0)
    {
        x = x + box;
    }

    return x;
}


// Subroutine to print a PDB of the coordinates
void print_pdb(double **coords, const int n_atoms, const int step)
{
    char filename[64];

    //snprintf(filename, 64, "output%d.pdb", step);

    snprintf(filename, 64, "harmonic_%00008d.pdb", step);

    FILE *f = fopen(filename, "w");

    fprintf(f, "CRYST1 %8.3f %8.3f %8.3f  90.00  90.00  90.00\n",
            box_size[0], box_size[1], box_size[2]);

    for (int i = 0; i < n_atoms; i = i + 1)
    {

//        coords[i][0] = wrap_into_box(coords[i][0], box_size[0]);
//        coords[i][1] = wrap_into_box(coords[i][1], box_size[1]);
//        coords[i][2] = wrap_into_box(coords[i][2], box_size[2]);

        fprintf(f, "ATOM  %5d  C   C     1    %8.3f%8.3f%8.3f  1.00  0.00          C\n",
                i+1, coords[i][0], coords[i][1], coords[i][2]);
        fprintf(f, "TER\n");
    }

    fclose(f);
}


// Function to assign a random set of velocities drawn from the kinetic energy distribution
double **assign_velocities(const double T, const double Kb)
{
    auto **vels = new double *[n_atoms];

    default_random_engine gen_random;
    normal_distribution<double> P_v(0,Kb*T);


    for (int i = 0; i < n_atoms; i++)
    {
        vels[i] = new double[3];

        vels[i][0] = P_v(gen_random);
        vels[i][1] = P_v(gen_random);
        vels[i][2] = P_v(gen_random);
    }

    return vels;

}



// function to return a random number between 'start' to 'end'
double rand(const double start, const double end)
{
    return (end-start) * (double(rand()) / RAND_MAX) + start;
}

// Calculate r_ij for a pair of atomic coordinates
const double get_radial(double **coords, const double *box_size, int i, int j) //Will be useful only in LJ potential example
{

    double delta_x = coords[j][0] - coords[i][0];
    double delta_y = coords[j][1] - coords[i][1];
    double delta_z = coords[j][2] - coords[i][2];

    // Apply periodic boundaries
    delta_x = make_periodic(delta_x, box_size[0]);
    delta_y = make_periodic(delta_y, box_size[1]);
    delta_z = make_periodic(delta_z, box_size[2]);

    const double r_ij = (delta_x * delta_x) + (delta_y * delta_y) +
                        (delta_z * delta_z);

    return r_ij;
}


double **calc_3Dforce(double **pos, const double k)
{
    auto **f = new double *[n_atoms];

    for (int i = 0; i < n_atoms; i++) {
        f[i] = new double[3];

        f[i][0] = -k * pos[i][0];
        f[i][1] = -k * pos[i][1];
        f[i][2] = -k * pos[i][2];
    }

    for (int i = 0; i < n_atoms-1; i = i + 1)
    {
        for (int j = i+1; j < n_atoms; j = j + 1)
        {

            double delta_x = pos[i][0] - pos[j][0];
            double delta_y = pos[i][1] - pos[j][1];
            double delta_z = pos[i][2] - pos[j][2];

            // Apply periodic boundaries
            delta_x = make_periodic(delta_x, box_size[0]);
            delta_y = make_periodic(delta_y, box_size[1]);
            delta_z = make_periodic(delta_z, box_size[2]);

            const double r2 = (delta_x*delta_x) + (delta_y*delta_y) +
                              (delta_z*delta_z);



            // F_LJ = 48*epsilon*1/r2[ (sigma/r)^12 - 0.5*(sigma/r)^6 ]
            const double sig2_over_r2 = 1.0 / r2;
            const double sig6_over_r6 = sig2_over_r2*sig2_over_r2*sig2_over_r2;
            //const double sig12_over_r12 = sig6_over_r6 * sig6_over_r6;

            const double f_x = 48.0 * sig2_over_r2 * sig6_over_r6 * ( sig6_over_r6 - 0.5 );

            cout << "force = " << f_x << "\t" << "radial comp = " << r2
            << "\t" << delta_x << "\t" << delta_y << "\t" << delta_z << endl;

            f[i][0] = f[i][0] + (f_x * delta_x);
            f[i][1] = f[i][1] + (f_x * delta_y);
            f[i][2] = f[i][2] + (f_x * delta_z); //update force components for atom i

            f[j][0] = f[j][0] - (f_x * delta_x);
            f[j][1] = f[j][1] - (f_x * delta_y);
            f[j][2] = f[j][2] - (f_x * delta_z); //update force components for atom j



        }
    }

    return f;

}

double **calc_bondforce(double **pos, const double kr, const double r_eq)
{

    auto **f = new double *[n_atoms];

    for (int t = 0; t < n_atoms; t++){

        f[t] = new double[3];

        f[t][0] = 0.0;
        f[t][1] = 0.0;
        f[t][2] = 0.0;

    }

    for (int i = 0; i < n_atoms-1; i++) {
        for (int j = i+1; j < n_atoms; j++){

            double delta_x = pos[i][0] - pos[j][0];
            double delta_y = pos[i][1] - pos[j][1];
            double delta_z = pos[i][2] - pos[j][2];

            // Apply periodic boundaries
            delta_x = make_periodic(delta_x, box_size[0]);
            delta_y = make_periodic(delta_y, box_size[1]);
            delta_z = make_periodic(delta_z, box_size[2]);

            const double r = sqrt((delta_x*delta_x) + (delta_y*delta_y) +
                              (delta_z*delta_z));

            const double fscalar = -2.0 * kr * (1 - (r_eq / r));

//            cout << "force = " << fscalar << "\t" << "radial comp = " << r
//            << "\t" << delta_x << "\t" << delta_y << "\t" << delta_z << endl;

            f[i][0] = f[i][0] + (fscalar * delta_x);
            f[i][1] = f[i][1] + (fscalar * delta_y);
            f[i][2] = f[i][2] + (fscalar * delta_z); //update force components for atom i

            f[j][0] = f[j][0] - (fscalar * delta_x);
            f[j][1] = f[j][1] - (fscalar * delta_y);
            f[j][2] = f[j][2] - (fscalar * delta_z); //update force components for atom j




        }

    }


    return f;
}

double calc_potential(double **pos, const double k)
{
    double V = 0.0;

    for (int i = 0; i < n_atoms; i++)
    {
        V += 0.5*k*pow(pos[i][0], 2) + 0.5*k*pow(pos[i][1], 2) + 0.5*k*pow(pos[i][2], 2);
    }

    return V;
}

double calc_bondpotential(double **pos, const double kr, const double r_eq){

    double V = 0.0;

    for (int i = 0; i < n_atoms-1; i++) {
        for (int j = i + 1; j < n_atoms; j++) {

            double delta_x = pos[i][0] - pos[j][0];
            double delta_y = pos[i][1] - pos[j][1];
            double delta_z = pos[i][2] - pos[j][2];

            // Apply periodic boundaries
            delta_x = make_periodic(delta_x, box_size[0]);
            delta_y = make_periodic(delta_y, box_size[1]);
            delta_z = make_periodic(delta_z, box_size[2]);

            const double r = sqrt((delta_x * delta_x) + (delta_y * delta_y) +
                                  (delta_z * delta_z));

            const double E_harmonic = kr * ((r - r_eq) * (r - r_eq));

            V += E_harmonic;

        }

    }

    return V;

}

// Subroutine that calculates the potential energies of the atoms
double calc_LJpot(double **coords, const double *box_size)
{
    // Loop over all pairs of atoms and calculate
    // the LJ energy
    double pot_energy = 0;

    for (int i = 0; i < n_atoms-1; i = i + 1)
    {
        for (int j = i+1; j < n_atoms; j = j + 1)
        {

            double r2 = get_radial(coords, box_size, i, j);

            // E_LJ = 4*epsilon[ (sigma/r)^12 - (sigma/r)^6 ]
            const double sig2_over_r2 = 1 / r2;
            const double sig6_over_r6 = sig2_over_r2*sig2_over_r2*sig2_over_r2;
            const double sig12_over_r12 = sig6_over_r6 * sig6_over_r6;

            const double e_lj = 4.0  * ( sig12_over_r12 - sig6_over_r6 );

            pot_energy = pot_energy + e_lj;
        }
    }

    // return the total energy of the atoms
    return pot_energy;
}

double calc_kinetic(double **pos, const int mass)
{
    double K = 0.0;

    for (int i = 0; i < n_atoms; i++)
    {
        K += 0.5*mass*pow(pos[i][0], 2) + 0.5*mass*pow(pos[i][1], 2) + 0.5*mass*pow(pos[i][2], 2);
    }

    return K;
}

//Set functions for Velocity Verlet update of positions and velocities

void position_3Dupdate(double **pos, double **vel, double **F, const int mass, const float dt, const double stepfraction=1.0)
{
    for (int i = 0; i < n_atoms; i++)
    {
        pos[i][0] = pos[i][0] + vel[i][0]*dt*stepfraction + (0.5*dt*dt/mass)*F[i][0];
        pos[i][1] = pos[i][1] + vel[i][1]*dt*stepfraction + (0.5*dt*dt/mass)*F[i][0];
        pos[i][2] = pos[i][2] + vel[i][2]*dt*stepfraction + (0.5*dt*dt/mass)*F[i][0];
    }

}

void velocity_3Dupdate(double **vel, double **F, const int mass, const float dt, const double stepfraction=1.0)
{


    for (int i=0; i < n_atoms; i++)
    {
        vel[i][0] = vel[i][0] + (0.5*dt*stepfraction/mass)*F[i][0];
        vel[i][1] = vel[i][1] + (0.5*dt*stepfraction/mass)*F[i][1];
        vel[i][2] = vel[i][2] + (0.5*dt*stepfraction/mass)*F[i][2];
    }
}

struct MDarray{
    double ** q_traj;
    double ** p_traj;
    double * e_pot;
    double * e_kin;
    double * e_tot;
};



MDarray VelocityVerletIntegrator(double **x, double **v, double **F, const float dt=0.1, const int mass=1, const double k=10.0)
{

    const int length = MDSteps;

    auto **x_trj = new double*[n_atoms];
    auto **v_trj = new double*[n_atoms];
    static double  Epot[length];
    static double  Ekin[length];
    static double  Etot[length];



    for (int step = 0; step <= length; step++)
    {

        //double U = calc_potential(x, k) + calc_LJpot(x, box_size);
        double U = calc_bondpotential(x, k, r_eq);
        Epot[step] = U;

        double K = calc_kinetic(v, mass);
        Ekin[step] = K;

        double Tot = U + K;
        Etot[step] = Tot;

        //F = calc_3Dforce(x, k);
        F = calc_bondforce(x, k, r_eq);
        velocity_3Dupdate(v, F, mass, dt, 0.5);
        position_3Dupdate(x, v, F, mass, dt, 1.0);



        //F = calc_3Dforce(x, k);
        F = calc_bondforce(x, k, r_eq);
        velocity_3Dupdate(v, F, mass, dt, 0.5);

        if (step % 1000 == 0)
        {
            print_pdb(x, n_atoms, step);

            for (int i = 0; i < n_atoms; i++)
            {
                cout << x[i][0] << "\t" << x[i][1] << "\t" << x[i][2] << endl;
            }

            cout << "\n" ;

            for (int i = 0; i < n_atoms; i++)
            {
                cout << v[i][0] << "\t" << v[i][1] << "\t" << v[i][2] << endl;

            }

            cout << "-----------------------------------------------" << endl;

        }



        for (int i = 0; i < n_atoms; i++){

            x_trj[i] = new double[3];
            v_trj[i] = new double[3];

            x_trj[i] = x[i];
            v_trj[i] = v[i];

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

    const double fconstant = 40.0;

    double **v_0 = assign_velocities(temperature, k_boltz);
    double **x_0 = new double *[n_atoms];

    for (int i = 0; i < n_atoms; i++)
    {
        x_0[i] = new double[3];

        x_0[i][0] = rand(0 , box_size[0]);
        x_0[i][1] = rand(0 , box_size[1]);
        x_0[i][2] = rand(0 , box_size[2]);
    }






    //double **f_0 = calc_3Dforce(x_0, fconstant);
    double **f_0 = calc_bondforce(x_0, fconstant, r_eq);

    const double kT = temperature * k_boltz;

    int numsteps = MDSteps;

    cout << "Initial coordinates and velocities:-" << endl;

    for (int i = 0; i < n_atoms; i++)
    {
        cout << x_0[i][0] << "\t" << x_0[i][1] << "\t" << x_0[i][2] << endl;
    }
    cout << "-----------------------------------------------------------------------" << endl;

    for (int i = 0; i < n_atoms; i++)
    {
        cout << v_0[i][0] << "\t" << v_0[i][1] << "\t" << v_0[i][2] << endl;
    }


    MDarray T;

    T = VelocityVerletIntegrator(x_0, v_0, f_0, 0.001, 1, fconstant);

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
        if (t % 1000 == 0)
        {
            cout << "Step " << t << "/" << numsteps << " - " << numsteps-t << " steps remaining. " << endl;

//        for (int i = 0; i < n_atoms; i++)
//        {
//            cout << T.q_traj[i][0] << "\t" << T.q_traj[i][1] << "\t" << T.q_traj[i][2] << endl;
//
//        }
//
//        cout << "\n" << endl;
//
//        for (int i = 0; i < n_atoms; i++)
//        {
//            cout << T.p_traj[i][0] << "\t" << T.p_traj[i][1] << "\t" << T.p_traj[i][2] << endl;
//
//        }


        cout << "                                                                  " << endl;
        cout << "Potential = " << T.e_pot[t] << "\t" << ", Kinetic = " << T.e_kin[t] << endl;
        cout << "Total Energy = " << T.e_tot[t] << endl;
        cout << "                                                                  " << endl;

        }

    }


    return 0;
}