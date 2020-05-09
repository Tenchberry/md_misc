//
// Created by Khaled Maksoud on 2019-04-08.
// Example for a Velocity Verlet integration scheme on a Lennard-Jones fluid

#include <cstdio>
#include <cstdlib>
#include <cmath>
//#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#define MDSteps 50000

using namespace std;

// Set initial functions for calculating energies + forces

// Set the number of atoms in the box
const int n_atoms = 200;

// Set the size of the box (in Angstroms)
double box_size[3] = { 100.0, 100.0, 100.0 };

// Simulation temperature
const double temperature = 298;   // kelvin

const double k_boltz = 1.987206504191549E-003;  // kcal mol-1 K-1
//const double k_boltz = 1.3806485279E-023; //J K-1

// Simulation pressure (atmospheres converted to internal
// units - kcal mol-1 A-3)
double pressure = 1 * 1.458397506863647E-005;   // atmospheres

// The maximum amount to change the volume - the
// best value is 10% of the number of atoms
//const double max_volume_change = 0.1 * n_atoms;   // Angstroms**3

// Give the Lennard Jones parameters for the atoms
// (these are the OPLS parameters for Krypton)
//const double sigma = 3.624;     // angstroms
//const double epsilon = 0.317;   // kcal mol-1
const double sigma = 1.5;
const double epsilon = 0.2;
const double mass = 1.00;     // Da = factor * Kg mol-1


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


// function to return a random number between 'start' to 'end'
double rand(const double start, const double end)
{
    return (end-start) * (double(rand()) / RAND_MAX) + start;
}

// Subroutine to print a PDB of the coordinates
void print_pdb(double **coords, const int n_atoms, const int step)
{
    char filename[64];

    //snprintf(filename, 64, "output%d.pdb", step);

    snprintf(filename, 64, "LJfluid_%000008d.pdb", step);

    FILE *f = fopen(filename, "w");

    fprintf(f, "CRYST1 %8.3f %8.3f %8.3f  90.00  90.00  90.00\n",
            box_size[0], box_size[1], box_size[2]);

    for (int i = 0; i < n_atoms; i = i + 1)
    {

        coords[i][0] = wrap_into_box(coords[i][0], box_size[0]);
        coords[i][1] = wrap_into_box(coords[i][1], box_size[1]);
        coords[i][2] = wrap_into_box(coords[i][2], box_size[2]);

        fprintf(f, "ATOM  %5d  C   C     1    %8.3f%8.3f%8.3f  1.00  0.00          Kr\n",
                i+1, coords[i][0], coords[i][1], coords[i][2]);
        fprintf(f, "TER\n");
    }

    fclose(f);
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


//void copy_coordinates(double **from, double **to)
//{
//    for (int i=0; i<n_atoms; ++i)
//    {
//        to[i][0] = from[i][0];
//        to[i][1] = from[i][1];
//        to[i][2] = from[i][2];
//    }
//}

double **initialisation_coords(const double *box_size)
{
    auto **coords = new double*[n_atoms];
    //double **old_coords = new double*[n_atoms];


    //Randomly generate the coordinates of the atoms in the box
    for (int i = 0; i < n_atoms; i = i + 1)
    {
        coords[i] = new double[3];
        //old_coords[i] = new double[3];

        // Note "rand(0,x)" would generate a random number
        // between 0 and $x
        coords[i][0] = rand(0, box_size[0]);
        coords[i][1] = rand(0, box_size[1]);
        coords[i][2] = rand(0, box_size[2]);
    }

    return coords;
}
//Initialise starting velocities for simulation - multiplication by
//the scalefactor is equivalent to Simulation.SetVelocitiestoTemperature() in OpenMM
double **initialisation_velocities(const double T, const double kb, const double mass)
{
    auto **vels = new double *[n_atoms];
    auto *v = new double [n_atoms]; //Velocity Vector - needed for calculating KE
    const double kT = kb * T;

    double sum_vx = 0.0;
    double sum_vy = 0.0;
    double sum_vz = 0.0;
    double sum_v2 = 0.0;

    //default_random_engine gen_random;
    //normal_distribution<double> P_v(0,0.15);

    for (int i = 0; i < n_atoms; i++)
    {
        vels[i] = new double[3];

        vels[i][0] = rand(-0.5, 0.5); //P_v(gen_random);
        vels[i][1] = rand(-0.5, 0.5); //P_v(gen_random);
        vels[i][2] = rand(-0.5, 0.5); //P_v(gen_random);

        v[i] = (vels[i][0]*vels[i][0]) +
                (vels[i][1]*vels[i][1]) + (vels[i][2]*vels[i][2]);

        sum_vx += vels[i][0];
        sum_vy += vels[i][1];
        sum_vz += vels[i][2];
        sum_v2 += v[i];
    }

    sum_vx = sum_vx/n_atoms;
    sum_vy = sum_vy/n_atoms;
    sum_vz = sum_vz/n_atoms; //Setting component-dependent center of mass

    const double scalefactor = sqrt( (3.0*n_atoms*kT) / (n_atoms*mass)*sum_v2); //scalefactor is T/(T(t))^1/2

    for (int j = 0; j < n_atoms; j++)
    {
        vels[j][0] = (vels[j][0] - sum_vx)*scalefactor;
        vels[j][1] = (vels[j][1] - sum_vy)*scalefactor;
        vels[j][2] = (vels[j][2] - sum_vz)*scalefactor;
    }

    return vels;

}

// Subroutine that calculates the potential energies of the atoms
double calculate_LJpot(double **coords, const double *box_size,
                       const double sigma, const double epsilon)
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
            const double sig2_over_r2 = (sigma*sigma) / r2;
            const double sig6_over_r6 = sig2_over_r2*sig2_over_r2*sig2_over_r2;
            const double sig12_over_r12 = sig6_over_r6 * sig6_over_r6;

            const double e_lj = 4.0 * epsilon * ( sig12_over_r12 - sig6_over_r6 );

            pot_energy = pot_energy + e_lj;
        }
    }

    // return the total energy of the atoms
    return pot_energy;
}

// Subroutine that calculates the force components of the atoms
double **calculate_LJforces(double **coords, const double *box_size,
                       const double sigma, const double epsilon)
{
    // Loop over all pairs of atoms and calculate
    // the LJ energy
    auto **f = new double *[n_atoms];

    const double r2_c = 64; //angstroms^2

    for (int t = 0; t < n_atoms; t++){

        f[t] = new double[3];

        f[t][0] = 0.0;
        f[t][1] = 0.0;
        f[t][2] = 0.0;

    }

    for (int i = 0; i < n_atoms-1; i = i + 1)
    {
        for (int j = i+1; j < n_atoms; j = j + 1)
        {

            double delta_x = coords[i][0] - coords[j][0];
            double delta_y = coords[i][1] - coords[j][1];
            double delta_z = coords[i][2] - coords[j][2];

            // Apply periodic boundaries
            delta_x = make_periodic(delta_x, box_size[0]);
            delta_y = make_periodic(delta_y, box_size[1]);
            delta_z = make_periodic(delta_z, box_size[2]);

            const double r2 = (delta_x*delta_x) + (delta_y*delta_y) +
                              (delta_z*delta_z);

            if (r2 < r2_c)
            {
                // F_LJ = 48*epsilon*1/r2[ (sigma/r)^12 - 0.5*(sigma/r)^6 ]
                const double one_over_r2 = 1.0 / r2;
                const double sig2_over_r2 = (sigma*sigma) / r2;
                const double sig6_over_r6 = sig2_over_r2*sig2_over_r2*sig2_over_r2;
                //const double sig12_over_r12 = sig6_over_r6 * sig6_over_r6;

                const double f_x = 48.0 * epsilon * one_over_r2 * sig6_over_r6 * ( sig6_over_r6 - 0.5 );

                //cout << "force = " << f_x << "\t" << "radial comp = " << r2
                //<< "\t" << delta_x << "\t" << delta_y << "\t" << delta_z << endl;

                f[i][0] = f[i][0] + (f_x * delta_x);
                f[i][1] = f[i][1] + (f_x * delta_y);
                f[i][2] = f[i][2] + (f_x * delta_z); //update force components for atom i

                f[j][0] = f[j][0] - (f_x * delta_x);
                f[j][1] = f[j][1] - (f_x * delta_y);
                f[j][2] = f[j][2] - (f_x * delta_z); //update force components for atom j

            }

        }
    }

    // return the total energy of the atoms
    return f;
}

double calc_kinetic(double **vels,  const double mass)
{
    double K = 0.0;
    auto *v = new double [n_atoms];

    for (int i = 0; i < n_atoms; i++)
    {
        v[i] = (vels[i][0]*vels[i][0]) +
               (vels[i][1]*vels[i][1]) + (vels[i][2]*vels[i][2]);

        K += 0.5*mass*v[i]/(3*n_atoms);
    }

    return K;
}

struct MDarray{
    double * e_pot;
    double * e_kin;
    double * e_tot;
};


MDarray VelocityVerletIntegrator(double **x, double **v, double **F,
        const double sigma, const double epsilon, const double *box_size, const double mass, const double dt, const int inc_print = 10)
{

    const int length = MDSteps;

    static double  Epot[length];
    static double  Ekin[length];
    static double  Etot[length];

    double vec[n_atoms];

    for (int step = 0; step <= length; step++)
    {

        double sum_v2 = 0.0, sum_v = 0.0;

        double U = calculate_LJpot(x, box_size, sigma, epsilon);
        Epot[step] = U;

        double K = calc_kinetic(v, mass);
        Ekin[step] = K;

        double Tot = U + K;
        Etot[step] = Tot;

        for (int atom = 0; atom < n_atoms; atom = atom + 1)
        {
            x[atom][0] = x[atom][0] + (v[atom][0] * dt) + ((0.5 * dt * dt * F[atom][0]) / mass);
            x[atom][1] = x[atom][1] + (v[atom][1] * dt) + ((0.5 * dt * dt * F[atom][1]) / mass);
            x[atom][2] = x[atom][2] + (v[atom][2] * dt) + ((0.5 * dt * dt * F[atom][2]) / mass);


            v[atom][0] = v[atom][0] + ((0.25 * dt * F[atom][0]) / mass);
            v[atom][1] = v[atom][1] + ((0.25 * dt * F[atom][1]) / mass);
            v[atom][2] = v[atom][2] + ((0.25 * dt * F[atom][2]) / mass);
        }

        F = calculate_LJforces(x, box_size, sigma, epsilon);

        for (int atom = 0; atom < n_atoms; atom = atom + 1 ) {

            v[atom][0] = v[atom][0] + ((0.25 * dt * F[atom][0]) / mass);
            v[atom][1] = v[atom][1] + ((0.25 * dt * F[atom][1]) / mass);
            v[atom][2] = v[atom][2] + ((0.25 * dt * F[atom][2]) / mass);

            vec[atom] = (v[atom][0]*v[atom][0])
                   + (v[atom][1]*v[atom][1]) + (v[atom][2]*v[atom][2]);

            sum_v2 += vec[atom];
            sum_v += v[atom][0] + v[atom][1] + v[atom][2];

        }

        //Print statements
        if (step % inc_print == 0) {

            const double temp_new = mass * sum_v2 / (3 * n_atoms * k_boltz);

            cout << "                                                          " << endl;
            cout << "Temperature - " << temp_new << " -- Center of mass for velocity - "
                 << sum_v << endl;

            cout << "Positions and velocities - Step " << step << " of " << length << endl;
            cout << "                                                                   " << endl;
            for (int i = 0; i < n_atoms; i++) {

                cout << x[i][0] << "\t" << x[i][1] << "\t" << x[i][2] << endl;
            }

            cout << "--------------------------------------------------------------" << endl;

            for (int i = 0; i < n_atoms; i++) {

                cout << v[i][0] << "\t" << v[i][1] << "\t" << v[i][2] << endl;
            }
            cout << "                                                                  " << endl;
            cout << "Potential = " << Epot[step] << " " << ",Kinetic = " << Ekin[step] << endl;
            cout << "Total Energy = " << Etot[step] << endl;

            print_pdb(x, n_atoms, step + 1);

        }

    }

    MDarray Traj;
    Traj.e_pot = Epot;
    Traj.e_kin = Ekin;
    Traj.e_tot = Etot;

    return Traj;
}

int main(int argc, const char **argv)
{
    //Set initial velocities and positions

    double **v_0 = initialisation_velocities(temperature, k_boltz, mass);
    double **x_0 = initialisation_coords(box_size);
    double **f_0 = calculate_LJforces(x_0, box_size, sigma, epsilon);

    print_pdb(x_0, n_atoms, 0);

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

    const double DeltaT = pow(2.0, -15);

    T = VelocityVerletIntegrator(x_0, v_0, f_0, sigma, epsilon, box_size, mass, DeltaT, 1000);

    return 0;
}
