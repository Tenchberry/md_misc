//
// Created by Khaled Maksoud on 2019-04-11.
//

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>

#define MDSteps 60000

using namespace std;

// Set initial functions for calculating energies + forces

// Set the number of atoms in the box
const int n_atoms = 25;

// Set the size of the box (in Angstroms)
double box_size[3] = { 25.0, 25.0, 25.0 };

// Simulation temperature
const double temperature = 100;   // kelvin

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
const double sigma = 3.624;     // angstroms
const double epsilon = 0.317;   // kcal mol-1


double **position_3Dupdate(double **pos, double **vel, double **f, const float dt, const double stepfraction=1.0)
{

    auto **pos_new = new double *[n_atoms];

    for (int i = 0; i < n_atoms; i++)
    {
        pos_new[i] = new double[3];

        pos_new[i][0] = pos[i][0] + (vel[i][0]*dt*stepfraction) + (0.5 * dt * dt * f[i][0]);
        pos_new[i][1] = pos[i][1] + (vel[i][1]*dt*stepfraction) + (0.5 * dt * dt * f[i][0]);
        pos_new[i][2] = pos[i][2] + (vel[i][2]*dt*stepfraction) + (0.5 * dt * dt * f[i][0]);

    }


    return pos_new;
}

double **velocity_3Dupdate(double **vel, double **F, const int mass, const float dt,
                           const double kb, const double stepfraction=1.0)
{

    auto **vel_new = new double *[n_atoms];
    auto *v = new double [n_atoms];
    double sum_v2 = 0.0;

    for (int i=0; i < n_atoms; i++)
    {
        vel_new[i] = new double[3];

        vel_new[i][0] = vel[i][0] + (0.5*dt*stepfraction*F[i][0])/mass;
        vel_new[i][1] = vel[i][1] + (0.5*dt*stepfraction*F[i][1])/mass;
        vel_new[i][2] = vel[i][2] + (0.5*dt*stepfraction*F[i][2])/mass;

        v[i] = (vel_new[i][0]*vel_new[i][0])
               + (vel_new[i][1]*vel_new[i][1]) + (vel_new[i][2]*vel_new[i][2]);

        sum_v2 += v[i];

    }

    const double temp_new = mass*sum_v2 /(3*n_atoms*kb);

    cout << "Temperature - " << temp_new << endl;

    return vel_new;
}

//        cout << "                                                          " << endl;
//        cout << "Printing recalculated forces for step " << step << "/" << length << endl;
//
//        for (int i = 0; i < n_atoms; i++)
//        {
//            cout << F[i][0] << "\t" << F[i][1] << "\t" << F[i][2] << endl;
//        }

