//
// Created by Khaled Maksoud on 2019-04-09.
//

#include <iostream>
#include <cstdlib>
#include <cstdio>

int main()
{
    char filename[64];

    snprintf(filename, 64, "output_%d.pdb", 0);

    puts(filename);

    FILE *f = fopen(filename, "w");

    fprintf(f, "CRYST1 %8.3f %8.3f %8.3f  90.00  90.00  90.00\n", 0.1, 0.2, 0.3);

    fclose(f);


    return 0;
}

