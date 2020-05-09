//
// Created by Khaled Maksoud on 2019-04-09.
//

#include <iostream>
#include <fstream>

using namespace std;

int main()
{
    ofstream file;
    file.open("output%d.pdb", 0);

    if (file.is_open())
    {
        file << "This is a line.\n";
        file << "This is another line.\n";
        file.close();
    }
    else cout << "Unable to open file";

    return 0;
}
