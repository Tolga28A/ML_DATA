#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <cmath>
#include <cstdlib>
#include <dirent.h>
#include <glob.h>
#include <bits/stdc++.h>

using namespace std;

void Write_File(vector<vector<double>> data, int nsim, double simlen)   {

    int natoms = data.size();

    ofstream newfile; 
    newfile.open ("/home/tolga/Documents/mlip/Aluminum/ML/Derived_LJ_based_training_data/training_set/AL_selection/AL_selected_poscars/POSCAR"+to_string(nsim));

    newfile << "Al_FCC\n" << 1 << endl;
    newfile << simlen << " " << 0.0000 << " " << 0.0000 << endl;
    newfile << 0.0000 << " " << simlen << " " << 0.0000 << endl;
    newfile << 0.0000 << " " << 0.0000 << " " << simlen << endl;

    newfile << natoms << "\nc" << endl;

    for (int i=0; i<natoms; i++)   {

        if (i==0)   {newfile << data[i][0] << " " << data[i][1] << " " << data[i][2] << " positions" << endl;}
        else {newfile << data[i][0] << " " << data[i][1] << " " << data[i][2] << endl;}
    }

    newfile.close();   
        
}

int Write_Poscar()   {

    // This function is just to write POSCAR from the .cfg database
    int natoms;
    double simlen, dum;
    string line;
    ifstream infile;
    infile.open("/home/tolga/Documents/mlip/Aluminum/ML/Derived_LJ_based_training_data/training_set/AL_selection/milp41_selected_qs1.cfg");

    vector<vector <double> > data;

    int count = 0;
    while (!infile.eof())  {

       infile >> line;
       //cout << line << endl;

       if ( line == "Size" )  { 
		   
                    
           infile >> natoms;
           //cout << natoms << endl;
           data.resize(natoms,vector<double>(3)); 
           infile >> line;

       }        

       if ( line == "SuperCell" )  { 
                    
           infile >> simlen;

       }

       if ( line == "fz" )  { 
		   		   
                    
           for (int i=0; i<natoms; i++)  {

                infile >> dum;                        
                infile >> dum;  

                infile >> data[i][0];
                infile >> data[i][1];
                infile >> data[i][2];

                infile >> dum;                        
                infile >> dum; 
                infile >> dum;                         

                //cout << data[i][0] << " " << data[i][1] << " " << data[i][2] << endl;            
               
            }
            Write_File(data,count,simlen);   
            count++;
       }
    }
    return count;
}


int main()   {

	int nenv = Write_Poscar();
	cout << nenv << endl;

}
