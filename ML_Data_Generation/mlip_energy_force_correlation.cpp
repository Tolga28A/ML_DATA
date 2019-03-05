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

void Write_Matlab_Data(vector<vector<double>> force_data,double energy,string source)  {
	
	ofstream matlab_force_file;
    matlab_force_file.open ("/home/tolga/Documents/mlip/Aluminum/ML/Derived_LJ_based_training_data/" + source + "_force_matlab.dat", ios_base::app);
    
    for (int i=0; i<force_data.size(); i++)  {
	
		for (int j=0; j<3; j++)  {
			
			matlab_force_file << force_data[i][j] << endl; 
			
		}
	}
	
	ofstream matlab_energy_file;
    matlab_energy_file.open ("/home/tolga/Documents/mlip/Aluminum/ML/Derived_LJ_based_training_data/" + source + "_energy_matlab.dat", ios_base::app);	
    
    matlab_energy_file << energy << endl;
	
}

void Read_Forces(string source)  {
	
	int natoms;
	double dum, energy;
    string line;
    vector<vector<double>> force_data;
    
    ifstream force;      
    force.open("/home/tolga/Documents/mlip/Aluminum/ML/Derived_LJ_based_training_data/" + source + ".cfg");	

    while (!force.eof())  {

        force >> line;

        if (line == "Size") {
			
            force >> natoms; 
            //cout << natoms << endl;
            force_data.resize(natoms,vector<double>(8));    
            
            while(true)  {
				
				force >> line;
				
				if (line == "fz")  {
					
					for (int i=0; i<natoms; i++)  {
						
						for (int j=0; j<5; j++)  {force >> dum;}
						for (int k=0; k<3; k++)  {force >> force_data[i][k];}					
						
					}
					
					force >> line;
					//cout << line << endl; 
					force >> energy;
					//cout << energy << endl; 
					break;	
									
				}
			}
			
			Write_Matlab_Data(force_data,energy,source);
			
        }                   
    }
	
}


int main()   {

	string file1 = "training_set/mlip_al_lj41";
	string file2 = "fitted_potential_forces/Al_trained_forces_F05S01_epoch1";

	//Read_Forces(file1);
	Read_Forces(file2);
	
}
