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
#include <algorithm>

using namespace std;

int Find_Center_Atom(vector<vector<double>> force_data,double length)  {
	
	vector<double> dist;
	double r, mid=length/2;
	for (int i=0; i<force_data.size(); i++)  {
		
		r = sqrt(pow((force_data[i][2]-mid),2) + pow((force_data[i][3]-mid),2) + pow((force_data[i][4]-mid),2));
		dist.push_back(r);
		
	}			
	
	auto min = min_element(dist.begin(),dist.end());
	double min_val = *min;
	return distance(begin(dist),min);
	
	
}

void Dist_Force_Write(vector<vector<double>> force_data,int ind)  {
	
	vector<double> fx;
	vector<double> fy;
	vector<double> fz;
	double r;
	
	//cout << ind << endl;
	
	ofstream matlab_force_file;
    matlab_force_file.open ("/home/tolga/Documents/mlip/Aluminum/ML/EAM_based_training_data/Bulk_Al/training_set/training_forces.dat", ios_base::app);
    
	
	for (int i=0; i<force_data.size(); i++)  {
		
		if (i == ind) {continue;}
		
		else {
			
			r = sqrt(pow((force_data[i][2]-force_data[ind][2]),2) + pow((force_data[i][3]-force_data[ind][3]),2) + pow((force_data[i][4]-force_data[ind][4]),2));
			
			matlab_force_file << r << " ";
			matlab_force_file << force_data[i][5] << " ";
			matlab_force_file << force_data[i][6] << " ";
			matlab_force_file << force_data[i][7] << endl;
			
		}
	}
}

void Training_Forces(string file)  {

	
	int natoms;
	double dum, length;
    string line;
	vector<vector<double>> force_data;    
    
    ifstream force;      
    force.open(file);		
    
    while (!force.eof())  {

        force >> line;

        if (line == "Size") {
			
            force >> natoms; 
            //cout << natoms << endl;
            force_data.resize(natoms,vector<double>(8));  
            
		}
            
        else if (line == "SuperCell") {
        
			force >> length;
			//cout << length << endl;
			
		}

		else if (line == "fz")  {
			
			for (int i=0; i<natoms; i++)  {
							
				for (int j=0; j<8; j++)  {
			  
					force >> force_data[i][j];
					
				}					
				
			}
		
			// Find the index of the atom that is closest to the center
			// Which is very likely the center atom of AE 
			int mind_ind = Find_Center_Atom(force_data,length);	
			// And save the results
			Dist_Force_Write(force_data,mind_ind);		
									
		}						
	} 	                
}


int main()   {

	//string file = "/home/tolga/Documents/mlip/doc/examples/2.fitting_error_calculation/Li_NVT300.cfg";
	string file = "/home/tolga/Documents/mlip/Aluminum/ML/EAM_based_training_data/Bulk_Al/training_set/training_set.cfg";	
	Training_Forces(file);
	
}
