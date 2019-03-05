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

vector<string> globVector(const string& pattern){
    glob_t glob_result;
    glob(pattern.c_str(),GLOB_TILDE,NULL,&glob_result);
    vector<string> files;
    for(unsigned int i=0;i<glob_result.gl_pathc;++i){
        files.push_back(string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
    return files;
}

void Configurations(vector<vector<double>>& data, double lat,int n)  {

	int crystal=4;
	int natoms = crystal*pow(n,3);
    data.resize(natoms,vector<double>(3));

    // Primitive cell for different crystal structures
    vector<vector <double> > pr;
	pr.resize(crystal,vector<double>(3));

	pr[0][0] = 0; pr[0][1] = 0; pr[0][2] = 0;
	pr[1][0] = 0.5; pr[1][1] = 0.5; pr[1][2] = 0;
	pr[2][0] = 0.5; pr[2][1] = 0; pr[2][2] = 0.5;
	pr[3][0] = 0; pr[3][1] = 0.5; pr[3][2] = 0.5;

    // Coordinates
    int t=0;
    for (int x=0; x<n; x++)  {
        for (int y=0; y<n; y++)  {
            for (int z=0; z<n; z++)  {
                for (int i=0; i<crystal; i++)  {

                    data[t][0] = (pr[i][0]+x)*lat;
                    data[t][1] = (pr[i][1]+y)*lat;
                    data[t][2] = (pr[i][2]+z)*lat;
                    t++;

                }
            }
        }   
    }   

	double max = 0.045*lat, min = -0.045*lat;

	for (int i=0; i<natoms; i++)  {
		for (int j=0; j<3; j++)   {

			data[i][j] += (max - min) * ( (double)rand() / (double)RAND_MAX ) + min; 
			if (data[i][j] < 0) { data[i][j] = abs(data[i][j]); }
			// cout << data[i][j] << endl;
		}         
	} 
}

void Write_POSCAR(double lat, vector<vector<double>> data, int numsim, int n)   {

    // Creates POSCAR file of a single cubic fcc cell of Al for VASP

    int natoms = data.size();
    int crystal=4;

    ofstream newfile; 
    newfile.open ("/home/tolga/Documents/mlip/Aluminum/ML/Bulk_Al/Data_108atoms/POSCAR/POSCAR"+to_string(numsim));

	if (crystal == 4) {newfile << "Al_FCC\n" << 1 << endl;}
    
    newfile << n*lat << " " << 0.0000 << " " << 0.0000 << endl;
    newfile << 0.0000 << " " << n*lat << " " << 0.0000 << endl;
    newfile << 0.0000 << " " << 0.0000 << " " << n*lat << endl;

    newfile << natoms << "\nc" << endl;

    for (int i=0; i<natoms; i++)   {

        if (i==0)   {newfile << data[i][0] << " " << data[i][1] << " " << data[i][2] << " positions" << endl;}
        else {newfile << data[i][0] << " " << data[i][1] << " " << data[i][2] << endl;}
    }

    newfile.close();   
        
}

int Read_Outcar(string filen)   {


    int natoms;
    string line, lined;
    double sxx, syy, szz, sxy, syz, sxz, toten, len;

    vector<vector <double> > mldata;
    vector<vector <double> > sphdata;
    
    ifstream readfile;  
    readfile.open(filen);

    while (true)  {

        readfile >> line;

        if (line == "type") {
			
			readfile >> line;
            readfile >> natoms;     
            break;    

        }
    }	
    return natoms;
}

void Create_DFT_Confs(vector<vector<double>>& data, double lat, double mass, int numsim)   {

	int crystal=4, n=3;
	
	for (int i=270; i<numsim+270; i++)  {
		
		Configurations(data,lat,n);
		Write_POSCAR(lat,data,i,n);			
			
	}
}


int main()   {

/*
    string fdump = "/home/tolga/Documents/github/lj_fitting/LJ_fitting/Diverse_Set/OUTCAR/OUTCAR*";
    vector<string> files = globVector(fdump);     
    
    int natoms;
    vector<string> outcars108;
		
	// Find the OUTCARs with 108 atoms
    for (int i=0; i<files.size(); i++)  {    

		natoms = Read_Outcar(files[i]);   
		
		if (natoms == 108) { outcars108.push_back(files[i]);}

    }
   
 
    // And take those files
    for (int i=0; i<outcars108.size()-1; i++)  {    

		//cout << "cp "+outcars108[i]+" /home/tolga/Documents/mlip/Aluminum/ML/Bulk_Al/Data_108atoms/" << endl;
		system(("cp "+outcars108[i]+" /home/tolga/Documents/mlip/Aluminum/ML/Bulk_Al/Data_108atoms/").c_str());	

    }  
    
    // Re-order
    for(int i=0; i<outcars108.size()-1; i++)  {
		
		system(("mv "+outcars108[i]+" /home/tolga/Documents/mlip/Aluminum/ML/Bulk_Al/Data_108atoms/OUTCAR/OUTCAR"+to_string(i)).c_str());
		
	}    
*/  


    double lat = 4.05, mass = 26.98; 
	int nsim = 30; // # of simulation for each crystal structure & size
    vector<vector <double> > data; 
    Create_DFT_Confs(data,lat,mass,nsim); 

  
}
