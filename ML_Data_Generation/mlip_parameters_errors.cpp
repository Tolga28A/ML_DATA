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

string Set_String_Precision(double value,int n)  {
	
	ostringstream out;
	out.precision(n);
	out << fixed << value;
	return out.str();
	
}

void Read_Errors(string er_file, double& epa_ave, double& epa_rms, double& fmax, double& frms, double& smax, double& srms)   {

    double dum;
    string line;
    
    ifstream error;      
    error.open(er_file);	

    while (!error.eof())  {
		
		error >> line;

        if (line == "per") {
			
			//error >> line;
			//cout << line << endl;
					
			for (int i=0; i<21; i++) {
							
				if (i == 15)  {error >> epa_ave;}
				else if (i == 20)  {error >> epa_rms; break;}
				else {error >> line;}
				
			}
        }
        
        else if (line == "Forces:") {
			
			for (int i=0; i<21; i++) {
							
				if (i == 14)  {error >> fmax;}
				else if (i == 19)  {error >> frms; break;}
				else {error >> line;}
				
			}			
        }
        
        else if (line == "GPa):") {
			
			for (int i=0; i<21; i++) {
							
				if (i == 14)  {error >> smax;}
				else if (i == 19)  {error >> srms; break;}
				else {error >> line;}
				
			}			
        }                                             
    }
    
    error.close();
	
}

int main()   {
	
	double epa_ave,epa_rms,fmax,frms,smax,srms;
	
	vector<string> f{"01","05","1"};
	
	ofstream energy_ave;
	energy_ave.open("energy_ave.dat");
	
	ofstream energy_rms;
	energy_rms.open("energy_rms.dat");
	
	ofstream fmaxf;
	fmaxf.open("fmax.dat");
	
	ofstream frmsf;
	frmsf.open("frms.dat");
		
	ofstream smaxf;
	smaxf.open("smax.dat");
	
	ofstream srmsf;
	srmsf.open("srms.dat");	
	
	

	for (int i=0; i<3; i++)  {
		for (int j=0; j<3; j++)  {
							
			string er_file = "/home/tolga/Documents/mlip/Aluminum/ML/Derived_LJ_based_training_data/error_files/lj_relaxed41_mlip_errors/errorsF"+f[i]+"S"+f[j]+".dat";	
			
			Read_Errors(er_file,epa_ave,epa_rms,fmax,frms,smax,srms);
			
			energy_ave << epa_ave << " ";
			energy_rms << epa_rms << " ";
			
			fmaxf << fmax << " ";
			frmsf << frms << " ";
			
			smaxf << smax << " ";
			srmsf << srms << " ";								
				
		}		
		
		energy_ave << endl;
		energy_rms << endl;
		fmaxf << endl;
		frmsf << endl;
		smaxf << endl;
		srmsf << endl;	
							
	}
	
	energy_ave.close();
	energy_rms.close();	
	fmaxf.close();
	frmsf.close();
	smaxf.close();
	srmsf.close();			

}
