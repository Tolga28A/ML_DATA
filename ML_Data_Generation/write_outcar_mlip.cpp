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


void Write_MLIP(vector<vector<double>> data, double len, double energy, double sxx, double syy, double szz, double sxy, double sxz, double syz)  {

    int numofa = data.size();
    
    ofstream newfile;
    newfile.open ("mlip_al_lj59.cfg", ios_base::app);

    newfile << "\n BEGIN_CFG" << endl;
    newfile << " Size" << endl;
    newfile << "\t" << numofa << endl;
    newfile << " SuperCell" << endl;
    newfile << " \t\t" << len << "\t 0.0 \t 0.0" << endl;
    newfile << " \t\t" << "0.0 \t " << len << "\t 0.0" << endl;
    newfile << " \t\t" << "0.0 \t " << "0.0 \t " << len << "\n" << endl;
    newfile << "AtomData: \t id \t type \t cartes_x \t cartes_y \t cartes_z \t fx \t fy \t fz" << endl;

    for (int i=0; i<numofa; i++) {
        newfile << "\t\t";
        for (int j=0; j<8; j++)  {

            newfile << "\t" << data[i][j] << "\t";

        }
        newfile << "\n";
    }

    newfile << " Energy" << endl;
    newfile << "\t" << energy << endl;
    newfile << " Stress: xx \t yy \t zz \t yz \t xz \t xy" << endl;
    newfile << "\t" << sxx << "\t" << syy << "\t" << szz << "\t" << syz << "\t" << sxz << "\t" << sxy << "\n" << endl;
    newfile << " Feature   from  database:mlip.cfg" << endl;

    newfile << "END_CFG \n" << endl;
//    newfile.close();

}

void Read_Outcar(int nsim)   {


    int natoms;
    string line1, line2, lined;
    double sxx, syy, szz, sxy, syz, sxz, toten, len;

    vector<vector <double> > mldata;
    
    ifstream readfile;  
    readfile.open("/home/tolga/Documents/mlip/Aluminum/ML/Derived_LJ_based_training_data/lj_relaxed59_dft_data/OUTCAR"+to_string(nsim));

    while (true)  {

        readfile >> line1;
        readfile >> line2;

        if (line1 == "type" && line2 == "=") {
			
            readfile >> natoms; 
            //cout << natoms << endl;
            mldata.resize(natoms,vector<double>(8));       
            break;    

        }
    }

    while (!readfile.eof())  {
		
        readfile >> line1;
        readfile >> line2;
        //cout << line1 << endl;

        if (line1 == "TOTEN" && line2 == "=") {

            readfile >> toten; 
            //cout << toten << endl;

        }

        if (line2 == "Total")  {

            readfile >> sxx;
            readfile >> syy;
            readfile >> szz;
            readfile >> sxy;
            readfile >> syz;
            readfile >> sxz;

        }


        if (line1 == "of" && line2 == "vectors")  {

            readfile >> len;
            

        }


        if (line1 == "(eV/Angst)")  {

            for (int i=0; i<natoms; i++)  {
                for (int j=2; j<=7; j++)  {
                
                    readfile >> mldata[i][j];
                    //cout << mldata[i][j] << endl;

                }
            }  

            break;
        }
    }
    readfile.close();

    for (int i=0; i<natoms; i++)   {

        mldata[i][0] = i+1;
        mldata[i][1] = 0;   

    }
    
	if (mldata[10][3] == 0 && mldata[50][4] == 0 && mldata[110][5] == 0)  {
		
		cout << nsim << endl;
		
	}
	
	else { Write_MLIP(mldata,len,toten,sxx,syy,szz,sxy,syz,sxz); }

/*
    cout << "numofa: " << natoms << endl;
    cout << "energy: " << toten << endl;
    cout << "length: " << len << endl;
    cout << "sxx: " << sxx << endl;
    cout << "syy: " << syy << endl;    
    cout << "szz: " << szz << endl;
    cout << "sxy: " << sxy << endl;
    cout << "sxz: " << sxz << endl;
    cout << "syz: " << syz << endl;
*/


}


int main()   {

    string fdump = "/home/tolga/Documents/mlip/Aluminum/ML/Derived_LJ_based_training_data/lj_relaxed59_dft_data/OUTCAR*";
    vector<string> files = globVector(fdump); 
    int nsim= files.size(); 
    cout << nsim << endl;

    for (int i=0; i<41; i++)  {    

		//if (i > 22 & i < 38) {continue;}
		//else if (i == 10 || i == 11) {continue;}	
		//else {Read_Outcar(i);}        
		Read_Outcar(i);

    }

}

