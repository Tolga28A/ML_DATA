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
#include <ctime>
#include <random>
#include <time.h>
#include <iomanip>
#include <bits/stdc++.h> 
#include <boost/accumulators/statistics/variance.hpp>


using namespace std;

int RandomNum(int nu) {
			
    mt19937 rng;
    rng.seed(std::random_device()());
    uniform_int_distribution<std::mt19937::result_type> dist6(1,nu); // distribution in range [1, 6]

    return dist6(rng);

}

vector<double> ReadDump(string file_name) {
	
	string line;
	ifstream read_dump;  
    read_dump.open(file_name);
    
    vector<double> ycoord;
    int natoms;
    
    while (true)  {

       read_dump >> line;

       if ( line == "ATOMS" )  { 
                    
           read_dump >> natoms;
           ycoord.resize(natoms);
           break;

        }
    }    
    
    while ( !read_dump.eof() )  {

        read_dump >> line;

        if ( line == "c_eng" )  { 

			for (int i=0; i<natoms; i++)  {
				for (int j=0; j<7; j++)  {
					
					if (j == 3)  { read_dump >> ycoord[i]; }
					else { read_dump >> line; }

				}
			}       
        }        
    }   
    
    read_dump.close();	
    
    return ycoord;
}

int main()   {

	// File name
	string file_name = "/home/tolga/Documents/github/cpp_analysis_codes-master/build/dump.2697";
	// Get the y-coordinates
	vector<double> ycoord = ReadDump(file_name);
	
	// GB and Surface atoms boundary
	int gb_ymin=-20, gb_ymax=5;
	int surf_ymin=-95, surf_ymax=95;
		
	// Indices for 3 different region	
	vector<int> gb_ind;
	vector<int> bulk_ind;
	vector<int> surf_ind;
		
	for (int i=0; i<ycoord.size(); i++)  {
		
		if (ycoord[i] > gb_ymin && ycoord[i] < gb_ymax) {gb_ind.push_back(i);}
		else if (ycoord[i] < surf_ymin || ycoord[i] > surf_ymax) {surf_ind.push_back(i);}
		else {bulk_ind.push_back(i);}		
		
		//cout << i << " " << ycoord[i] << endl;
		
	}
	
	// # of environments for each 3 cases
	int ngb = 100, nbulk = 10, nsurf = 15;
	
	vector<int> sel_gb_ind;
	vector<int> sel_bulk_ind;
	vector<int> sel_surf_ind;		
	
	vector<int> selected_ids;



	
	// Randomly choose AE indices
	int r;
	
	for (int i=0; i<ngb; i++) {
		
		r = RandomNum(gb_ind.size());
		selected_ids.push_back(r);
		
	}
/*
	for (int i=0; i<nbulk; i++) {
		
		r = RandomNum(bulk_ind.size());
		selected_ids.push_back(r);
		
	}
	
	for (int i=0; i<nsurf; i++) {
		
		r = RandomNum(surf_ind.size());
		selected_ids.push_back(r);
		
	}	
*/	
	// Now find the relaxed POSCARs 
	for (int i=0; i<selected_ids.size(); i++)  {
	
		cout << selected_ids[i] << endl;
		system(("cp ~/Documents/github/cpp_analysis_codes-master/build/Poscars/LJ/TRY2/POSCAR"+to_string(selected_ids[i])+" ~/Documents/mlip/Aluminum/ML/cpp_codes/POSCAR"+to_string(i)).c_str());
	
	}
	
	// Compress them and delete the leftovers
	system("tar czvf new_selected_part2.tar.gz POSCAR*");
	system("rm POSCAR*");
}
