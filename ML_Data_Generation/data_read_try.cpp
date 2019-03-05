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


int Atom_Remove()  {

    vector<double> maxen;
    vector<double> numofa;
    vector<double> id;
    vector<double> en;

    int numv;
    int max;

    string line;
    ifstream readen;  
    readen.open("/home/tolga/Documents/github/cpp_analysis_codes-master/build/atomenergy.dat");

    while (true)  {

       readen >> line;

       if ( line == "ATOMS" )  { 
                    
           readen >> numv;
           numofa.push_back(numv);
           break;

        }
    }

    id.resize(numv);  en.resize(numv);

    while ( !readen.eof() )  {

       readen >> line;

       if ( line == "c_atomenergyc[*]" )  { 

           for (int i=0; i<numv; i++)  {

                readen >> id[i];
                readen >> en[i];
                
           }       
        }        
    }   
    
    readen.close();
    
    for (int i=0; i<numv; i++)  {
    
		cout << i << endl;
		cout << id[i] << endl;
		cout << en[i] << endl;
                
	}


    auto greatest = max_element(en.begin(), en.end());
    maxen.push_back(*greatest);
    max = distance(begin(en), greatest);

    return id[max];

}


int main()   {

    int atomrem = Atom_Remove();
	
}
