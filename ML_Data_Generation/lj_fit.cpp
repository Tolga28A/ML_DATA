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

//*****************************************************************************************************************************************************//
//*****************************************************************************************************************************************************//
//*****************************************************************************************************************************************************//

void nelmin ( double fn ( double x[] ), int n, double start[], double xmin[], 
  double *ynewlo, double reqmin, double step[], int konvge, int kcount, 
  int *icount, int *numres, int *ifault )   {

  double ccoeff = 0.5;
  double del;
  double dn;
  double dnn;
  double ecoeff = 2.0;
  double eps = 0.001;
  int i;
  int ihi;
  int ilo;
  int j;
  int jcount;
  int l;
  int nn;
  double *p;
  double *p2star;
  double *pbar;
  double *pstar;
  double rcoeff = 1.0;
  double rq;
  double x;
  double *y;
  double y2star;
  double ylo;
  double ystar;
  double z;
//
//  Check the input parameters.
//
  if ( reqmin <= 0.0 )
  {
    *ifault = 1;
    return;
  }

  if ( n < 1 )
  {
    *ifault = 1;
    return;
  }

  if ( konvge < 1 )
  {
    *ifault = 1;
    return;
  }

  p = new double[n*(n+1)];
  pstar = new double[n];
  p2star = new double[n];
  pbar = new double[n];
  y = new double[n+1];

  *icount = 0;
  *numres = 0;

  jcount = konvge; 
  dn = ( double ) ( n );
  nn = n + 1;
  dnn = ( double ) ( nn );
  del = 1.0;
  rq = reqmin * dn;
//
//  Initial or restarted loop.
//
  for ( ; ; )
  {
    for ( i = 0; i < n; i++ )
    { 
      p[i+n*n] = start[i];
    }
    y[n] = fn ( start );
    *icount = *icount + 1;

    for ( j = 0; j < n; j++ )
    {
      x = start[j];
      start[j] = start[j] + step[j] * del;
      for ( i = 0; i < n; i++ )
      {
        p[i+j*n] = start[i];
      }
      y[j] = fn ( start );
      *icount = *icount + 1;
      start[j] = x;
    }
//                    
//  The simplex construction is complete.
//                    
//  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
//  the vertex of the simplex to be replaced.
//                
    ylo = y[0];
    ilo = 0;

    for ( i = 1; i < nn; i++ )
    {
      if ( y[i] < ylo )
      {
        ylo = y[i];
        ilo = i;
      }
    }
//
//  Inner loop.
//
    for ( ; ; )
    {
      if ( kcount <= *icount )
      {
        break;
      }
      *ynewlo = y[0];
      ihi = 0;

      for ( i = 1; i < nn; i++ )
      {
        if ( *ynewlo < y[i] )
        {
          *ynewlo = y[i];
          ihi = i;
        }
      }
//
//  Calculate PBAR, the centroid of the simplex vertices
//  excepting the vertex with Y value YNEWLO.
//
      for ( i = 0; i < n; i++ )
      {
        z = 0.0;
        for ( j = 0; j < nn; j++ )
        { 
          z = z + p[i+j*n];
        }
        z = z - p[i+ihi*n];  
        pbar[i] = z / dn;
      }
//
//  Reflection through the centroid.
//
      for ( i = 0; i < n; i++ )
      {
        pstar[i] = pbar[i] + rcoeff * ( pbar[i] - p[i+ihi*n] );
      }
      ystar = fn ( pstar );
      *icount = *icount + 1;
//
//  Successful reflection, so extension.
//
      if ( ystar < ylo )
      {
        for ( i = 0; i < n; i++ )
        {
          p2star[i] = pbar[i] + ecoeff * ( pstar[i] - pbar[i] );
        }
        y2star = fn ( p2star );
        *icount = *icount + 1;
//
//  Check extension.
//
        if ( ystar < y2star )
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = pstar[i];
          }
          y[ihi] = ystar;
        }
//
//  Retain extension or contraction.
//
        else
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = p2star[i];
          }
          y[ihi] = y2star;
        }
      }
//
//  No extension.
//
      else
      {
        l = 0;
        for ( i = 0; i < nn; i++ )
        {
          if ( ystar < y[i] )
          {
            l = l + 1;
          }
        }

        if ( 1 < l )
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = pstar[i];
          }
          y[ihi] = ystar;
        }
//
//  Contraction on the Y(IHI) side of the centroid.
//
        else if ( l == 0 )
        {
          for ( i = 0; i < n; i++ )
          {
            p2star[i] = pbar[i] + ccoeff * ( p[i+ihi*n] - pbar[i] );
          }
          y2star = fn ( p2star );
          *icount = *icount + 1;
//
//  Contract the whole simplex.
//
          if ( y[ihi] < y2star )
          {
            for ( j = 0; j < nn; j++ )
            {
              for ( i = 0; i < n; i++ )
              {
                p[i+j*n] = ( p[i+j*n] + p[i+ilo*n] ) * 0.5;
                xmin[i] = p[i+j*n];
              }
              y[j] = fn ( xmin );
              *icount = *icount + 1;
            }
            ylo = y[0];
            ilo = 0;

            for ( i = 1; i < nn; i++ )
            {
              if ( y[i] < ylo )
              {
                ylo = y[i];
                ilo = i;
              }
            }
            continue;
          }
//
//  Retain contraction.
//
          else
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
        }
//
//  Contraction on the reflection side of the centroid.
//
        else if ( l == 1 )
        {
          for ( i = 0; i < n; i++ )
          {
            p2star[i] = pbar[i] + ccoeff * ( pstar[i] - pbar[i] );
          }
          y2star = fn ( p2star );
          *icount = *icount + 1;
//
//  Retain reflection?
//
          if ( y2star <= ystar )
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
          else
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = pstar[i];
            }
            y[ihi] = ystar;
          }
        }
      }
//
//  Check if YLO improved.
//
      if ( y[ihi] < ylo )
      {
        ylo = y[ihi];
        ilo = ihi;
      }
      jcount = jcount - 1;

      if ( 0 < jcount )
      {
        continue;
      }
//
//  Check to see if minimum reached.
//
      if ( *icount <= kcount )
      {
        jcount = konvge;

        z = 0.0;
        for ( i = 0; i < nn; i++ )
        {
          z = z + y[i];
        }
        x = z / dnn;

        z = 0.0;
        for ( i = 0; i < nn; i++ )
        {
          z = z + pow ( y[i] - x, 2 );
        }

        if ( z <= rq )
        {
          break;
        }
      }
    }
//
//  Factorial tests to check that YNEWLO is a local minimum.
//
    for ( i = 0; i < n; i++ )
    {
      xmin[i] = p[i+ilo*n];
    }
    *ynewlo = y[ilo];

    if ( kcount < *icount )
    {
      *ifault = 2;
      break;
    }

    *ifault = 0;

    for ( i = 0; i < n; i++ )
    {
      del = step[i] * eps;
      xmin[i] = xmin[i] + del;
      z = fn ( xmin );
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] - del - del;
      z = fn ( xmin );
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] + del;
    }

    if ( *ifault == 0 )
    {
      break;
    }
//
//  Restart the procedure.
//
    for ( i = 0; i < n; i++ )
    {
      start[i] = xmin[i];
    }
    del = eps;
    *numres = *numres + 1;
  }
  delete [] p;
  delete [] pstar;
  delete [] p2star;
  delete [] pbar;
  delete [] y;

  return;

}

//*****************************************************************************************************************************************************//
//*****************************************************************************************************************************************************//
//*****************************************************************************************************************************************************//
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

void Calculate_LJ(double eps, double sigma, int nsim)   {

    //Return LJ energies
    double rcut = 8;
    
    ofstream infile;
    infile.open ("lj_calc.in");
    
    infile << "units	     metal" << endl;
    infile << "atom_style    atomic" << endl;
    infile << "dimension     3" << endl;
    infile << "boundary      p p p" << endl;
    infile << "timestep      0.001" << endl;
    infile << "log	         /home/tolga/Documents/github/lj_fitting/ljfitting_lammps_files/lj_fit"+to_string(nsim)+".lammps" << endl;
    infile << "read_data     /home/tolga/Documents/github/lj_fitting/ljfitting_lammps_files/lj_fit"+to_string(nsim)+".dat" << endl;
    infile << "pair_style lj/cut " << rcut << endl;
    infile << "pair_coeff    * * " << eps << " " << sigma << endl;
    infile << "neighbor      0.5 bin" << endl;
    infile << "neigh_modify  every 1 delay 0 check yes" << endl;
//    infile << "velocity      all create 1.0 4928459 dist gaussian" << endl;
    infile << "fix           NVE all nve" << endl;
    infile << "compute       lje all pe/atom" << endl;
    infile << "dump          ljenergy all custom 1 lje.dat c_lje[*]" << endl;
    infile << "dump_modify   ljenergy sort id" << endl;
    infile << "run           0" << endl;

    infile.close();

    system("nohup /home/tolga/Documents/LAMMPS/lmp_mpi < lj_calc.in"); 
}

void Write_Data_LJ(vector<vector<double>> data, double lat, int nsim, int n, double mass)   {

    //Run LAMMPS to calculate LJ energies    
    int natoms = data.size();

    ofstream newfile; 
    newfile.open ("lj_fit"+to_string(nsim)+".dat");

    newfile << "#Created by Tolga" << endl;
    newfile << "\t" << natoms << " atoms" << endl;
    newfile << "\t" << 0 << " bonds" << endl;
    newfile << "\t" << 0 << " angles" << endl;
    newfile << "\t" << 0 << " dihedrals" << endl;
    newfile << "\t" << 0 << " impropers \n" << endl;

    newfile << "\t" << 1 << " atom types" << endl;
    newfile << "\t" << 0 << " bond types" << endl;
    newfile << "\t" << 0 << " angle types" << endl;
    newfile << "\t" << 0 << " dihedral types" << endl;
    newfile << "\t" << 0 << " improper types \n" << endl;

    newfile << "\t" << 0 << " " << n*lat << " xlo xhi" << endl;  
    newfile << "\t" << 0 << " " << n*lat << " ylo yhi" << endl; 
    newfile << "\t" << 0 << " " << n*lat << " zlo zhi \n" << endl;  

    newfile << "Masses \n" << endl;
    newfile << 1 << " " << mass << "\n" << endl;

    newfile << "Atoms \n\n";

    for (int i=0; i<natoms; i++)  {
        newfile << i+1 << " " << 1 << " ";
        for (int j=0; j<3; j++)  {

            newfile << data[i][j] << " ";

        }
        newfile << "\n";
    }

    newfile.close();   
    
}

void Read_LJ(vector<double>& lje_data)   {

    int natoms = lje_data.size();
    
    string line;
    ifstream readfile;  
    readfile.open("lje.dat");    

    while (!readfile.eof())  {
    
        readfile >> line;

        if (line == "c_lje[*]")  {
            for (int i=0; i<natoms; i++)   {

                readfile >> lje_data[i];
                //cout << "lj: " << lje_data[i] << endl;
                                              
            }
            break;
        }
    }

}


void Write_POSCAR(double lat, vector<vector<double>> data, int numsim, int n, int crystal)   {

    // Creates POSCAR file of a single cubic fcc cell of Al for VASP

    int natoms = data.size();

    ofstream newfile; 
    newfile.open ("POSCAR"+to_string(numsim));

	if (crystal == 1) {newfile << "Al_SC\n" << 1 << endl;}
	if (crystal == 2) {newfile << "Al_BCC\n" << 1 << endl;}
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

void Configurations(vector<vector<double>>& data, double lat, int crystal, int n, int flag)  {

	int natoms = crystal*pow(n,3);
    data.resize(natoms,vector<double>(3));

    // Primitive cell for different crystal structures
    vector<vector <double> > pr;
	pr.resize(crystal,vector<double>(3));
	
	if (crystal == 4)  {

		pr[0][0] = 0; pr[0][1] = 0; pr[0][2] = 0;
		pr[1][0] = 0.5; pr[1][1] = 0.5; pr[1][2] = 0;
		pr[2][0] = 0.5; pr[2][1] = 0; pr[2][2] = 0.5;
		pr[3][0] = 0; pr[3][1] = 0.5; pr[3][2] = 0.5;
	}

	else if (crystal == 2)  {
	
		pr[0][0] = 0; pr[0][1] = 0; pr[0][2] = 0;
		pr[1][0] = 0.5; pr[1][1] = 0.5; pr[1][2] = 0.5;		

	}

	else if (crystal == 1)  {

		pr[0][0] = 0; pr[0][1] = 0; pr[0][2] = 0;

	}

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
    
    if (flag != 0)  {

        double max = 0.025*lat, min = -0.025*lat;

	    for (int i=0; i<natoms; i++)  {
		    for (int j=0; j<3; j++)   {

			    data[i][j] += (max - min) * ( (double)rand() / (double)RAND_MAX ) + min; 
                if (data[i][j] < 0) { data[i][j] = abs(data[i][j]); }
			    // cout << data[i][j] << endl;
		    }         
	    }
    }    
}

void Create_DFT_Confs(vector<vector<double>>& data, int flag, double lat, double mass)   {

	// Now I will have 9 times nsim simulations correspondin to 
	// 3 different crystal structures x 3 different sizes
	// Also the total # of simulations
	int count = 3*3*flag; 
	// 3 different crystal structures are of interest for now
	// I'm labeling them depending on the # of atoms in a primitive cell
	// 1 = SC, 2 = BCC, 4 = FCC 
	for (int n=2; n<=4; n++)  {
		// Now iterate over # of cells in one direction which will
		// give different volume of simulations
		for (int crystal=1; crystal<5; (crystal=2*crystal))  {
			
			Configurations(data,lat,crystal,n,flag);
			Write_POSCAR(lat,data,count,n,crystal);
			Write_Data_LJ(data,lat,count,n,mass);
			count++;			
				
		}
	}

}

void Read_DFT_Vasp(vector<double>& dft_data, int nsim, int& natoms)  {

    // Read number of atoms and energy from VASP-OUTCAR
    double toten;
    string line1, line2;
    ifstream readfile;  
    readfile.open("/home/tolga/Documents/github/lj_fitting/LJ_fitting/Diverse_Set/OUTCAR/OUTCAR"+to_string(nsim));

    while (true)  {

        readfile >> line1;
        readfile >> line2;

        if (line1 == "type" && line2 == "=") {
			
            readfile >> natoms; 
            dft_data.resize(natoms);       
            break;    
        }
    }


    while (!readfile.eof())  {

        readfile >> line1;
        readfile >> line2;
        if (line1 == "TOTEN" && line2 == "=") {readfile >> toten;}

    }


    for (int i=0; i<natoms; i++)  {

	    dft_data[i] = toten/natoms;
        //cout << "dft: " << dft_data[i] << endl;

    }
}

/*
double Objective_Function(double ljp[2])   {

    double lat = 4.05;
    double n = 4*pow(2,3);
    double natoms = 4;  
    double nsim = 10;  

    double eps = ljp[0];
    double sigma = ljp[1];

    // Atomic configuration data
    vector<vector <double> > data;
    data.resize(natoms,vector<double>(3));

    // LJ energy data
    vector<double> lje_data;
    lje_data.resize(natoms);

    // DFT energy data
    vector<double> dft_data;
    dft_data.resize(natoms);

    // Objective function to minimize 
    vector<double> objf;
    objf.resize(natoms*nsim);

    for (int i=0; i<nsim; i++)   {

	    // Get DFT energy results
	    Create_DFT_Confs(data,i,lat);
	    Write_POSCAR(lat,data,i);
	    ////// Run VASP //////
	    Read_DFT_Vasp(dft_data,i);

	    // Get LJ energy results
	    Write_Data_LJ(data,lat,nsim);
	    Calculate_LJ(eps,sigma,nsim);
	    Read_LJ(lje_data);

	    int k=0;
	    for (int j=0; j<natoms; j++)   {

	        objf[k] = (lje_data[j] - dft_data[j])*(lje_data[j] - dft_data[j]);
	        k++;
	        //cout << objf[i] << "  " << i << endl;       

	    }
    }

    double er2 = accumulate(objf.begin(), objf.end(), 0.0)/objf.size(); 

    return er2;
    
}
*/


double Objective_Function2(double ljp[2])   {

	// First find the # of files in the directory 
	// which is the total # of simulations
    string fdump = "/home/tolga/Documents/github/lj_fitting/LJ_fitting/Diverse_Set/OUTCAR/OUTCAR*";
    vector<string> files = globVector(fdump); 
    int nsim= files.size(); 
    cout << nsim << endl;

    double lat = 4.05;

    double eps = ljp[0];
    double sigma = ljp[1];

    // LJ energy data
    vector<double> lje_data;

    // DFT energy data
    vector<double> dft_data;

    // Objective function to minimize of its error-square
    vector<double> objf;

    ofstream enfile; 
    enfile.open ("energies.dat");

	int natoms;
	int k=0;
    for (int i=0; i<nsim; i++)   {

        Read_DFT_Vasp(dft_data,i,natoms);
	    Calculate_LJ(eps,sigma,i);
	    	    
		lje_data.resize(natoms);
		objf.resize(natoms*nsim);

	    Read_LJ(lje_data);	    

	    for (int j=0; j<natoms; j++)   {

	        objf[k] = (lje_data[j] - dft_data[j])*(lje_data[j] - dft_data[j]);
            enfile << lje_data[j] << " " << dft_data[j] << " " << objf[k] << endl;
	        k++;      
	    }  	     
    } 

    enfile.close();   

    double er2 = accumulate(objf.begin(), objf.end(), 0.0)/objf.size(); 
    return er2;
    
}

void find_min(double ljp[2])   {

    int i;
    int icount;
    int ifault;
    int kcount;
    int konvge;
    int n;
    int numres;
    double reqmin;
    double *start;
    double *step;
    double *xmin;
    double ynewlo;
  
    n = 2;

    start = new double[n];
    step = new double[n];
    xmin = new double[n];

    start[0] = ljp[0];
    start[1] = ljp[1];

    reqmin = 1.0E-08;

    step[0] = 1.0;
    step[1] = 1.0;

    konvge = 10;
    kcount = 5000;

    cout << "\n";
    cout << "  Starting point X:\n";
    cout << "\n";
    for ( i = 0; i < n; i++ )
    {
    cout << "  " << setw(14) << start[i] << "\n";
    }

    ynewlo = Objective_Function2(start);

    cout << "\n";
    cout << "  F(X) = " << ynewlo << "\n";

    nelmin (Objective_Function2, n, start, xmin, &ynewlo, reqmin, step,
    konvge, kcount, &icount, &numres, &ifault );

    cout << "\n";
    cout << "  Return code IFAULT = " << ifault << "\n";
    cout << "\n";
    cout << "  Estimate of minimizing value X*:\n";
    cout << "\n";
    for ( i = 0; i < n; i++ )
    {
    cout << "  " << setw(14) << xmin[i] << "\n";
    }

    cout << "\n";
    cout << "  F(X*) = " << ynewlo << "\n";

    cout << "\n";
    cout << "  Number of iterations = " << icount << "\n";
    cout << "  Number of restarts =   " << numres << "\n";

    delete [] start;
    delete [] step;
    delete [] xmin;

    return;

}


int main()  {


    double ljp[2];
//    ljp[0] = 3.92; //epsilon in eV (@ 1K) [Halicioglu&Pound]
//    ljp[1] = 2.62; //sigma in A [Halicioglu&Pound]

    ljp[0] = 1; //epsilon in eV 
    ljp[1] = 1; //sigma in A 

    //double er2 = Objective_Function2(ljp);
    //cout << "\nAverage error square: " << er2 << endl;

    find_min(ljp);

/*

	// Lattice constant and atomic mass is the only input
    double lat = 4.05, mass = 26.98; 
	int nsim = 30; // # of simulation for each crystal structure & size
    vector<vector <double> > data;
    //data.resize(natoms,vector<double>(3));

	srand(time(0));

    for (int i=0; i<nsim; i++)    {

        Create_DFT_Confs(data,i,lat,mass);

    }
 
*/


}











