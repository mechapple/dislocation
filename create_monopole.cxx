/*
 * create_monopole.cxx
 * 
 * Copyright 2017 Anirban <anirban@ZeroPointEnergy>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * USAGE:
 * ./a.out <glide plane index(1 or 2 for b1 and b2 respectively)> <disl_type=(e or s for edge/screw)> <disl_sign=(+/-1 for +/- dislocation)>  <nx> <ny> <nz>
 */


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <vector>
#include <stdlib.h>
#include <algorithm>
#include <cmath>

#define NUM_TYPE_LINES 5
#define NUM_BOUND_LINES 4

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

std::vector<std::string> split_into_words(std::string line) {
	std::stringstream ss(line);
	std::istream_iterator<std::string> begin(ss);
	std::istream_iterator<std::string> end;
	std::vector<std::string> vstrings(begin, end);
	return vstrings;
}

std::vector<double> split_into_doubles(std::string line) {
	std::stringstream ss(line);
	std::istream_iterator<std::string> begin(ss);
	std::istream_iterator<std::string> end;
	std::vector<std::string> vstrings(begin, end);
	std::vector<double> vnumbers;
	
	for(std::vector<std::string>::iterator it = vstrings.begin(); it != vstrings.end(); ++it)
	{
		std::stringstream geek(*it); double num; geek >> num; vnumbers.push_back(num);
	}
	return vnumbers;
}

std::vector<int> split_into_integers(std::string line) {
	std::stringstream ss(line);
	std::istream_iterator<std::string> begin(ss);
	std::istream_iterator<std::string> end;
	std::vector<std::string> vstrings(begin, end);
	std::vector<int> vnumbers;
	
	for(std::vector<std::string>::iterator it = vstrings.begin(); it != vstrings.end(); ++it)
	{
		std::stringstream geek(*it); int num; geek >> num; vnumbers.push_back(num);
	}
	return vnumbers;
}

class CoMCell {
		double a,b,c,xy,xz,yz;
	public:
		double **mols,**lcom;
		std::vector<double> bounds;
		int num_mols;
		void initialize_arrays();
		void replicate_cell(double **,std::vector<double>,int,int,int,int);
		void add_dislocation(double **,std::vector<double>,int,int,int,int,char,int);
};

void CoMCell::replicate_cell(double **com,std::vector<double> bounds_cell,int nmols,int nx,int ny,int nz){
	int num_cells = nx*ny*nz;
	a = bounds_cell[1]-bounds_cell[0]; b = bounds_cell[3]-bounds_cell[2]; c = bounds_cell[5]-bounds_cell[4];
	xy = bounds_cell[6]; xz = bounds_cell[7]; yz = bounds_cell[8];
	
	bounds = bounds_cell;
	bounds[1] += a*(nx-1); bounds[3] += b*(ny-1); bounds[5] += c*(nz-1);
	bounds[6] = xy*ny; bounds[7] = xz*nz; bounds[8] = yz*nz;
	
	mols = (double **)malloc(nmols*num_cells*sizeof(double *));	for(int i=0;i<nmols*num_cells;i++) mols[i] = (double *)malloc(3*sizeof(double));
	int count_mols=0;
	for(size_t i=0;i<nx;i++) {
		for(size_t k=0;k<nz;k++) {
			for(size_t j=0;j<ny;j++) {
				for(size_t j1=0; j1<nmols; j1++) {
					for(size_t k1=0; k1<3; k1++) 
						mols[count_mols][k1] = com[j1][k1]+(k1==0)*(a*i+xz*k+xy*j)+(k1==1)*(b*j+yz*k)+(k1==2)*(c*k);
					
					count_mols++;
				}
			}
		}
	}
	
	for(int i=0;i<nmols*num_cells;i++) { 
		std::cout << i+1 << " "; 
		for(size_t j=0;j<3;j++) std::cout << mols[i][j] << " "; 
		std::cout << std::endl; 
	}
	
	num_mols = count_mols;
}

void CoMCell::add_dislocation(double **com,std::vector<double> bounds_cell,int nmols,int nx,int ny,int nz,char type, int P){
	int num_cells = nx*ny*nz;
	a = bounds_cell[1]-bounds_cell[0]; b = bounds_cell[3]-bounds_cell[2]; c = bounds_cell[5]-bounds_cell[4];
	xy = bounds_cell[6]; xz = bounds_cell[7]; yz = bounds_cell[8];
	
	lcom = (double **)malloc(nmols*sizeof(double *));
	for(int i=0;i<nmols;i++) {lcom[i] = (double *)malloc(3*sizeof(double)); for(int j=0;j<3;j++) lcom[i][j]=com[i][j];}
	
	for(size_t j1=0; j1<nmols; j1++) { //ensure all mols are above the P plane
		if(lcom[j1][2]<(P-1)*c/nmols) {
				lcom[j1][0] += xz; lcom[j1][1] += yz; lcom[j1][2] += c;
		} 
	}
	
	bounds = bounds_cell;
	bounds[1] += a*(nx-1); bounds[3] += b*(ny-1); bounds[5] += c*(nz-1);
	bounds[6] = xy*ny; bounds[7] = xz*nz; bounds[8] = yz*nz;
	
	mols = (double **)malloc(nmols*num_cells*sizeof(double *));	for(int i=0;i<nmols*num_cells;i++) mols[i] = (double *)malloc(3*sizeof(double));
	int count_mols=0;
	for(size_t i=0;i<nx;i++) {
		for(size_t k=0;k<nz;k++) {
			for(size_t j=0;j<ny;j++) {
				for(size_t j1=0; j1<nmols; j1++) {
					for(size_t k1=0; k1<3; k1++) 
						mols[count_mols][k1] = lcom[j1][k1]+(k1==0)*(a*i+xz*k+xy*j)+(k1==1)*(b*j+yz*k)+(k1==2)*(c*k);
					
					//apply field
					count_mols++;
				}
			}
		}
	}
	
	for(int i=0;i<nmols*num_cells;i++) { 
		std::cout << i+1 << " "; 
		for(size_t j=0;j<3;j++) std::cout << mols[i][j] << " "; 
		std::cout << std::endl; 
	}
	
	num_mols = count_mols;
}

class UnitCell {
		double **atoms;
		int ***topology; //bonds, angles, dihedrals, impropers
	public:
		std::vector<int> members;
		double **com; int num_mols;
		void read_cell(char *);
		void initialize_arrays();
		void write_cell(char *);
		void write_headers(char *,char *, int);
		void write_bounds(char *,std::vector<double>,int);
		void write_atoms(char *, double **, int);
		void write_coeffs(char *);
		void write_topology(char *, int);
		std::vector<std::string> types,coeffs;
		std::vector<int> num_types,num_entries;
		std::vector<double> bounds,mass;
};

void UnitCell::initialize_arrays()
{
	atoms = (double **)malloc(num_entries[0]*sizeof(double *));		for(int i=0;i<num_entries[0];i++) atoms[i] = (double *)malloc(7*sizeof(double));
	int mem[] = {4,5,6,6};
	members.assign(mem,mem+4); // 4 entries in bonds, 5 entries in angles, etc
	
	topology = (int ***)malloc(4*sizeof(int **));
	for(size_t k=0; k<4; k++) {
		topology[k] = (int **)malloc(num_entries[k+1]*sizeof(int *)); for(int i=0;i<num_entries[k+1];i++) topology[k][i] = (int *)malloc(members[k]*sizeof(int));
	}
}

void UnitCell::read_cell(char *input) {
	
	std::ifstream inputfile;
	std::string line;
	inputfile.open(input);
	if(inputfile.is_open())
	{
		for(size_t i=0; i <2; i++) getline(inputfile, line);
		
		for(size_t i=0; i <NUM_TYPE_LINES; i++) {
			getline(inputfile, line); std::vector<std::string> words = split_into_words(line); 
			std::stringstream geek(words[0]); int num; geek >> num; num_entries.push_back(num);
		}
		
		getline(inputfile, line);
		for(size_t i=0; i <NUM_TYPE_LINES; i++) {
			getline(inputfile, line); types.push_back(line); std::vector<std::string> words = split_into_words(line); 
			std::stringstream geek(words[0]); int num; geek >> num; num_types.push_back(num);
		}
		
		getline(inputfile, line); getline(inputfile, line);
				
		for(size_t i=0; i<NUM_BOUND_LINES; i++) {
			getline(inputfile, line); std::vector<std::string> words = split_into_words(line); 
			{ std::stringstream geek(words[0]); double num; geek >> num; bounds.push_back(num); }
			{ std::stringstream geek(words[1]); double num; geek >> num; bounds.push_back(num); }
			if (i==(NUM_BOUND_LINES-1)) { std::stringstream geek(words[2]); double num; geek >> num; bounds.push_back(num); }
		}
		
		getline(inputfile, line);
		int num_coeffs = 0; for(std::vector<int>::iterator it = num_types.begin(); it != num_types.end(); ++it) num_coeffs += *it;
		
		for(size_t i=0; i <(num_coeffs+3*num_types.size()); i++) {
			getline(inputfile, line); coeffs.push_back(line); 
			if(i>1 && i<6) {
				std::vector<std::string> words = split_into_words(line); 
				std::stringstream geek1(words[1]); double num; geek1 >> num; mass.push_back(num);
			}	
		}
		
		initialize_arrays();
		
		num_mols=0; double total_mass = 0;
		//READING ATOMS
		getline(inputfile, line); getline(inputfile, line);
		for(size_t i=0; i<num_entries[0]; i++) {
			getline(inputfile, line); std::vector<double> nwords = split_into_doubles(line); 
			for(size_t j=0; j<nwords.size(); j++) atoms[i][j] = nwords[j];
			num_mols = (num_mols < atoms[i][1])? atoms[i][1] : num_mols;
			total_mass += mass[(atoms[i][2]-1)];
		}
		getline(inputfile, line);
		
		//CALCULATING COMs
		com = (double **)malloc(num_mols*sizeof(double *));	double mol_mass = total_mass/num_mols;
		for(int i=0;i<num_mols;i++) {com[i] = (double *)malloc(3*sizeof(double)); for(int j=0;j<3;j++) com[i][j]=0.0;}
		for(size_t i=0; i<num_entries[0]; i++) {
			int mol_id = atoms[i][1], type = atoms[i][2];
			for(size_t j=0; j<3; j++) com[mol_id-1][j] += atoms[i][j+4]*mass[type-1]/mol_mass;
		}
		
		//READING TOPOLOGY
		for(size_t k=0; k<4; k++) {
			getline(inputfile, line); getline(inputfile, line);
			for(size_t i=0; i<num_entries[k+1]; i++) {
				getline(inputfile, line); std::vector<int> nwords = split_into_integers(line); 
				for(size_t j=0; j<nwords.size(); j++) topology[k][i][j] = nwords[j];
				//for(size_t j=0; j<nwords.size(); j++) std::cout << topology[k][i][j] << " "; std::cout << std::endl;
			}
			getline(inputfile, line);
		}
	}
	
	inputfile.close();
}

void UnitCell::write_headers(char *outfile, char *infile, int num_cells) {
	std::ofstream outputfile(outfile, std::ios::out | std::ios::app);
	std::cout << "Writing headers\n" ;
	
	outputfile << "LAMMPS Description " << infile << " " << num_cells << "\n";
	std::string str[] = { " atoms", " bonds", " angles", " dihedrals", " impropers" };
	
	outputfile << std::endl;
	for(size_t i=0; i<num_entries.size(); i++) {
		outputfile << " " << num_entries[i]*num_cells << str[i] << std::endl;
	}
	
	outputfile << std::endl;
	for(std::vector<std::string>::iterator it = types.begin(); it != types.end(); ++it) {
		outputfile << *it << std::endl;
	}
	outputfile.close();
}

void UnitCell::write_bounds(char *outfile, std::vector<double> bounds2, int tilt_correct) {
	std::ofstream outputfile(outfile, std::ios::out | std::ios::app);
	outputfile << "\n# Cell: triclinic" << std::endl;
	
	// CORRECT FOR EXTREME TILTS > 0.5*BOX LENGTH (LAMMPS REQUIREMENT)
	if(tilt_correct==1)
	{
		//std::cout << "\nCorrecting tilt factors\n" ;
		double lx = bounds2[1]-bounds2[0], ly = bounds2[3]-bounds2[2], lz = bounds2[5]-bounds2[4];
		while(fabs(bounds2[6])>lx*0.5) bounds2[6]-=sgn(bounds2[6])*lx; 
		while(fabs(bounds2[7])>lx*0.5) bounds2[7]-=sgn(bounds2[7])*lx; 
		while(fabs(bounds2[8])>ly*0.5) bounds2[8]-=sgn(bounds2[8])*ly; 
	}
	
	outputfile << " " << bounds2[0] << " " << bounds2[1] << " xlo xhi" << std::endl;
	outputfile << " " << bounds2[2] << " " << bounds2[3] << " ylo yhi" << std::endl;
	outputfile << " " << bounds2[4] << " " << bounds2[5] << " zlo zhi" << std::endl;
	outputfile << " " << bounds2[6] << " " << bounds2[7] << " " << bounds2[8] << " xy xz yz" << std::endl << std::endl;
	
	outputfile.close();
}

void UnitCell::write_coeffs(char *outfile) {
	std::ofstream outputfile(outfile, std::ios::out | std::ios::app);
		
	for(std::vector<std::string>::iterator it = coeffs.begin(); it != coeffs.end(); ++it) {
		outputfile << *it << std::endl;
	}
	
	outputfile.close();
}

void UnitCell::write_atoms(char *outfile, double **mols, int nmols) {
	std::ofstream outputfile(outfile, std::ios::out | std::ios::app);
	std::cout << "Writing atoms\n" ;
	//std::cout << "Num_cells " << nmols/num_mols << " Num_mols " <<  nmols;
	
	outputfile << "Atoms\n";
	for(size_t i=0; i<nmols/num_mols; i++) {
		for(size_t j=0; j<num_entries[0]; j++) {
			
			int mol_id = atoms[j][1], type = atoms[j][2]; double disp[3];
			for(size_t k=0; k<3; k++) disp[k] = mols[i*num_mols+mol_id-1][k]-com[mol_id-1][k];
			
			outputfile << std::endl << j+i*num_entries[0]+1 << " " << atoms[j][1] << " " << atoms[j][2] << " " << atoms[j][3];
			outputfile << " " << atoms[j][4]+disp[0] << " " << atoms[j][5]+disp[1] << " " << atoms[j][6]+disp[2];
		}
	}
	outputfile.close();
}

void UnitCell::write_topology(char *outfile, int num_cells) {
	std::ofstream outputfile(outfile, std::ios::out | std::ios::app);
	std::cout << "Writing topology\n" ;
	
	std::string str[] = { "Bonds", "Angles", "Dihedrals", "Impropers" };
	
	for(size_t k=0; k<4; k++) {
		if(num_entries[k+1]*num_cells>0) outputfile << "\n\n" << str[k] << "\n";
		
		for(size_t m=0; m<num_cells; m++) {
		
			for(size_t i=0; i<num_entries[k+1]; i++) {
				outputfile << std::endl;
				for(size_t j=0; j<members[k]; j++) 
					outputfile << topology[k][i][j]+(j<1)*m*num_entries[k+1]+(j>1)*m*num_entries[0] << " ";
			}
		}
	}

	outputfile.close();
}

void UnitCell::write_cell(char *outfile) {
	std::ofstream outputfile(outfile, std::ios::out | std::ios::app);
		

	outputfile.close();
}

int main(int argc, char **argv)
{
	UnitCell cell;
	CoMCell cell_mols;
	cell.read_cell(argv[1]);
	
	int nx=10, ny=5, nz=10;
	int num_cells = nx*ny*nz;
	char type='s'; int P = 1;
	
	
	//cell_mols.replicate_cell(cell.com,cell.bounds,cell.num_mols,nx,ny,nz);
	cell_mols.add_dislocation(cell.com,cell.bounds,cell.num_mols,nx,ny,nz,type,P);
	
	char outfile[]="output.lmp";
	std::ofstream outputfile(outfile); outputfile.close(); //clear contents of output.lmp
	
	cell.write_headers(outfile,argv[1],num_cells);
	cell.write_bounds(outfile, cell_mols.bounds,1); // tilt_correct = 1
	cell.write_coeffs(outfile);
	cell.write_atoms(outfile, cell_mols.mols, cell_mols.num_mols);
	cell.write_topology(outfile,num_cells); //num_cells >= 1
	
	return 0;
}

