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
 * COMPILE:
 * g++ -std=c++11 -o create_monopole create_monopole.cxx
 * USAGE:
 * ./a.out <data_file> <glide plane index(1 or 2 for b1 and b2 respectively)> <disl_char=(1/0 for edge/screw)> <disl_sign=(+/-1 for +/- dislocation)>  <nx> <ny> <nz> <pole=(1/0 for dipole/monopole)>
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
#include <string.h>

#define NUM_TYPE_LINES 5
#define NUM_BOUND_LINES 4
#define HALF_PARTIAL_SEPARATION 2.0
#define IMAGES 10
#define v 0.1
#define PI 3.14159265
#define VACUUM 100.0

template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }
double inbox(double x0) { return x0 - floor(x0); }

#include <vector.h>
#include <rotate.h>
#include <split.h>
#include <comcell.h>
#include <unitcell.h>

int main(int argc, char **argv)
{
	UnitCell cell;
	CoMCell cell_mols;
	Configuration cfg;
	cell.read_cell(argv[1]);
	
	cfg.nx=20; cfg.ny=1; cfg.nz=40;
	int num_cells;
	cfg.P = 1; cfg.disl_char = 0.0; cfg.dip_type = 1; cfg.dipole = 0;
	
	if(argc==9) {
		cfg.P = atoi(argv[2]);
		cfg.disl_char = atof(argv[3]);
		cfg.dip_type = atoi(argv[4]);
		cfg.nx = atoi(argv[5]);
		cfg.ny = atoi(argv[6]);
		cfg.nz = atoi(argv[7]);
		cfg.dipole = atoi(argv[8]);
	}
	//cell_mols.replicate_cell(cell.com,cell.bounds,cell.num_mols,nx,ny,nz);
	cell_mols.add_dislocation(cell.com,cell.bounds,cell.num_mols,cfg);
	num_cells = cell_mols.num_mols/2;
	cell_mols.write_dumps();
	
	char outfile[]="output.lmp"; char intro[100];
	std::ofstream outputfile(outfile); outputfile.close(); //clear contents of output.lmp
	
	sprintf(intro,"%s P%d dc%.0f dt%d nx=%d ny=%d nz=%d cells=%d",argv[1],cfg.P,cfg.disl_char,cfg.dip_type,cfg.nx,cfg.ny,cfg.nz,num_cells);
	
	cell.write_headers(outfile,intro,num_cells);
	cell.write_bounds(outfile, cell_mols.bounds);
	cell.write_coeffs(outfile);
	cell.write_atoms(outfile, cell_mols.mols, num_cells, cell_mols.Rotation);
	cell.write_topology(outfile,num_cells); //num_cells >= 1
	
	return 0;
}

