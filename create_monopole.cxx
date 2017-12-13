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
 * ./a.out <data_file> <glide plane index(1 or 2 for b1 and b2 respectively)> <disl_char=(1/0 for edge/screw)> <disl_sign=(+/-1 for +/- dislocation)>  <nx> <ny> <nz>
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
#define HALF_PARTIAL_SEPARATION 2.0
#define IMAGES 10
#define v 0.1
#define PI 3.14159265
#define VACUUM 100.0

template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }
double inbox(double x0) { return x0 - floor(x0); }

typedef struct{ double i,j,k; } Vector;
Vector crossProduct(Vector a,Vector b) { Vector c = {a.j*b.k - a.k*b.j, a.k*b.i - a.i*b.k, a.i*b.j - a.j*b.i};	return c; }
double dotProduct(Vector a, Vector b) {	return a.i*b.i+a.j*b.j+a.k*b.k; }
void printVector(Vector a) {std::cout << a.i << " " << a.j << " " << a.k << std::endl;}

void rotate(double R[3][3],double vec[3])
{
	double vec2[3];
	vec2[0] = R[0][0]*vec[0] + R[0][1]*vec[1] + R[0][2]*vec[2];
	vec2[1] = R[1][0]*vec[0] + R[1][1]*vec[1] + R[1][2]*vec[2];
	vec2[2] = R[2][0]*vec[0] + R[2][1]*vec[1] + R[2][2]*vec[2];
	
	vec[0] = vec2[0]; vec[1] = vec2[1]; vec[2] = vec2[2];
}

double angle_dot(double x1,double x2,double x3,double y1,double y2,double y3)
{
	double xmag = sqrt(x1*x1+x2*x2+x3*x3);
	double ymag = sqrt(y1*y1+y2*y2+y3*y3);
	
	double dot = x1*y1+x2*y2+x3*y3;
	double theta = acos(dot/(xmag*ymag));
	
	return theta;
}

void unit_cross(double x1,double x2,double x3,double y1,double y2,double y3,double n[3])
{
	n[0] = (x2*y3-x3*y2);
	n[1] = (x3*y1-x1*y3);
	n[2] = (x1*y2-x2*y1);
	
	double mag = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
	n[0] = n[0]/mag; n[1] = n[1]/mag; n[2] = n[2]/mag;
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
		double **mols,**lcom,**mols0;
		Vector av,bv,cv,rav,rbv,rcv; double box_vol;
		std::vector<double> bounds;
		int num_mols;
		double Rotation[3][3];
		//void initialize_arrays();
		void reciprocal_vectors();
		Vector real_to_fractional(double,double,double);
		Vector fractional_to_real(double,double,double);
		void tilt_correct();
		void replicate_cell(double **,std::vector<double>,int,int,int,int);
		void add_dislocation(double **,std::vector<double>,int,int,int,int,double,int,int);
		void move_into_box();
		void wrap_into_box();

};

void CoMCell::reciprocal_vectors() {
	av = {bounds[1]-bounds[0],0,0};
	bv = {bounds[6],bounds[3]-bounds[2],0};
	cv = {bounds[7],bounds[8],bounds[5]-bounds[4]};
	
	rav = crossProduct(bv,cv), rbv = crossProduct(cv,av), rcv = crossProduct(av,bv);
	box_vol = dotProduct(av,rav);
	std::cout << "\nVolume " << box_vol << "\n";
}

Vector CoMCell::real_to_fractional(double x,double y,double z) {
	Vector vec = {x,y,z};
	Vector frac = { dotProduct(vec,rav)/box_vol, dotProduct(vec,rbv)/box_vol, dotProduct(vec,rcv)/box_vol};
	return frac;
}

Vector CoMCell::fractional_to_real(double x,double y,double z) {
	Vector vec = {x,y,z};
	Vector real = { x*av.i+y*bv.i+z*cv.i, x*av.j+y*bv.j+z*cv.j, x*av.k+y*bv.k+z*cv.k};
	return real;
}

void CoMCell::tilt_correct() {
	double lx = bounds[1]-bounds[0], ly = bounds[3]-bounds[2], lz = bounds[5]-bounds[4];
	while(fabs(bounds[6])>lx*0.5) bounds[6]-=sgn(bounds[6])*lx; 
	while(fabs(bounds[7])>lx*0.5) bounds[7]-=sgn(bounds[7])*lx; 
	while(fabs(bounds[8])>ly*0.5) bounds[8]-=sgn(bounds[8])*ly; 
}

void CoMCell::replicate_cell(double **com,std::vector<double> bounds_cell,int nmols,int nx,int ny,int nz){
	int num_cells = nx*ny*nz;
	a = bounds_cell[1]-bounds_cell[0]; b = bounds_cell[3]-bounds_cell[2]; c = bounds_cell[5]-bounds_cell[4];
	xy = bounds_cell[6]; xz = bounds_cell[7]; yz = bounds_cell[8];
	
	bounds = bounds_cell;
	bounds[1] += a*(nx-1); bounds[3] += b*(ny-1); bounds[5] += c*(nz-1);
	bounds[6] = xy*ny; bounds[7] = xz*nz; bounds[8] = yz*nz;
	
	tilt_correct();
	reciprocal_vectors();
	
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
	num_mols = count_mols;
	
	wrap_into_box();
}

void CoMCell::add_dislocation(double **com,std::vector<double> bounds_cell,int nmols,int nx,int ny,int nz,double disl_char, int P, int dip_type){
	
	double pfrac = (P-1)*1.0/nmols;
	int num_cells = nx*ny*nz;
	a = bounds_cell[1]-bounds_cell[0]; b = bounds_cell[3]-bounds_cell[2]; c = bounds_cell[5]-bounds_cell[4];
	xy = bounds_cell[6]; xz = bounds_cell[7]; yz = bounds_cell[8];
	
	lcom = (double **)malloc(nmols*sizeof(double *));
	for(int i=0;i<nmols;i++) {lcom[i] = (double *)malloc(3*sizeof(double)); for(int j=0;j<3;j++) lcom[i][j]=com[i][j];}
	
	//ensure all mol coms are above the P plane
	for(size_t j1=0; j1<nmols; j1++) { 
		if(lcom[j1][2]<pfrac*c) {
				lcom[j1][0] += xz; lcom[j1][1] += yz; lcom[j1][2] += c;
		} 
	}
	
	bounds = bounds_cell;
	bounds[1] += a*(nx-1); bounds[3] += b*(ny-1); bounds[5] += c*(nz-1);
	bounds[6] = xy*ny; bounds[7] = xz*nz; bounds[8] = yz*nz;
	if(disl_char==1.0) bounds[1] -= 0.5*a;
	
	reciprocal_vectors();
	
	double coreZ = (nz/2)*c + pfrac*c, coreX = a*nx*0.5+coreZ*bounds[7]/(bounds[5]-bounds[4]);
	double coreX1 = coreX-HALF_PARTIAL_SEPARATION*a;
	double coreX2 = coreX+HALF_PARTIAL_SEPARATION*a;
	int MAG=IMAGES;
	std::cout << "\n Core " << coreX << " " << coreZ << "\n";
	
	mols = (double **)malloc(nmols*num_cells*sizeof(double *));	for(int i=0;i<nmols*num_cells;i++) mols[i] = (double *)malloc(3*sizeof(double));
	mols0 = (double **)malloc(nmols*num_cells*sizeof(double *));	for(int i=0;i<nmols*num_cells;i++) mols0[i] = (double *)malloc(3*sizeof(double));
	int count_mols=0;
	for(size_t i=0;i<nx;i++) {
		for(size_t k=0;k<nz;k++) {
			
			int cond;
			if(dip_type==1) cond = (k<=(nz-nz/2-1));
			if(dip_type==-1) cond = (k>(nz-nz/2-1));
			if((i==0)&&(cond==1)&&(disl_char==1.0)) continue;
			
			for(size_t j=0;j<ny;j++) {
				for(size_t j1=0; j1<nmols; j1++) {
					for(size_t k1=0; k1<3; k1++) {
						mols[count_mols][k1] = lcom[j1][k1]+(k1==0)*(a*i+xz*k+xy*j)+(k1==1)*(b*j+yz*k)+(k1==2)*(c*k);
						mols0[count_mols][k1] = mols[count_mols][k1];
					}
					//apply field
					double x,y,z;
					double com_dx = 0,com_dy = 0,com_dz = 0;
					for(int k2=-MAG; k2<=MAG; k2++)
					{
						///// 1st partial
						z = mols[count_mols][2] - coreZ;
						x = mols[count_mols][0] - (coreX1 + a*nx*k2 + mols[count_mols][1]*bounds[6]/(bounds[3]-bounds[2]) );
						
						double bmag = sqrt(b*b+xy*xy);
						double burger = 0.5*dip_type*bmag*(1.0-disl_char) + 0.5*dip_type*a*disl_char;
						double lam = (1-v)*(x*x + z*z);
						
						//std::cout << "\n Burgers " << burger << "\n";
						com_dy += (1.0-disl_char)*(burger/(2*PI))*atan2(z,x)*b/bmag;
						com_dz += 0.0;
						com_dx += (1.0-disl_char)*(burger/(2*PI))*atan2(z,x)*xy/bmag;
						
						com_dz += disl_char*(-burger/(2*PI))*( ((1-2*v)/(4*(1-v)))*log(x*x + z*z) + ((x*x-z*z)/(4*lam)) );
						com_dy += 0.0;
						com_dx += disl_char*(burger/(2*PI))*( atan2(z,x) + (x*z/(2*lam)) );
						
						///// 2nd partial
						x = mols[count_mols][0] - (coreX2 + a*nx*k2 + mols[count_mols][1]*bounds[6]/(bounds[3]-bounds[2]) );
						lam = (1-v)*(x*x + z*z);
						com_dy += (1.0-disl_char)*(burger/(2*PI))*atan2(z,x)*b/bmag;
						com_dz += 0.0;
						com_dx += (1.0-disl_char)*(burger/(2*PI))*atan2(z,x)*xy/bmag;
						
						com_dz += disl_char*(-burger/(2*PI))*( ((1-2*v)/(4*(1-v)))*log(x*x + z*z) + ((x*x-z*z)/(4*lam)) );
						com_dy += 0.0;
						com_dx += disl_char*(burger/(2*PI))*( atan2(z,x) + (x*z/(2*lam)) );
					}
					
					com_dy -= 0.5*MAG*b*dip_type*sgn(com_dy)*(1.0-disl_char);
					com_dx -= 0.5*MAG*xy*dip_type*sgn(com_dy)*(1.0-disl_char);
					if((mols[count_mols][2]>coreZ)) com_dx -= MAG*a*dip_type*disl_char;
					
					mols[count_mols][0] += com_dx; mols[count_mols][1] += com_dy; mols[count_mols][2] += com_dz;
					Vector disp = {com_dx,com_dy,com_dz}; //printVector(disp);
					count_mols++;
				}
			}
		}
	}
	num_mols = count_mols;
	//std::cout << std::endl;
	
	move_into_box();
	tilt_correct();
	
	bounds[5] += VACUUM*0.5; bounds[4] -= VACUUM*0.5; // add vacuum padding for monopoles
}

void CoMCell::wrap_into_box() {
	for(int i=0;i<num_mols;i++) { 
		Vector fm = real_to_fractional(mols[i][0]-bounds[0],mols[i][1]-bounds[2],mols[i][2]-bounds[4]);
		Vector rm = fractional_to_real(inbox(fm.i),inbox(fm.j),inbox(fm.k));
		mols[i][0] = bounds[0] + rm.i; mols[i][1] = bounds[2] + rm.j; mols[i][2] = bounds[4] + rm.k;
		std::cout << i+1 << " " << mols[i][0] << " " << mols[i][1] << " " << mols[i][2] << std::endl;
	}
}

void CoMCell::move_into_box() {
	
	int idlz=0,idhz=0;
	for(int i=0;i<num_mols;i++) { 
		if(mols0[i][2]<=mols0[idlz][2] && mols0[i][1]<=b) idlz = i;
		if(mols0[i][2]>=mols0[idhz][2] && mols0[i][1]<=b) idhz = i;
	}
	//printf("idlz idhz %d %d\n",idlz+1,idhz+1);
	
	int idll=idlz,idlr=idlz;
	int idul=idhz,idur=idhz;
	
	for(int i=0;i<num_mols;i++)
	{
		if(abs(mols0[i][2]-mols0[idlz][2])<1.0 && mols0[i][0]<mols0[idll][0] && mols0[i][1]<=b) idll=i;
		if(abs(mols0[i][2]-mols0[idlz][2])<1.0 && mols0[i][0]>mols0[idlr][0] && mols0[i][1]<=b) idlr=i;
		
		if(abs(mols0[i][2]-mols0[idhz][2])<1.0 && mols0[i][0]<mols0[idul][0] && mols0[i][1]<=b) idul=i;
		if(abs(mols0[i][2]-mols0[idhz][2])<1.0 && mols0[i][0]>mols0[idur][0] && mols0[i][1]<=b) idur=i;
	}
	
	//printf("ll lr %d %d\nul ur %d %d\n",idll+1,idlr+1,idul+1,idur+1);
	
	double u[3]; // unit vector for rotation axis
	
	double u0[3] = {mols0[idll][0],mols0[idll][1],mols0[idll][2]}; 
	double u1[3] = {mols0[idlr][0]-mols0[idll][0],mols0[idlr][1]-mols0[idll][1],mols0[idlr][2]-mols0[idll][2]};
	double u2[3] = {mols0[idul][0]-mols0[idll][0],mols0[idul][1]-mols0[idll][1],mols0[idul][2]-mols0[idll][2]};
	
	double v0[3] = {mols[idll][0],mols[idll][1],mols[idll][2]};
	double v1[3] = {mols[idlr][0]-mols[idll][0],mols[idlr][1]-mols[idll][1],mols[idlr][2]-mols[idll][2]};
	double v2[3] = {mols[idul][0]-mols[idll][0],mols[idul][1]-mols[idll][1],mols[idul][2]-mols[idll][2]};
	
	double theta = angle_dot(u1[0],u1[1],u1[2],v1[0],v1[1],v1[2]);
	unit_cross(v1[0],v1[1],v1[2],u1[0],u1[1],u1[2],u);
	
	double R[3][3] = \
	{{cos(theta)+u[0]*u[0]*(1-cos(theta)),      u[0]*u[1]*(1-cos(theta))-u[2]*sin(theta),    u[0]*u[2]*(1-cos(theta))+u[1]*sin(theta)},\
	{u[0]*u[1]*(1-cos(theta))+u[2]*sin(theta),  cos(theta)+u[1]*u[1]*(1-cos(theta)),         u[1]*u[2]*(1-cos(theta))-u[0]*sin(theta)},\
	{u[0]*u[2]*(1-cos(theta))-u[1]*sin(theta),  u[1]*u[2]*(1-cos(theta))+u[0]*sin(theta),    cos(theta)+u[2]*u[2]*(1-cos(theta))}};
		
	double v1r[3] = {v1[0],v1[1],v1[2]};
	double v2r[3] = {v2[0],v2[1],v2[2]};
	rotate(R,v2r); rotate(R,v1r);
	//printf("\n v1r = %lf %lf %lf",v1r[0],v1r[1],v1r[2]);
	//printf("\n v2r = %lf %lf %lf\n",v2r[0],v2r[1],v2r[2]);
	double thetaX = angle_dot(0.0,u2[1],u2[2],0.0,v2r[1],v2r[2]);
	unit_cross(0,v2r[1],v2r[2],0.0,u2[1],u2[2],u);
	
	double Rx[3][3] = \
	{{cos(thetaX)+u[0]*u[0]*(1-cos(thetaX)),      u[0]*u[1]*(1-cos(thetaX))-u[2]*sin(thetaX),    u[0]*u[2]*(1-cos(thetaX))+u[1]*sin(thetaX)},\
	{u[0]*u[1]*(1-cos(thetaX))+u[2]*sin(thetaX),  cos(thetaX)+u[1]*u[1]*(1-cos(thetaX)),         u[1]*u[2]*(1-cos(thetaX))-u[0]*sin(thetaX)},\
	{u[0]*u[2]*(1-cos(thetaX))-u[1]*sin(thetaX),  u[1]*u[2]*(1-cos(thetaX))+u[0]*sin(thetaX),    cos(thetaX)+u[2]*u[2]*(1-cos(thetaX))}};
	
	v2r[0]=0.0;
	rotate(Rx,v2r);
	
	for(int i=0;i<3; i++) for(int k=0;k<3;k++) 
		Rotation[i][k] = Rx[i][0]*R[0][k] + Rx[i][1]*R[1][k] + Rx[i][2]*R[2][k];
	
	//printf("theta %lf thetaX %lf\n",theta*180/PI, thetaX*180/PI);
	
	for(int i=0;i<num_mols;i++)
	{
		double r[3];
		
		for(int j=0;j<3;j++) r[j] = mols[i][j]-1.0*v0[j];
		rotate(Rotation,r);
		for(int j=0;j<3;j++) mols[i][j] = r[j]+1.0*u0[j];
	}
	
	double span0[3] = { 0.5*(mols0[idlr][0]-mols0[idll][0] + mols0[idur][0]-mols0[idul][0]), 
						0.5*(mols0[idul][2]-mols0[idll][2] + mols0[idur][2]-mols0[idlr][2]), 
						0.5*(mols0[idul][0]-mols0[idll][0] + mols0[idur][0]-mols0[idlr][0]) };
	
	double span[3] = { 0.5*(mols[idlr][0]-mols[idll][0] + mols[idur][0]-mols[idul][0]), 
						0.5*(mols[idul][2]-mols[idll][2] + mols[idur][2]-mols[idlr][2]), 
						0.5*(mols[idul][0]-mols[idll][0] + mols[idur][0]-mols[idlr][0]) };
	
	//double span[3] = { mols[idlr][0] - mols[idll][0], mols[idul][2] - mols[idll][2], mols[idul][0] - mols[idll][0] };
	
	printf("Final SpanX = %lf, Absolute change = %lf, Fraction of initial %lf\n",span[0],span[0]-span0[0],span[0]/span0[0]);
	printf("Final SpanZ = %lf, Absolute change = %lf, Fraction of initial %lf\n",span[1],span[1]-span0[1],span[1]/span0[1]);
	printf("Final SpanXZ = %lf, Absolute change = %lf, Fraction of initial %lf\n",span[2],span[2]-span0[2],span[2]/span0[2]);
	
	bounds[1] *= span[0]/span0[0]; bounds[5] *= span[1]/span0[1]; 
	if(fabs(bounds[7])>1e-16) bounds[7] *= span[2]/span0[2];
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
		void write_bounds(char *,std::vector<double>);
		void write_atoms(char *, double **,int,double [3][3]);
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
	
	outputfile << "LAMMPS Description " << infile << "\n";
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

void UnitCell::write_bounds(char *outfile, std::vector<double> bounds2) {
	std::ofstream outputfile(outfile, std::ios::out | std::ios::app);
	outputfile << "\n# Cell: triclinic" << std::endl;
	
	// CORRECT FOR EXTREME TILTS > 0.5*BOX LENGTH (LAMMPS REQUIREMENT)
	//if(tilt_correct==1)
	//{
		////std::cout << "\nCorrecting tilt factors\n" ;
		//double lx = bounds2[1]-bounds2[0], ly = bounds2[3]-bounds2[2], lz = bounds2[5]-bounds2[4];
		//while(fabs(bounds2[6])>lx*0.5) bounds2[6]-=sgn(bounds2[6])*lx; 
		//while(fabs(bounds2[7])>lx*0.5) bounds2[7]-=sgn(bounds2[7])*lx; 
		//while(fabs(bounds2[8])>ly*0.5) bounds2[8]-=sgn(bounds2[8])*ly; 
	//}
	
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

void UnitCell::write_atoms(char *outfile, double **mols, int ncells, double Rotation[3][3]) {
	std::ofstream outputfile(outfile, std::ios::out | std::ios::app);
	std::cout << "Writing atoms\n" ;
	//std::cout << "Num_cells " << nmols/num_mols << " Num_mols " <<  nmols;
	
	outputfile << "Atoms\n";
	for(size_t i=0; i<ncells; i++) {
		for(size_t j=0; j<num_entries[0]; j++) {
			
			int mol_id = atoms[j][1], type = atoms[j][2]; double disp[3];
			for(size_t k=0; k<3; k++) disp[k] = mols[i*num_mols+mol_id-1][k]-com[mol_id-1][k];
			
			outputfile << std::endl << j+i*num_entries[0]+1 << " " << atoms[j][1]+i*2 << " " << atoms[j][2] << " " << atoms[j][3];
			
			double r[3] = {atoms[j][4]-com[mol_id-1][0], atoms[j][5]-com[mol_id-1][1], atoms[j][6]-com[mol_id-1][2]};
			rotate(Rotation,r);
			
			outputfile << " " << r[0]+mols[i*num_mols+mol_id-1][0] << " " << r[1]+mols[i*num_mols+mol_id-1][1] << " " << r[2]+mols[i*num_mols+mol_id-1][2];
			//outputfile << " " << atoms[j][4]+disp[0] << " " << atoms[j][5]+disp[1] << " " << atoms[j][6]+disp[2];
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
	
	int nx=20, ny=1, nz=40;
	int num_cells = nx*ny*nz;
	int P = 1; double disl_char = 0.0; int dip_type = 1;
	
	if(argc==8) {
		P = atoi(argv[2]);
		disl_char = atof(argv[3]);
		dip_type = atoi(argv[4]);
		nx = atoi(argv[5]);
		ny = atoi(argv[6]);
		nz = atoi(argv[7]);
	}
	//cell_mols.replicate_cell(cell.com,cell.bounds,cell.num_mols,nx,ny,nz);
	cell_mols.add_dislocation(cell.com,cell.bounds,cell.num_mols,nx,ny,nz,disl_char,P,dip_type);
	num_cells = cell_mols.num_mols/2;
	
	char outfile[]="output.lmp"; char intro[100];
	std::ofstream outputfile(outfile); outputfile.close(); //clear contents of output.lmp
	
	sprintf(intro,"%s P%d dc%.0f dt%d nx=%d ny=%d nz=%d cells=%d",argv[1],P,disl_char,dip_type,nx,ny,nz,num_cells);
	
	cell.write_headers(outfile,intro,num_cells);
	cell.write_bounds(outfile, cell_mols.bounds);
	cell.write_coeffs(outfile);
	cell.write_atoms(outfile, cell_mols.mols, num_cells, cell_mols.Rotation);
	cell.write_topology(outfile,num_cells); //num_cells >= 1
	
	return 0;
}

