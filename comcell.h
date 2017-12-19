typedef struct{ 
	int nmols,nx,ny,nz,P,dip_type,dipole;
	double disl_char,pfrac;
} Configuration;

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
		void add_dislocation(double **,std::vector<double>,int,Configuration);
		void define_cores(Configuration,double [],double [],double [],double []);
		void write_maps(Configuration,double [],double []);
		void write_dumps();
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

void CoMCell::add_dislocation(double **com,std::vector<double> bounds_cell,int nmols,Configuration cfg) {
	
	int P = cfg.P;
	double disl_char = cfg.disl_char;
	int dip_type = cfg.dip_type;
	int nx = cfg.nx;
	int ny = cfg.ny;
	int nz = cfg.nz;
	int dipole = cfg.dipole;
	
	double pfrac = (P-1)*1.0/nmols; cfg.pfrac = pfrac;
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
	if(dipole==1 && disl_char==1.0) bounds[1] -= 0.5*a*dip_type;
	reciprocal_vectors();
	
	double coreX[2],coreZ[2],coreX1[2],coreX2[2];
	define_cores(cfg,coreX,coreZ,coreX1,coreX2);
	
	int MAG=IMAGES;
	
	mols = (double **)malloc(nmols*num_cells*sizeof(double *));	for(int i=0;i<nmols*num_cells;i++) mols[i] = (double *)malloc(3*sizeof(double));
	mols0 = (double **)malloc(nmols*num_cells*sizeof(double *)); for(int i=0;i<nmols*num_cells;i++) mols0[i] = (double *)malloc(3*sizeof(double));
	int count_mols=0;
	for(size_t i=0;i<nx;i++) {
		for(size_t k=0;k<nz;k++) {
			
			int cond=0;
			if(dipole==0) cond = (k<=(nz/2-1))*(dip_type==1) + (k>(nz/2-1))*(dip_type==-1);
			if(dipole==1) cond = (k<=(nz/4-1) || k>(nz*3/4-1) )*(dip_type==1) + (k>(nz/4-1) && k<=(nz*3/4-1) )*(dip_type==-1);
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
					for(int k1=-MAG*dipole; k1<=MAG*dipole; k1++) {
					for(int k2=-MAG; k2<=MAG; k2++)	{
						//DISLOCATION 1
						
							///// 1st partial
							z = mols[count_mols][2] - (coreZ[0] + c*nz*k1);
							x = mols[count_mols][0] - (coreX1[0] + a*nx*k2 + xz*nz*k1 + mols[count_mols][1]*bounds[6]/(bounds[3]-bounds[2]) );
							
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
							x = mols[count_mols][0] - (coreX1[1] + a*nx*k2 + xz*nz*k1 + mols[count_mols][1]*bounds[6]/(bounds[3]-bounds[2]) );
							lam = (1-v)*(x*x + z*z);
							com_dy += (1.0-disl_char)*(burger/(2*PI))*atan2(z,x)*b/bmag;
							com_dz += 0.0;
							com_dx += (1.0-disl_char)*(burger/(2*PI))*atan2(z,x)*xy/bmag;
							
							com_dz += disl_char*(-burger/(2*PI))*( ((1-2*v)/(4*(1-v)))*log(x*x + z*z) + ((x*x-z*z)/(4*lam)) );
							com_dy += 0.0;
							com_dx += disl_char*(burger/(2*PI))*( atan2(z,x) + (x*z/(2*lam)) );
						
						
						if(dipole==1) //DISLOCATION 2 if dipole
						{
							///// 1st partial
							z = mols[count_mols][2] - (coreZ[1] + c*nz*k1);
							x = mols[count_mols][0] - (coreX2[0] + a*nx*k2 + xz*nz*k1 + mols[count_mols][1]*bounds[6]/(bounds[3]-bounds[2]) );
							
							double bmag = sqrt(b*b+xy*xy);
							double burger = -0.5*dip_type*bmag*(1.0-disl_char) - 0.5*dip_type*a*disl_char;
							double lam = (1-v)*(x*x + z*z);
							
							//std::cout << "\n Burgers " << burger << "\n";
							com_dy += (1.0-disl_char)*(burger/(2*PI))*atan2(z,x)*b/bmag;
							com_dz += 0.0;
							com_dx += (1.0-disl_char)*(burger/(2*PI))*atan2(z,x)*xy/bmag;
							
							com_dz += disl_char*(-burger/(2*PI))*( ((1-2*v)/(4*(1-v)))*log(x*x + z*z) + ((x*x-z*z)/(4*lam)) );
							com_dy += 0.0;
							com_dx += disl_char*(burger/(2*PI))*( atan2(z,x) + (x*z/(2*lam)) );
							
							///// 2nd partial
							x = mols[count_mols][0] - (coreX2[1] + a*nx*k2 + xz*nz*k1 + mols[count_mols][1]*bounds[6]/(bounds[3]-bounds[2]) );
							lam = (1-v)*(x*x + z*z);
							com_dy += (1.0-disl_char)*(burger/(2*PI))*atan2(z,x)*b/bmag;
							com_dz += 0.0;
							com_dx += (1.0-disl_char)*(burger/(2*PI))*atan2(z,x)*xy/bmag;
							
							com_dz += disl_char*(-burger/(2*PI))*( ((1-2*v)/(4*(1-v)))*log(x*x + z*z) + ((x*x-z*z)/(4*lam)) );
							com_dy += 0.0;
							com_dx += disl_char*(burger/(2*PI))*( atan2(z,x) + (x*z/(2*lam)) );
						}
					
					}
					}
				
					if(dipole==0) {
						com_dy -= 0.5*MAG*b*sgn(com_dy)*(1.0-disl_char);
						com_dx -= 0.5*MAG*xy*sgn(com_dy)*(1.0-disl_char);
						if((mols[count_mols][2]>coreZ[0])) com_dx -= MAG*a*dip_type*disl_char;
					}
					
					if(dipole==1) {
						if((mols[count_mols][2]<coreZ[1])&&(mols[count_mols][2]>coreZ[0])) {
							com_dy -= 1.0*MAG*b*dip_type*(1.0-disl_char);
							com_dx -= 1.0*MAG*xy*dip_type*(1.0-disl_char);
							com_dx -= 1.0*MAG*a*dip_type*disl_char;
						}
					}
					
					mols[count_mols][0] += com_dx; mols[count_mols][1] += com_dy; mols[count_mols][2] += com_dz;
					Vector disp = {com_dx,com_dy,com_dz}; printVector(disp);
					count_mols++;
				}
			}
		}
	}
	num_mols = count_mols;
	//std::cout << std::endl;
	write_maps(cfg,coreX,coreZ);
	
	for(int i=0;i<3; i++) for(int k=0;k<3;k++) Rotation[i][k] = (i==k);
	move_into_box();
	tilt_correct();
	
	if(dipole==0) {bounds[5] += VACUUM*0.5; bounds[4] -= VACUUM*0.5;} // add vacuum padding for monopoles
}

void CoMCell::write_maps(Configuration cfg,double coreX[],double coreZ[]) {
	int **planeL,**planeU,**planeS,***mapLU;
	int count_S[2]={0,0},count_U[2]={0,0},count_L[2]={0,0},count_mid[2]={0,0};
	int num_cores = 1+cfg.dipole;
	double dis_reg[2] = {a*0.5,a*0.5};
	
	FILE *fcomg,*map[2],*fcoms;
	map[0] = fopen("COM_map1.txt","w"); if(cfg.dipole==1) map[1] = fopen("COM_map2.txt","w");
	fcomg = fopen("Fixed_molecules_glide.RDX","w");
	fcoms = fopen("Fixed_molecules_surf.RDX","w");
	
	planeS = (int **)malloc(2*sizeof(int *)); for(int i=0;i<2;i++) planeS[i] = (int *)malloc(2*cfg.nx*cfg.nz*sizeof(int));
	planeU = (int **)malloc(num_cores*sizeof(int *)); for(int i=0;i<num_cores;i++) planeU[i] = (int *)malloc(2*cfg.nx*cfg.nz*sizeof(int));
	planeL = (int **)malloc(num_cores*sizeof(int *)); for(int i=0;i<num_cores;i++) planeL[i] = (int *)malloc(2*cfg.nx*cfg.nz*sizeof(int));
	
	mapLU = (int ***)malloc(num_cores*sizeof(int **)); for(int i=0;i<num_cores;i++) { 
		mapLU[i] = (int **)malloc(4*cfg.nx*cfg.ny*sizeof(int *));
		for(int j=0;j<4*cfg.nx*cfg.ny;j++) mapLU[i][j] = (int *)malloc(2*sizeof(int));
	}
	
	for(int i=0;i<num_mols;i++)	{
		for(int j=0;j<num_cores;j++) {
			
			if(abs(mols0[i][2]-coreZ[j])<1.1*c && ((mols0[i][1]-yz*(int)(mols0[i][2]/c))<(b+yz)) ) {
				if(mols0[i][2]>coreZ[j]) planeU[j][count_U[j]++] = i;
				else planeL[j][count_L[j]++] = i;
			}		
		}
		for(int j=0;j<2;j++) {if(abs(mols0[i][2]-cfg.nz*c*j)<1.1*c) planeS[j][count_S[j]++] = i;}
	}
	for(int j=0;j<2;j++) printf("count_S%d %d ",j,count_S[j]);
	for(int j=0;j<num_cores;j++) printf("count_U%d %d ",j,count_U[j]);
	for(int j=0;j<num_cores;j++) printf("count_L%d %d ",j,count_L[j]);
	printf("\n");
	
	char gname[20];
	for(int k=0;k<2;k++) {
		if(k==0) strcpy(gname, "bottom\0"); else strcpy(gname, "top\0");
		for(int i=0;i<count_S[k];i++) {
			if(i%8==0) fprintf(fcoms,"\ngroup %s molecule",gname);
			fprintf(fcoms," %d",planeS[k][i]+1);
		}
		
	}
		
	for(int k=0;k<num_cores;k++) {
		for(int i=0;i<count_L[k];i++) for(int j=0;j<count_U[k];j++) 
		{
			double com_dx = (mols0[planeU[k][j]][0] - mols0[planeL[k][i]][0])*((xz>=0.0)-(xz<0.0)); 
			if(com_dx>=0 && com_dx<dis_reg[k]) dis_reg[k]=com_dx; 
		} dis_reg[k]+=a*0.01;
		printf("disreg%d %lf ",k,dis_reg[k]);
	} printf("\n");
	
	for(int k=0;k<num_cores;k++) {
		printf("upper2 %d %d\n",planeU[k][0]+1,planeU[k][1]+1);
		printf("lower2 %d %d\n",planeL[k][0]+1,planeL[k][1]+1);
		
		for(int i=0;i<count_L[k];i++) {
			for(int j=0;j<count_U[k];j++) {
			
				double com_dx = (mols0[planeU[k][j]][0] - mols0[planeL[k][i]][0])*((xz>=0.0)-(xz<0.0));
				double com_dy = (mols0[planeU[k][j]][1] - mols0[planeL[k][i]][1])*((yz>=0.0)-(yz<0.0));
				double com_dz = mols0[planeU[k][j]][2] - mols0[planeL[k][i]][2];
				
				if(com_dx<=dis_reg[k] && com_dx>=0.0) {
					mapLU[k][count_mid[k]][0] = planeL[k][i];
					mapLU[k][count_mid[k]][1] = planeU[k][j];
				
					fprintf(map[k],"%d %d\n",mapLU[k][count_mid[k]][0]+1,mapLU[k][count_mid[k]][1]+1);
					
					if(count_mid[k]%4==0) fprintf(fcomg,"\ngroup mid%d molecule ",k+1);
					fprintf(fcomg,"%d %d ",mapLU[k][count_mid[k]][0]+1,mapLU[k][count_mid[k]][1]+1);
					count_mid[k]++;
				}
			}
		}
		
		if(cfg.disl_char==1.0) {
			
			if(cfg.dip_type*sgn(0.5-k*1.0)==1) fprintf(fcomg,"%d %d ",planeU[k][0]+1,planeU[k][1]+1);
			if(cfg.dip_type*sgn(0.5-k*1.0)==-1) fprintf(fcomg,"%d %d ",planeL[k][0]+1,planeL[k][1]+1);
			
			for(int i=count_mid[k];i<2*count_mid[k];i++) {
				if(cfg.dip_type*sgn(0.5-k*1.0)==-1) {
					mapLU[k][i][0] = planeL[k][(i-count_mid[k])%2];
					if((i-count_mid[k])>=2) mapLU[k][i][0] = mapLU[k][i-count_mid[k]-2][0];
					mapLU[k][i][1] = mapLU[k][i-count_mid[k]][1];
				}
				if(cfg.dip_type*sgn(0.5-k*1.0)==1) {
					mapLU[k][i][1] = planeU[k][(i-count_mid[k])%2];
					if((i-count_mid[k])>=2) mapLU[k][i][1] = mapLU[k][i-count_mid[k]-2][1];
					mapLU[k][i][0] = mapLU[k][i-count_mid[k]][0];
				}
				
				fprintf(map[k],"%d %d\n",mapLU[k][i][0]+1,mapLU[k][i][1]+1);
			}
		}
	}
	
	fclose(fcomg); fclose(fcoms); fclose(map[0]);
	if(cfg.dipole==1) fclose(map[1]);
}

void CoMCell::define_cores(Configuration cfg,double coreX[],double coreZ[],double coreX1[],double coreX2[]) {
	if(cfg.dipole==0) {
		coreZ[0] = (cfg.nz/2)*c + cfg.pfrac*c; coreX[0] = a*cfg.nx*0.5+coreZ[0]*bounds[7]/(bounds[5]-bounds[4]);
		coreX1[0] = coreX[0]-HALF_PARTIAL_SEPARATION*a;
		coreX1[1] = coreX[0]+HALF_PARTIAL_SEPARATION*a;
		std::cout << "\n Core " << coreX[0] << " " << coreZ[0] << "\n";
	}
	
	if(cfg.dipole==1) {
		double centerX = a*cfg.nx*0.5+((cfg.nz*0.5+cfg.pfrac)*c)*bounds[7]/(bounds[5]-bounds[4]);
		if(cfg.disl_char==0.0) { // screw
			coreZ[0] = (cfg.nz/4)*c + cfg.pfrac*c;
			coreZ[1] = (cfg.nz*3/4)*c + cfg.pfrac*c;
			coreX[0] = centerX;
			coreX[1] = centerX;
		}
		if(cfg.disl_char==1.0) { // edge
			coreZ[0] = (cfg.nz/4)*c + cfg.pfrac*c;
			coreZ[1] = (cfg.nz*3/4)*c + cfg.pfrac*c;
			coreX[0] = centerX - ((xz>=0.0)-(xz<0.0))*(coreZ[1]-coreZ[0])*0.5;;
			coreX[1] = centerX + ((xz>=0.0)-(xz<0.0))*(coreZ[1]-coreZ[0])*0.5;;
		}
		
		coreX1[0] = coreX[0]-HALF_PARTIAL_SEPARATION*a;
		coreX1[1] = coreX[0]+HALF_PARTIAL_SEPARATION*a;
		
		coreX2[0] = coreX[1]-HALF_PARTIAL_SEPARATION*a;
		coreX2[1] = coreX[1]+HALF_PARTIAL_SEPARATION*a;
		
		std::cout << "Core1 " << coreX[0] << " " << coreZ[0] << "\n";
		std::cout << "Core2 " << coreX[1] << " " << coreZ[1] << "\n";
	}
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
	
	printf("Rotation matrix :\n");
	for(int i=0;i<3; i++) {for(int k=0;k<3;k++) printf(" %lf", Rotation[i][k]); printf("\n");}
	
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

void CoMCell::write_dumps() {
	FILE *fcom = fopen("COM.lammpstrj","w");
	fprintf(fcom,"ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS xy xz yz pp pp pp",num_mols);
	fprintf(fcom,"\n%.3f %.3f %.3f\n%.3f %.3f %.3f\n%.3f %.3f %.3f\nITEM: ATOMS id mol type q xu yu zu\n", \
				 bounds[0] + std::min({0.0, bounds[6], bounds[7], bounds[6]+bounds[7]}), bounds[1] + std::max({0.0, bounds[6], bounds[7], bounds[6]+bounds[7]}), bounds[6], \
				 bounds[2] + std::min({0.0, bounds[8]}),bounds[3] + std::max({0.0, bounds[8]}),bounds[7], \
				 bounds[4],bounds[5],bounds[8]);
	for(int i=0;i<num_mols;i++)	{
		fprintf(fcom,"%d 1 1 0.00 %lf %lf %lf\n",i+1,mols0[i][0],mols0[i][1],mols0[i][2]);
	}
	
	fprintf(fcom,"ITEM: TIMESTEP\n1\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS xy xz yz pp pp pp",num_mols);
	fprintf(fcom,"\n%.3f %.3f %.3f\n%.3f %.3f %.3f\n%.3f %.3f %.3f\nITEM: ATOMS id mol type q xu yu zu\n", \
				 bounds[0] + std::min({0.0, bounds[6], bounds[7], bounds[6]+bounds[7]}), bounds[1] + std::max({0.0, bounds[6], bounds[7], bounds[6]+bounds[7]}), bounds[6], \
				 bounds[2] + std::min({0.0, bounds[8]}),bounds[3] + std::max({0.0, bounds[8]}),bounds[7], \
				 bounds[4],bounds[5],bounds[8]);
	for(int i=0;i<num_mols;i++)	{
		fprintf(fcom,"%d 1 1 0.00 %lf %lf %lf\n",i+1,mols[i][0],mols[i][1],mols[i][2]);
	}
	
	fclose(fcom);
}
