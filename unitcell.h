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
