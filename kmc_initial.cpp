#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "kmc_initial.h"
#include "kmc_par.h"

using namespace std;

void class_initial::ltc_constructor(){
	double (*ptr_vbra)[3]; 
	int   (*ptr_v1nbr)[3];
	int   (*ptr_v2nbr)[3];
	int   (*ptr_v3nbr)[3];
    int    *ptr_v1sp;
    int    *ptr_v2sp;
			
	// coordinate vectors of bravais lattices
		
	// BCC
	double vbra_bcc[3][3]= {{-0.5,  0.5,  0.5}, { 0.5, -0.5,  0.5}, { 0.5,  0.5, -0.5}};
	
    int   v1nbr_bcc[8][3]= {{ 1,  0,  0}, { 0,  1,  0}, { 0,  0,  1}, { 1,  1,  1},
                            {-1,  0,  0}, { 0, -1,  0}, { 0,  0, -1}, {-1, -1, -1}};
	
    int   v2nbr_bcc[6][3]= {{ 0,  1,  1}, { 1,  0,  1}, { 1,  1,  0},
                            { 0, -1, -1}, {-1,  0, -1}, {-1, -1,  0}};
    
    int  v3nbr_bcc[12][3]= {{ 2, 1, 1}, { 1, 2, 1}, { 1, 1, 2}, { 1, 0,-1}, { 1,-1, 0}, {0, 1,-1},
                            {-2,-1,-1}, {-1,-2,-1}, {-1,-1,-2}, {-1, 0, 1}, {-1, 1, 0}, {0,-1, 1}};

    int v1sp_bcc[8][6][3]= {{{ 0,-1,-1}, { 0,-1, 0}, { 0, 0,-1}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}},
                            {{-1, 0,-1}, {-1, 0, 0}, { 0, 0,-1}, {0, 1, 1}, {1, 1, 0}, {1, 1, 1}},
                            {{-1,-1, 0}, {-1, 0, 0}, { 0,-1, 0}, {0, 1, 1}, {1, 0, 1}, {1, 1, 1}},
                            {{ 0, 0, 1}, { 0, 1, 0}, { 0, 1, 1}, {1, 0, 0}, {1, 0, 1}, {1, 1, 0}},
                            {{-1,-1,-1}, {-1,-1, 0}, {-1, 0,-1}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}},
                            {{-1,-1,-1}, {-1,-1, 0}, { 0,-1,-1}, {0, 0, 1}, {1, 0, 0}, {1, 0, 1}},
                            {{-1,-1,-1}, {-1, 0,-1}, { 0,-1,-1}, {0, 1, 0}, {1, 0, 0}, {1, 1, 0}},
                            {{-1,-1, 0}, {-1, 0,-1}, {-1, 0, 0}, {0,-1,-1}, {0,-1, 0}, {0, 0,-1}}};

    int v2sp_bcc[8][6][3]= {{{-1,-1,-1}, { 0, 0, 1}, { 0, 1, 0}, {1,-1, 0}, {1, 0,-1}, {2, 1, 1}},
                            {{-1,-1,-1}, {-1, 1, 0}, { 0, 0, 1}, {0, 1,-1}, {1, 0, 0}, {1, 2, 1}},
                            {{-1,-1,-1}, {-1, 0, 1}, { 0,-1, 1}, {0, 1, 0}, {1, 0, 0}, {1, 1, 2}},
                            {{-1, 0, 0}, { 0,-1, 0}, { 0, 0,-1}, {1, 1, 2}, {1, 2, 1}, {2, 1, 1}},
                            {{-2,-1,-1}, {-1, 0, 1}, {-1, 1, 0}, {0,-1, 0}, {0, 0,-1}, {1, 1, 1}},
                            {{-1,-2,-1}, {-1, 0, 0}, { 0,-1, 1}, {0, 0,-1}, {1,-1, 0}, {1, 1, 1}},
                            {{-1,-1,-2}, {-1, 0, 0}, { 0,-1, 0}, {0, 1,-1}, {1, 0,-1}, {1, 1, 1}},
                            {{-2,-1,-1}, {-1,-2,-1}, {-1,-1,-2}, {0, 0, 1}, {0, 1, 0}, {1, 0, 0}}};

    // FCC                   //fh  have the FCC neighbor indexes available
	double vbra_fcc[3][3]= {{0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};

	int   v1nbr_fcc[12][3]= {{ 1,  0,  0}, { 0,  1,  0}, { 0,  0,  1}, { 1, -1,  0}, { 1,  0, -1}, { 0,  1, -1},
                             {-1,  0,  0}, { 0, -1,  0}, { 0,  0, -1}, {-1,  1,  0}, {-1,  0,  1}, { 0, -1,  1}};

	int   v2nbr_fcc[6][3]=  {{ 1,  1, -1}, { 1, -1,  1}, {-1,  1,  1},
                             {-1, -1,  1}, {-1,  1, -1}, { 1, -1, -1}};

	// Choose ltc structure

	ptr_vbra = vbra_bcc;
	n1nbr = 8; ptr_v1nbr = v1nbr_bcc;
	n2nbr = 6; ptr_v2nbr = v2nbr_bcc;
	n3nbr = 12; ptr_v3nbr = v3nbr_bcc;
	n1sp = 6; ptr_v1sp = &v1sp_bcc[0][0][0];
	n2sp = 6; ptr_v2sp = &v2sp_bcc[0][0][0];
	
			
	// assign array values
	for(int i=0; i<3; i ++)      // Bravice vectors
		for(int j=0; j<3; j ++)
			vbra[i][j]= (*(ptr_vbra+i))[j];
    
    n12nbr=  n1nbr+n2nbr;
    n123nbr= n1nbr+n2nbr+n3nbr;
	for(int i=0; i<n1nbr; i ++){ // v1nbr
		for(int j=0; j<3; j ++){
			v1nbr[i][j]= (*(ptr_v1nbr+i))[j];
            if(i>=n1nbr/2) if(v1nbr[i][j] != -v1nbr[i-n1nbr/2][j]) error(1, "(ltc_constructor) v1nbr isn't symmetry");
            v12nbr[i][j]= v1nbr[i][j];
            v123nbr[i][j]= v1nbr[i][j];

            for(int k=0; k<n1sp; k ++) v1sp[i][k][j]= *(ptr_v1sp+i*n1sp*3+k*3+j);
            for(int k=0; k<n2sp; k ++) v2sp[i][k][j]= *(ptr_v2sp+i*n2sp*3+k*3+j);
		}
	}
	for(int i=0; i<n2nbr; i ++){ // v2nbr
		for(int j=0; j<3; j ++){
			v2nbr[i][j]= (*(ptr_v2nbr+i))[j];
            if(i>=n2nbr/2) if(v2nbr[i][j] != -v2nbr[i-n2nbr/2][j]) error(1, "(ltc_constructor) v2nbr isn't symmetry");
            v12nbr[i+n1nbr][j]= v2nbr[i][j];
            v123nbr[i+n1nbr][j]= v2nbr[i][j];
		}
	}
	for(int i=0; i<n3nbr; i ++){ // v3nbr
		for(int j=0; j<3; j ++){
			v3nbr[i][j]= (*(ptr_v3nbr+i))[j];
            if(i>=n3nbr/2) if(v3nbr[i][j] != -v3nbr[i-n3nbr/2][j]) error(1, "(ltc_constructor) v3nbr isn't symmetry");
            v123nbr[i+n1nbr+n2nbr][j]= v3nbr[i][j];
		}
	}
}

void class_initial::init_states_array(double compV, double compA, int nMlayer){
	// STATE 0: vacancy, 1: A atom, -1: B atom, 4: Vacuum

    double pV, pA;
    bool is_1vcc; // if compV>1, is 1 defect               //fh: this designates that you only will have 1 vacancy
    if(compV >1.0){ is_1vcc= true;  pV= 0;     pA= compA; }
    else          { is_1vcc= false; pV= compV; pA= compA*(1-compV); }

    for(int i=0; i<nx; i ++){	
	    for(int j=0; j<ny; j ++){	
	        for(int k=0; k<nz; k ++){
                srf[i][j][k]= false;  
				                               //fh: to get rid of surface, will have to get rid of this
                if(i<nMlayer || i>(nx-nMlayer-1))   states[i][j][k]= 4; // vacuum layers     //fh line is where vacuum layers initiated
		        else{
			            double ran= ran_generator();

                        if(ran < pV)           
                        {
                        	states[i][j][k]= 0;
                        	int ltcp= i*ny*nz+j*nz+k;
                        	ltcp_array[ltcp][0]=0;
                        	ltcp_array[ltcp][1]=0;
                        	ltcp_array[ltcp][2]=0;
                        	ltcp_array[ltcp][3]=0;
                        }    
                        else if(ran < (pV+pA))   
                        {
                        	states[i][j][k]= 1;
                        	int ltcp= i*ny*nz+j*nz+k;
                        	ltcp_array[ltcp][0]=1;
                        	ltcp_array[ltcp][1]=0;
                        	ltcp_array[ltcp][2]=0;
                        	ltcp_array[ltcp][3]=0;
                        }   
			            else         
			            {
			            	states[i][j][k]=-1;
			            	int ltcp= i*ny*nz+j*nz+k;
			            	ltcp_array[ltcp][0]=-1;
                        	ltcp_array[ltcp][1]=0;
                        	ltcp_array[ltcp][2]=0;
                        	ltcp_array[ltcp][3]=0;
			            }               
		        }
    }}}
    if(is_1vcc){
        states[0][0][0]= par_typeD;
        if(nMlayer != 0) error(1, "(init_states_array) nM != 0 & nV=1 is not allowed (vcc in vacuum)");
        cout << "compV gt 1: only 1 defect: " << par_typeD << endl;
    }

	nV= 0; nA= 0; nB= 0; nAA= 0; nBB= 0; nAB= 0; nM= 0;
	
	////////// CHECK //////////
	for(int i=0; i<nx; i ++){ 
	    for(int j=0; j<ny; j ++){ 
	        for(int k=0; k<nz; k ++){
		        switch(states[i][j][k]){
                    case  0:
			            list_vcc.push_back(vcc());
			            list_vcc[nV].x= i;
			            list_vcc[nV].y= j;
			            list_vcc[nV].z= k;
			            list_vcc[nV].ix= 0;
			            list_vcc[nV].iy= 0;
			            list_vcc[nV].iz= 0;
                        nV ++; break;
		            case  1:
                        nA ++; break;
                    case -1:
                        nB ++; break;
                    case  4:
                        for(int a=0; a<n1nbr; a ++){ // mark srf atoms
                            int x= pbc(i+v1nbr[a][0], nx);
                            int y= pbc(j+v1nbr[a][1], ny);
                            int z= pbc(k+v1nbr[a][2], nz);

                            if(     states[x][y][z] == 0) error(1, "(init_states_array) vcc adjacent to vacuum");
                            else if(states[x][y][z] != 4) srf[x][y][z]= true; //fh srf seems to be an array that marks surface atoms
                        } 
                        nM ++; break;
                    default: // itl
                        list_itl.push_back(itl());
	                    list_itl.back().x= i;
	                    list_itl.back().y= j;
	                    list_itl.back().z= k;
	                    list_itl.back().dir= (int) (ran_generator()*n1nbr);
	                    list_itl.back().head= 1; // choose the ltcp[0] because dir is randomly selected
	                    list_itl.back().ix= 0; 
	                    list_itl.back().iy= 0; 
	                    list_itl.back().iz= 0; 
	                    if(2==states[i][j][k])      nAA ++;
                        else if(3==states[i][j][k]) nAB ++;
                        else                        nBB ++;
                }
    }}}
	int nAtotal= nA + nB;
	//if(abs((double) nA/nAtotal-compA) > 0.01) error(1, "(init_states_array) the composition of generated conf is inconsistent of compA", 1, nA);
	////////// CHECK //////////
	
	cout << "The random solution configuration has been generated!" << endl;
	cout << "Vacancy: " << nV << endl;
    cout << "Vacuum:  " << nM << endl;
	cout << "Atype A: " << nA << ", pct: " << 100* (double) nA/(nAtotal) << "%" << endl;
	cout << "Atype B: " << nB << ", pct: " << 100* (double) nB/(nAtotal) << "%" << endl;
    
    cout << "number AB: " << nAB << endl;
} 


//fh new read restart function
void class_initial::read_restart(char name_restart[], char name_correction_factors[], long long int &ts_initial, double &time_initial, double &time_initial_corrected_A, double &time_initial_corrected_B)
{
	ifstream if_re(name_restart, ios::in);
	if (!if_re.is_open()) error(1, "(read_restart) the file is not opened!");

	long long int timestep;
	double time;
	double time_corrected_A;      //fh
	double time_corrected_B;      //fh

	int ntotal;
	if_re >> ntotal;

	char c_ltcp[5];
	if_re >> c_ltcp >> timestep >> time >> time_corrected_A >> time_corrected_B;   //fh addition at the end
	if (strcmp(c_ltcp, "ltcp") != 0) error(1, "(read_restart) please input a ltcp file (ltcp at the second line)"); // check
	ts_initial = timestep;
	time_initial = time;
	time_initial_corrected_A = time_corrected_A;   //fh addiition
	time_initial_corrected_B = time_corrected_B;   //fh addiition

	for (int i = 0; i<nx; i++) {
		for (int j = 0; j<ny; j++) {
			for (int k = 0; k<nz; k++) {
				states[i][j][k] = 1;
				srf[i][j][k] = false;
			}
		}
	}

	nV = 0; nA = 0; nB = 0; nAA = 0; nBB = 0; nAB = 0;
	for (int index = 0; index<ntotal; index++) {
		// caution: conf file X contains A atoms
		int type, i, j, k, is_srf, ix, iy, iz, dir, head;

		if_re >> type >> i >> j >> k;

		if (-1 == type || 1 == type) {
			if_re >> is_srf;
			srf[i][j][k] = is_srf;

			if (-1 == type) nB++;
		}
		else if (4 == type) nM++;
		else if (0 == type) {
			if_re >> ix >> iy >> iz;

			list_vcc.push_back(vcc());
			list_vcc[nV].x = i;
			list_vcc[nV].y = j;
			list_vcc[nV].z = k;
			list_vcc[nV].ix = ix;
			list_vcc[nV].iy = iy;
			list_vcc[nV].iz = iz;

			nV++;
		}
		else {
			if_re >> ix >> iy >> iz >> dir >> head;

			list_itl.push_back(itl());
			list_itl[nAA + nAB + nBB].x = i;
			list_itl[nAA + nAB + nBB].y = j;
			list_itl[nAA + nAB + nBB].z = k;
			list_itl[nAA + nAB + nBB].ix = ix;
			list_itl[nAA + nAB + nBB].iy = iy;
			list_itl[nAA + nAB + nBB].iz = iz;
			list_itl[nAA + nAB + nBB].dir = dir;
			list_itl[nAA + nAB + nBB].head = head;

			if (2 == type)         nAA++;
			else if (-2 == type)  nBB++;
			else if (3 == type)  nAB++;
			else error(1, "(read_restart) cant identify the type", 1, type);
		}

		states[i][j][k] = type;
	}
	nA = nx * ny*nz - nB - nV - nAA - nAB - nBB - nM;

	for (int i = 0; i<nx; i++) { // set up srf array
		for (int j = 0; j<ny; j++) {
			for (int k = 0; k<nz; k++) {
				if (states[i][j][k] != 4) continue;

				for (int a = 0; a<n1nbr; a++) { // mark srf atoms
					int x = pbc(i + v1nbr[a][0], nx);
					int y = pbc(j + v1nbr[a][1], ny);
					int z = pbc(k + v1nbr[a][2], nz);

					if (states[x][y][z] == 0) error(1, "(init_states_array) vcc adjacent to vacuum");
					else if (states[x][y][z] != 4) srf[x][y][z] = true;
				}
			}
		}
	}


	///////////////////////////////////////////////////////////////////////////////////////////////////
	//fh additions
	ifstream correction_factors(name_correction_factors, ios::in);
	if (!correction_factors.is_open()) error(1, "(read_restart) the file is not opened!");
	int input_time;
	int vacancy_x;
	int vacancy_y;
	int vacancy_z;
	int solute_x;
	int solute_y;
	int solute_z;
	int itl_x;
	int itl_y;
	int itl_z;
	correction_factors >>input_time>> vacancy_x >> vacancy_y >> vacancy_z >> solute_x >> solute_y >> solute_z >> itl_x >> itl_y >> itl_z;
	vacancy_correction_x = vacancy_x;
	vacancy_correction_y = vacancy_y;
	vacancy_correction_z = vacancy_z;
	solute_correction_x = solute_x;
	solute_correction_y = solute_y;
	solute_correction_z = solute_z;
	itl_correction_x = itl_x;
	itl_correction_y = itl_y;
	itl_correction_z = itl_z;


	cout << " vacancy correction: " << vacancy_correction_x << " " << vacancy_correction_y << " " << vacancy_correction_z;

	correction_factors.close();



	/////////////////////////////////////////////////////////////////////////////////////////////////////

	cout << "The configuration has been generated from the restart file!" << endl;
	cout << "Vacancy: " << nV << endl;
	cout << "Atype A: " << nA << ", pct: " << 100 * (double)nA / (nx*ny*nz) << "%" << endl;
	cout << "Atype B: " << nB << ", pct: " << 100 * (double)nB / (nx*ny*nz) << "%" << endl;
	cout << " Itl AA: " << nAA << endl;
	cout << " Itl AB: " << nAB << endl;
	cout << " Itl BB: " << nBB << endl;

	if_re.close();
}

void class_initial::init_par(){ // cal ABVI Ising model pars; see JPCM Huang et al.
	// 1st-nn class 1
	c1_44= ( (eAA1AA -8*eAA1A -8*eAA1B +12*eAA1AB +2*eAA1BB +12*eAB1BB -8*eA1BB -8*eB1BB +eBB1BB)
	       -12*eAB1AB + (-48*eA1V +48*eV1V -48*eV1B) + (16*eA1A +32*eA1B +16*eB1B) )/576; 
	 
	c1_42= ( (-eAA1AA +20*eAA1A +20*eAA1B -36*eAA1AB -2*eAA1BB -36*eAB1BB +20*eA1BB +20*eB1BB -eBB1BB)
	       +36*eAB1AB + (216*eA1V -216*eV1V +216*eV1B) + (-64*eA1A -128*eA1B -64*eB1B) )/576; 

	c1_22= ( (eAA1AA -32*eAA1A -32*eAA1B +60*eAA1AB +2*eAA1BB +60*eAB1BB -32*eA1BB -32*eB1BB +eBB1BB)
	       -60*eAB1AB + (-960*eA1V +960*eV1V -960*eV1B) + (256*eA1A +512*eA1B +256*eB1B) )/576; 
	
	// 1st-nn class 2
	c1_43= ( (eAA1AA -6*eAA1A -2*eAA1B +6*eAA1AB -6*eAB1BB +2*eA1BB +6*eB1BB -eBB1BB) +
	       (-12*eA1V +12*eV1B) + (8*eA1A -8*eB1B) )/288; 
	 
	c1_41= ( (-eAA1AA +12*eAA1A -4*eAA1B -6*eAA1AB +6*eAB1BB +4*eA1BB -12*eB1BB +eBB1BB) +
	       (48*eA1V -48*eV1B) + (-32*eA1A +32*eB1B) )/288; 
	
	c1_32= ( (-eAA1AA +18*eAA1A +14*eAA1B -30*eAA1AB +30*eAB1BB -14*eA1BB -18*eB1BB +eBB1BB) +
	       (60*eA1V -60*eV1B) + (-32*eA1A +32*eB1B) )/288; 
	
	c1_21= ( (eAA1AA -24*eAA1A -8*eAA1B +30*eAA1AB -30*eAB1BB +8*eA1BB +24*eB1BB -eBB1BB) +
	       (-240*eA1V +240*eV1B) + (128*eA1A -128*eB1B) )/288; 
	
	// 1st-nn class 3
	c1_33= ( (eAA1AA -4*eAA1A +4*eAA1B -2*eAA1BB +4*eA1BB -4*eB1BB +eBB1BB) +
	       (4*eA1A -8*eA1B +4*eB1B) )/144; 
	
	c1_31= ( (-eAA1AA +10*eAA1A -10*eAA1B +2*eAA1BB -10*eA1BB +10*eB1BB -eBB1BB) +
	       (-16*eA1A +32*eA1B -16*eB1B) )/144; 
	
	c1_11= ( (eAA1AA -16*eAA1A +16*eAA1B -2*eAA1BB +16*eA1BB -16*eB1BB +eBB1BB) +
	       (64*eA1A -128*eA1B +64*eB1B) )/144; 
	
	// 2nd-nn class 1
	c2_44= ( (eAA2AA -8*eAA2A -8*eAA2B +12*eAA2AB +2*eAA2BB +12*eAB2BB -8*eA2BB -8*eB2BB +eBB2BB)
	       -12*eAB2AB + (-48*eA2V +48*eV2V -48*eV2B) + (16*eA2A +32*eA2B +16*eB2B) )/576; 
	 
	c2_42= ( (-eAA2AA +20*eAA2A +20*eAA2B -36*eAA2AB -2*eAA2BB -36*eAB2BB +20*eA2BB +20*eB2BB -eBB2BB)
	       +36*eAB2AB + (216*eA2V -216*eV2V +216*eV2B) + (-64*eA2A -128*eA2B -64*eB2B) )/576; 

	c2_22= ( (eAA2AA -32*eAA2A -32*eAA2B +60*eAA2AB +2*eAA2BB +60*eAB2BB -32*eA2BB -32*eB2BB +eBB2BB)
	       -60*eAB2AB + (-960*eA2V +960*eV2V -960*eV2B) + (256*eA2A +512*eA2B +256*eB2B) )/576; 
	
	// 2nd-nn class 2
	c2_43= ( (eAA2AA -6*eAA2A -2*eAA2B +6*eAA2AB -6*eAB2BB +2*eA2BB +6*eB2BB -eBB2BB) +
	       (-12*eA2V +12*eV2B) + (8*eA2A -8*eB2B) )/288; 
	 
	c2_41= ( (-eAA2AA +12*eAA2A -4*eAA2B -6*eAA2AB +6*eAB2BB +4*eA2BB -12*eB2BB +eBB2BB) +
	       (48*eA2V -48*eV2B) + (-32*eA2A +32*eB2B) )/288; 
	
	c2_32= ( (-eAA2AA +18*eAA2A +14*eAA2B -30*eAA2AB +30*eAB2BB -14*eA2BB -18*eB2BB +eBB2BB) +
	       (60*eA2V -60*eV2B) + (-32*eA2A +32*eB2B) )/288; 
	
	c2_21= ( (eAA2AA -24*eAA2A -8*eAA2B +30*eAA2AB -30*eAB2BB +8*eA2BB +24*eB2BB -eBB2BB) +
	       (-240*eA2V +240*eV2B) + (128*eA2A -128*eB2B) )/288; 
	
	// 2nd-nn class 3
	c2_33= ( (eAA2AA -4*eAA2A +4*eAA2B -2*eAA2BB +4*eA2BB -4*eB2BB +eBB2BB) +
	       (4*eA2A -8*eA2B +4*eB2B) )/144; 
	
	c2_31= ( (-eAA2AA +10*eAA2A -10*eAA2B +2*eAA2BB -10*eA2BB +10*eB2BB -eBB2BB) +
	       (-16*eA2A +32*eA2B -16*eB2B) )/144; 
	
	c2_11= ( (eAA2AA -16*eAA2A +16*eAA2B -2*eAA2BB +16*eA2BB -16*eB2BB +eBB2BB) +
	       (64*eA2A -128*eA2B +64*eB2B) )/144; 

	// 1st-nn constants (for itl jump only)
	c1_40= ( ( eAA1AB +eAB1BB) + (-4*eA1AB  +6*eAB1AB  -4*eAB1B) + (-4*eA1V  +6*eV1V  -4*eV1B) )/24; 
	c1_30= ( ( eAA1AB -eAB1BB) + (-2*eA1AB             +2*eAB1B) + (-2*eA1V           +2*eV1B) )/12; 
	c1_20= ( (-eAA1AB -eAB1BB) + (16*eA1AB -30*eAB1AB +16*eAB1B) + (16*eA1V -30*eV1V +16*eV1B) )/24; 
	c1_10= ( (-eAA1AB +eAB1BB) + ( 8*eA1AB             -8*eAB1B) + ( 8*eA1V           -8*eV1B) )/12; 
	
	// 2nd-nn constants (for itl jump only)
	c2_40= ( ( eAA2AB +eAB2BB) + (-4*eA2AB  +6*eAB2AB  -4*eAB2B) + (-4*eA2V  +6*eV2V  -4*eV2B) )/24; 
	c2_30= ( ( eAA2AB -eAB2BB) + (-2*eA2AB             +2*eAB2B) + (-2*eA2V           +2*eV2B) )/12; 
	c2_20= ( (-eAA2AB -eAB2BB) + (16*eA2AB -30*eAB2AB +16*eAB2B) + (16*eA2V -30*eV2V +16*eV2B) )/24; 
	c2_10= ( (-eAA2AB +eAB2BB) + ( 8*eA2AB             -8*eAB2B) + ( 8*eA2V           -8*eV2B) )/12; 

	// 1st-nn C00 (for itl jump only) (C00 are calculated on the fly(change with time))
	c1_00   = eAB1AB +eV1V;
	c1_0_ABA= eA1AB -eA1V -0.5*eAB1AB + 0.5*eV1V;
	c1_0_ABB= eAB1B -eV1B -0.5*eAB1AB + 0.5*eV1V;
	c1_0_A  = n1nbr*(-eA1AB +0.5*eAB1AB);
	c1_0_B  = n1nbr*(-eAB1B +0.5*eAB1AB);
	c1_0_V  = n1nbr*(-0.5*eAB1AB);
	c1_0_AA = n1nbr*( 0.5*eV1V);
	c1_0_AB = n1nbr*(-0.5*eV1V);
	c1_0_BB = n1nbr*( 0.5*eV1V);

	// 2nd-nn C00 (for itl jump only) (C00 are calculated on the fly(change with time))
	c2_00   = eAB2AB +eV2V;
	c2_0_ABA= eA2AB -eA2V -0.5*eAB2AB + 0.5*eV2V;
	c2_0_ABB= eAB2B -eV2B -0.5*eAB2AB + 0.5*eV2V;
	c2_0_A  = n2nbr*(-eA2AB +0.5*eAB2AB);
	c2_0_B  = n2nbr*(-eAB2B +0.5*eAB2AB);
	c2_0_V  = n2nbr*(-0.5*eAB2AB);
	c2_0_AA = n2nbr*( 0.5*eV2V);
	c2_0_AB = n2nbr*(-0.5*eV2V);
	c2_0_BB = n2nbr*( 0.5*eV2V);

	// 1st-nn Vacuum
    c1_MM   = eM1M -eV1V;
    c1_MA   = eM1A -eA1V;
    c1_MV   = eM1V -eV1V;
    c1_MB   = eM1B -eV1B;
	
    // 1st-nn Vacuum
    c2_MM   = eM2M -eV2V;
    c2_MA   = eM2A -eA2V;
    c2_MV   = eM2V -eV2V;
    c2_MB   = eM2B -eV2B;
    
    if(0==c2_44 && 0==c2_43 && 0==c2_42 && 0==c2_41 && 0==c2_33 && 0==c2_32 && 0==c2_31 && 0==c2_22 && 0==c2_21 && 0==c2_11) is_e2nbr= false;
	else is_e2nbr= true;

	// print out the parameters to log file
	cout << "\n##### Energy calculation parameters #####" << endl; 
	
	cout << "temperature= " << temp << ", beta= " << beta << endl;
	cout << "correlation factor= " << corrfac<< endl;
	printf("Vacancy mu= %f %f\n", muvA, muvB);
	printf("Interstitial mu= %f %f\n", muiA, muiB);
	printf("Vacancy Em= %f %f\n", emvA, emvB);
	printf("Interstitial Em= %f %f\n", emiA, emiB);
	printf("Rotation Er(AA, AB, BB)= %f %f %f\n", erAA, erAB, erBB);
	printf("Interstitial Em(AA AB BB)= %f %f %f (Dubey, CMS 2015)\n", emiAA, emiAB, emiBB);
	
	cout << "\n##### Input epsilons: #####" << endl;
	cout << "(1st neigbor)" << endl;
	printf("AA-AA: %f, AA-A: %f, AA-AB: %f, AA-B: %f, AA-BB: %f\n", eAA1AA, eAA1A, eAA1AB, eAA1B, eAA1BB);
	printf("A-A:   %f, A-V:  %f, A-AB: %f, A-B:   %f, A-BB: %f\n", eA1A, eA1V, eA1AB, eA1B, eA1BB);
	printf("V-V:   %f, V-B:  %f\n", eV1V, eV1B);
	printf("AB-AB: %f, AB-B: %f, AB-BB: %f\n", eAB1AB, eAB1B, eAB1BB);
	printf("B-B:   %f, B-BB: %f\n", eB1B, eB1BB);
	printf("BB-BB: %f\n", eBB1BB);
	cout << "(2nd neigbor)" << endl;
	printf("AA-AA: %f, AA-A: %f, AA-AB: %f, AA-B: %f, AA-BB: %f\n", eAA2AA, eAA2A, eAA2AB, eAA2B, eAA2BB);
	printf("A-A:   %f, A-V:  %f, A-AB: %f, A-B:   %f, A-BB: %f\n", eA2A, eA2V, eA2AB, eA2B, eA2BB);
	printf("V-V:   %f, V-B:  %f\n", eV2V, eV2B);
	printf("AB-AB: %f, AB-B: %f, AB-BB: %f\n", eAB2AB, eAB2B, eAB2BB);
	printf("B-B:   %f, B-BB: %f\n", eB2B, eB2BB);
	printf("BB-BB: %f\n", eBB2BB);
	
	cout << "\n##### Ising formulation constants: #####" << endl;
	cout << "(1st neighbor)" << endl;
	printf("Class 1\nC44: %f, C42: %f, C22: %f\n", c1_44, c1_42, c1_22);
	printf("Class 2\nC43: %f, C41: %f, C32: %f, C21: %f\n", c1_43, c1_41, c1_32, c1_21);
	printf("Class 3\nC33: %f, C31: %f, C11: %f\n", c1_33, c1_31, c1_11);
	printf("Class 0\nC40: %f, C30: %f, C20: %f, C10: %f\n", c1_40, c1_30, c1_20, c1_10);
	printf("C00: %f, C0_ABA: %f, C0_ABB: %f, C0_AA: %f, C0_A: %f, C0_AB: %f, C0_B: %f, C0_BB: %f\n",c1_00, c1_0_ABA, c1_0_ABB, c1_0_AA, c1_0_A, c1_0_AB, c1_0_B, c1_0_BB);
	printf("Vacuum \nCMM: %f, CMA: %f, CMV: %f, CMB: %f\n", c1_MM, c1_MA, c1_MV, c1_MB);
	cout << "(2nd neighbor)" << endl;
	printf("Class 1\nC44: %f, C42: %f, C22: %f\n", c2_44, c2_42, c2_22);
	printf("Class 2\nC43: %f, C41: %f, C32: %f, C21: %f\n", c2_43, c2_41, c2_32, c2_21);
	printf("Class 3\nC33: %f, C31: %f, C11: %f\n", c2_33, c2_31, c2_11);
	printf("Class 0\nC40: %f, C30: %f, C20: %f, C10: %f\n", c2_40, c2_30, c2_20, c2_10);
	printf("C00: %f, C0_ABA: %f, C0_ABB: %f, C0_AA: %f, C0_A: %f, C0_AB: %f, C0_B: %f, C0_BB: %f\n",c2_00, c2_0_ABA, c2_0_ABB, c2_0_AA, c2_0_A, c2_0_AB, c2_0_B, c2_0_BB);
	printf("Vacuum \nCMM: %f, CMA: %f, CMV: %f, CMB: %f\n", c2_MM, c2_MA, c2_MV, c2_MB);

	if(is_e2nbr) 	cout << "\n2nd nn parameters are non-zero" << endl;
	else            cout << "\n2nd nn are 0, skip 2nd-nn calculations" << endl;
}

void class_initial::try_interpenetrating()
{

}

