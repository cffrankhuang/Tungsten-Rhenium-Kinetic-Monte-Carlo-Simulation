#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <vector>
#include "kmc_global.h"

using namespace std;

////////// GLOBAL VARIABLES //////////
long long int timestep;
double totaltime;
double totaltime_vacancy_adjust_A;     //fh addition
double totaltime_vacancy_adjust_B;     //fh addition

#define MAX_NNBR 20
double vbra[3][3];
int n1nbr, n2nbr, n3nbr;
int v1nbr[MAX_NNBR][3];
int v2nbr[MAX_NNBR][3];
int v3nbr[MAX_NNBR][3];
int n1sp, n2sp;
int v1sp[MAX_NNBR][MAX_NNBR][3];
int v2sp[MAX_NNBR][MAX_NNBR][3];
int n12nbr, n123nbr;
int v12nbr[MAX_NNBR*2][3];
int v123nbr[MAX_NNBR*3][3];

int nA, nB, nV, nAA, nBB, nAB, nM;
int sum_mag;
vector<vector<vector<int>>> states(nx, vector<vector<int>>(ny, vector<int>(nz)));
vector<vector<vector<bool>>> srf(nx, vector<vector<bool>>(ny, vector<bool>(nz)));

int frenkel_pairs_generated;     //fh addition
int vacancy_correction_x=0;   //fh addition
int vacancy_correction_y=0;   //fh addition
int vacancy_correction_z=0;   //fh addition
int solute_correction_x = 0;
int solute_correction_y = 0;
int solute_correction_z = 0;
int itl_correction_x = 0;
int itl_correction_y = 0;
int itl_correction_z = 0;
double vacancy_formation_energy = 0.0;
int ltcp_array[nx*ny*nz][4];

FILE * his_sol;
FILE * his_def;

FILE* his_dis;               // fh edit
FILE* his_dis_correction; //fh addition
FILE* his_dis_itl;           //fh addition
FILE* correction_factors_file;

FILE * his_srf;
FILE * out_engy;
FILE * out_vdep;
FILE * out_sro;
FILE * out_msd;

vector <vcc> list_vcc;
vector <itl> list_itl;
vector <int> list_sink;

int N_genr= 0;
int njump[10]= {0};
long long int Vja[2]= {0};
long long int Ija[2]= {0};
double acc_dsroV=  0;
double acc_dsroG=  0;
double acc_dsroRi= 0;
double acc_dsroRv= 0;

double h0;
double c1_44, c1_43, c1_42, c1_41, c1_33, c1_32, c1_31, c1_22, c1_21, c1_11, c1_40, c1_30, c1_20, c1_10, c1_00;
double c1_0_ABA, c1_0_ABB, c1_0_A, c1_0_B, c1_0_V, c1_0_AA, c1_0_AB, c1_0_BB;
double c2_44, c2_43, c2_42, c2_41, c2_33, c2_32, c2_31, c2_22, c2_21, c2_11, c2_40, c2_30, c2_20, c2_10, c2_00;
double c2_0_ABA, c2_0_ABB, c2_0_A, c2_0_B, c2_0_V, c2_0_AA, c2_0_AB, c2_0_BB;
double c1_MM, c1_MA, c1_MV, c1_MB;
double c2_MM, c2_MA, c2_MV, c2_MB;
double unc1_44, unc1_43, unc1_42, unc1_41, unc1_33, unc1_32, unc1_31, unc1_22, unc1_21, unc1_11, unc1_40, unc1_30, unc1_20, unc1_10, unc1_00;
double unc2_44, unc2_43, unc2_42, unc2_41, unc2_33, unc2_32, unc2_31, unc2_22, unc2_21, unc2_11, unc2_40, unc2_30, unc2_20, unc2_10, unc2_00;
bool is_e2nbr;

////////// GLOBAL FUNCTIONS //////////
void error(int nexit, string errinfo, int nnum, double num1, double num2){
	// exit number represents:
	// 0: main; 1: class_system; 2: class_events
	
	cout << "\nError: ";
	switch(nexit){
		case 0:  cout << "In main function "; break;
		case 1:  cout << "In class_initial ";  break;
		case 2:  cout << "In class_events ";  break;
		default: cout << "UNKNOWN NEXIT" << endl;
	}

	cout << errinfo;
	switch(nnum){
		case 0:  cout << endl; break;
		case 1:  cout << ": " << num1 << endl; break;
		case 2:  cout << ": " << num1 << " " << num2 << endl; break;
		default: cout << "!!!ERROR FUNCTION MALFUNCTION!!! WRONG NNUM!!!" << endl;
	}
	cout << endl;
	exit(1); 
}

void error(int nexit, string errinfo, char c[]){
	// exit number represents:
	// 0: main; 1: class_system; 2: class_events
	
	cout << "\nError: ";
	switch(nexit){
		case 0:  cout << "In main function "; break;
		case 1:  cout << "In class_system ";  break;
		case 2:  cout << "In class_events ";  break;
		default: cout << "UNKNOWN NEXIT" << endl;
	}

	cout << errinfo << " " << c << endl;
	exit(1); 
}

//fh: edit
void error(int nexit, string errinfo, string c, int number) {
	// exit number represents:
	// 0: main; 1: class_system; 2: class_events

	cout << "\nError: ";
	switch (nexit) {
	case 0:  cout << "In main function "; break;
	case 1:  cout << "In class_system ";  break;
	case 2:  cout << "In class_events ";  break;
	default: cout << "UNKNOWN NEXIT" << endl;
	}

	cout << errinfo << " " << c << endl;
	exit(1);
}

double ran_generator(){
	static bool first= true;
	if(first){
		srand(time(NULL));
		first= false;
	}
	
	return rand()/((double) RAND_MAX+1.0);
}

int pbc(int x_, int nx_){ // Periodic Boundary Condition
	if	(x_<-nx_ || x_>=2*nx_)
		error(1, "(pbc) input x is out of bound", 2, x_, nx_);

	if	(x_<0)          return (x_ + nx_);
	else if	(x_<nx_)    return  x_;
	else                return (x_ - nx_);
}

int pbc_change(int x_, int nx_, int axis, int &crossed) { // Periodic Boundary Condition
	if (x_ < -nx_ || x_ >= 2 * nx_)
	{
	
	error(1, "(pbc) input x is out of bound", 2, x_, nx_);
	}

	if (x_ < 0)
	{
		
		switch (axis) 
		{
		case 0:
			vacancy_correction_x -= nx_;
			crossed += nx_;
			break;
		case 1:
			vacancy_correction_y -= nx_;
			crossed += nx_;
			break;
		case 2:
			vacancy_correction_z -= nx_;
			crossed += nx_;
			break;

		}
		return (x_ + nx_);
	}
	else if (x_<nx_)    return  x_;
	else
	{
		
		switch (axis)
		{
		case 0:
			vacancy_correction_x += nx_;
			crossed -= nx_;
			break;
		case 1:
			vacancy_correction_y += nx_;
			crossed -= nx_;
			break;
		case 2:
			vacancy_correction_z += nx_;
			crossed -= nx_;
			break;

		}
		return (x_ - nx_);
	}
}

int pbc_change_itl(int x_, int nx_, int axis, int &crossed) { // Periodic Boundary Condition
	if (x_ < -nx_ || x_ >= 2 * nx_)
	{

		error(1, "(pbc) input x is out of bound", 2, x_, nx_);
	}

	if (x_ < 0)
	{
		
		switch (axis)
		{
		case 0:
			itl_correction_x -= nx_;
			crossed += nx_;
			break;
		case 1:
			itl_correction_y -= nx_;
			crossed += nx_;
			break;
		case 2:
			itl_correction_z -= nx_;
			crossed += nx_;
			break;

		}
		return (x_ + nx_);
	}
	else if (x_<nx_)    return  x_;
	else
	{
		
		switch (axis)
		{
		case 0:
			itl_correction_x += nx_;
			crossed -= nx_;
			break;
		case 1:
			itl_correction_y += nx_;
			crossed -= nx_;
			break;
		case 2:
			itl_correction_z += nx_;
			crossed -= nx_;
			break;

		}
		return (x_ - nx_);
	}
}

void write_conf(int flag){
// A atoms are omitted
// flag: 0: t0; 1: timestep; 2: time; 3: RESTART
	ofstream of_xyz;
	ofstream of_ltcp;

	// determine the names of conf files
	char name_xyz[40], name_ltcp[40];
	
    if(0==flag){      // t0
		of_xyz.open("t0.xyz");
		of_ltcp.open("t0.ltcp");
	}
	else if(1==flag){ // step

		
		sprintf(name_xyz, "%lld.xyz", timestep);    of_xyz.open(name_xyz);
		sprintf(name_ltcp, "%lld.ltcp", timestep);  of_ltcp.open(name_ltcp);
		
	}
    else if(2==flag){ // time

		sprintf(name_xyz, "time%.2f.xyz", totaltime);   of_xyz.open(name_xyz);
		sprintf(name_ltcp, "time%.2f.ltcp", totaltime); of_ltcp.open(name_ltcp);
		
    }
    else if(3==flag) of_ltcp.open("RESTART");
    else error(0, "wrong flag: ", 1, flag);

	if(flag != 3 && !of_xyz.is_open()) error(1, "(write_conf) xyz file is not opened!");// check
	if(!of_ltcp.is_open()) error(1, "(write_conf) ltcp file is not opened!");           // check
	
	// write out data
	of_ltcp << nx*ny*nz-nA << "\n" << "ltcp " << timestep << " ";
	of_ltcp << setprecision(15) << totaltime << " " << totaltime_vacancy_adjust_A << " " << totaltime_vacancy_adjust_B << "\n";        //fh adjusted
    if(flag != 3){
	    of_xyz << nx*ny*nz-nA << "\n" << "xyz " << timestep << " ";
	    of_xyz << setprecision(15) << totaltime << " " << totaltime_vacancy_adjust_A << " " << totaltime_vacancy_adjust_B << "\n";       //fh adjusted
    }

	for(int i=0; i<nx; i ++){
		for(int j=0; j<ny; j ++){
			for(int k=0; k<nz; k ++){
				double x= i*vbra[0][0] + j*vbra[1][0] + k*vbra[2][0];
				double y= i*vbra[0][1] + j*vbra[1][1] + k*vbra[2][1];
				double z= i*vbra[0][2] + j*vbra[1][2] + k*vbra[2][2];
		
				if(1==states[i][j][k]) continue;
                else if(-1==states[i][j][k]){   // B atom
					if(flag != 3) of_xyz << states[i][j][k] << " " << x << " " << y << " " << z << endl;
					of_ltcp << states[i][j][k] << " " << i << " " << j << " " << k << " ";
				    if(srf[i][j][k]) of_ltcp << "1" << endl;
                    else             of_ltcp << "0" << endl;
                }
				else if(0==states[i][j][k]){    // vcc
					int id; for(id=0; list_vcc[id].x != i && list_vcc[id].y != j && list_vcc[id].z != k; id ++);
					if(flag != 3) of_xyz << states[i][j][k] << " " << x << " " << y << " " << z << " " << endl;
					of_ltcp << states[i][j][k] << " " << i << " " << j << " " << k << " " 
					        << list_vcc[id].ix << " " << list_vcc[id].iy << " " << list_vcc[id].iz << endl;
				}
                else if(4==states[i][j][k]){    // srf
					if(flag != 3) of_xyz << states[i][j][k] << " " << x << " " << y << " " << z << endl;
					of_ltcp << states[i][j][k] << " " << i << " " << j << " " << k << endl;
				}
				else{                           // itl
					int id; for(id=0; list_itl[id].x != i && list_itl[id].y != j && list_itl[id].z != k; id ++);
					int type= states[i][j][k];
					if(flag != 3) of_xyz  << type << " " << x << " " << y << " " << z << " " << endl; 
					of_ltcp << type << " " << i << " " << j << " " << k << " "
					        << list_itl[id].ix << " " << list_itl[id].iy << " " << list_itl[id].iz << " "
					        << list_itl[id].dir << " " << list_itl[id].head << endl;
				}
	}}}
	
	if(flag != 3) of_xyz.close();
	of_ltcp.close();
}

void write_hissol(){
	int ncheck= 0;
    vector <vector<int>> list_srf; // A list store srf info

    // OUTPUT his_sol
	fprintf(his_sol, "%d\n", nB);
	fprintf(his_sol, "T: %lld %e\n", timestep, totaltime);
	for(int i=0; i<nx; i++){
	    for(int j=0; j<ny; j++){
	        for(int k=0; k<nz; k++){
		        if(-1== states[i][j][k]){
			        ncheck ++;
                    int ltcp= i*ny*nz+j*nz+k;
			        fprintf(his_sol, "%d\n", ltcp);
                }

		        if(srf[i][j][k]) list_srf.push_back({i, j, k}); // know # before write
    }}}
	if(ncheck != nB) error(0, "(write_hissol) nB inconsistent", 2, ncheck, nB);
    
    // OUTPUT his_srf
    if(list_srf.size()==0) return;
	
    fprintf(his_srf, "%lu\n", list_srf.size());
	fprintf(his_srf, "T: %lld %e\n", timestep, totaltime);
    for(int i=0; i<list_srf.size(); i++){
        int x= list_srf[i][0];
        int y= list_srf[i][1];
        int z= list_srf[i][2];
		fprintf(his_srf, "%d %d %d %d\n", states[x][y][z], x, y, z);
	}

    fflush(his_sol);
    fflush(his_srf);
}

void write_hisdef(){
	// OUTPUT his_def
    fprintf(his_def, "%lu\n", list_vcc.size()+list_itl.size());         //fh print out number of vacancies + number interstitials
	fprintf(his_def, "T: %lld %e\n", timestep, totaltime);             //fh print out timestep and total time
    
    for(int i=0; i<list_vcc.size(); i++)                            //fh print out custom index
	{
        int ltcp= list_vcc[i].x*ny*nz + list_vcc[i].y*nz + list_vcc[i].z;
		fprintf(his_def, "0 %d %d %d %d\n", ltcp, list_vcc[i].ix, list_vcc[i].iy, list_vcc[i].iz);   //fh print out type, custom index, coordinates
	}
	for(int i=0; i<list_itl.size(); i++)   //fh loop through list of interstitials
	{
        int ltcp= list_itl[i].x*ny*nz + list_itl[i].y*nz + list_itl[i].z;       //fh calculate custom index 
		int type= states[list_itl[i].x][list_itl[i].y][list_itl[i].z];          //fh retrieve type
		fprintf(his_def, "%d %d %d %d %d\n", type, ltcp, list_itl[i].ix, list_itl[i].iy, list_itl[i].iz);    //fh print out type, custom index, and coordinates
	}
    
    fflush(his_def);
}

void write_hisdis()                          //fh function to acquire r2
{
	// OUTPUT his_def+
	fprintf(his_dis, "%lu\n", list_vcc.size() + nB);  //fh output number of vacancies and B atoms
	fprintf(his_dis, "T: %lld %e %e %e\n", timestep, totaltime, totaltime_vacancy_adjust_A, totaltime_vacancy_adjust_B);         //fh changed it
	

	for (int i = 0; i < list_vcc.size(); i++)       //fh loop through vacancies
	{
		int ltcp = list_vcc[i].x * ny * nz + list_vcc[i].y * nz + list_vcc[i].z;
		fprintf(his_dis, "0 %d %d %d %d\n", ltcp, list_vcc[i].x+ltcp_array[ltcp][1], list_vcc[i].y+ltcp_array[ltcp][2], list_vcc[i].z+ltcp_array[ltcp][3]);
	
	}
	//fh: need to print out all the B atoms
	//int first = 0;
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) 
			{

				if (states[i][j][k] == -1)
				{
					int ltcp = i * ny * nz + j * nz + k;
					int type = -1;

					fprintf(his_dis, "%d %d %d %d %d\n", type, ltcp, i+ltcp_array[ltcp][1], j+ltcp_array[ltcp][2], k+ltcp_array[ltcp][3]);
					
				}
				
			}
		}
	}
	fflush(his_dis);
}

void write_hisdis_correction()                          //fh function to acquire r2
{
	// OUTPUT his_def+
	fprintf(his_dis_correction, "%lu\n", list_vcc.size() + nB);  // output number of vacancies and B atoms
	

	for (int i = 0; i < list_vcc.size(); i++)       //loop through vacancies
	{
		int ltcp = list_vcc[i].x * ny * nz + list_vcc[i].y * nz + list_vcc[i].z;
		
		fprintf(his_dis_correction, "0 %d %d %d %d\n", ltcp, list_vcc[i].ix + vacancy_correction_x, list_vcc[i].iy + vacancy_correction_y, list_vcc[i].iz + vacancy_correction_z);
	}
	// need to print out all the B atoms
	int first = 0;
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++)
			{

				if (states[i][j][k] == -1)
				{
					int ltcp = i * ny * nz + j * nz + k;
					int type = -1;
					if (first == 0)
					{
						fprintf(his_dis_correction, "%d %d %d %d %d\n", type, ltcp, i + solute_correction_x, j + solute_correction_y, k + solute_correction_z);
					}
					else
					{
						fprintf(his_dis_correction, "%d %d %d %d %d\n", type, ltcp, i, j, k);
					}

				}
				first++;
			}
		}
	}
	fflush(his_dis_correction);
}

//fh addition
void write_hisdis_itl()
{
	// OUTPUT his_def+
	fprintf(his_dis_itl, "%lu\n", list_itl.size() + nB);  //output number of vacancies and B atoms
	fprintf(his_dis_itl, "T: %lld %e\n", timestep, totaltime);

	for (int i = 0; i < list_itl.size(); i++)       //loop through vacancies
	{
		int ltcp = list_itl[i].x * ny * nz + list_itl[i].y * nz + list_itl[i].z;
		
		int type =states[list_itl[i].x][list_itl[i].y][list_itl[i].z];
		fprintf(his_dis_itl, "%d %d %d %d %d\n", type, ltcp, list_itl[i].ix + itl_correction_x, list_itl[i].iy + itl_correction_y, list_itl[i].iz + itl_correction_z);
	}
	//need to print out all the B atoms
	int first = 0;
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++)
			{

				if (states[i][j][k] == -1)
				{
					int ltcp = i * ny * nz + j * nz + k;
					int type = -1;
					if (first == 0)
					{
						fprintf(his_dis_itl, "%d %d %d %d %d\n", type, ltcp, i + solute_correction_x, j + solute_correction_y, k + solute_correction_z);
					}
					else
					{
						fprintf(his_dis_itl, "%d %d %d %d %d\n", type, ltcp, i, j, k);
					}

				}
				first++;
			}
		}
	}
	fflush(his_dis_itl);
}

int cal_Bnbr(int N_Bnbr, int x, int y, int z){
    if(0==N_Bnbr){
        int n= 0;
        if(-1==states[x][y][z]) n ++;
        for(int a=0; a<n1nbr; a ++){ // search 1st-nn for B
		    int x1= pbc(x+v1nbr[a][0], nx);
		    int y1= pbc(y+v1nbr[a][1], ny);
		    int z1= pbc(z+v1nbr[a][2], nz);
		    if(-1==states[x1][y1][z1]) n ++;
	    }
	    for(int b=0; b<n2nbr; b ++){ // search 2nd-nn for B
		    int x2= pbc(x+v2nbr[b][0], nx);
		    int y2= pbc(y+v2nbr[b][1], ny);
		    int z2= pbc(z+v2nbr[b][2], nz);
		    if(-1==states[x2][y2][z2]) n ++;
        }

        return n;
    }
    else return N_Bnbr;
}

double cal_sro(){
    double cA= nA*1.0/(nx*ny*nz);
    double sro= 0;

    int ncheck= 0;
	for(int i=0; i<nx; i ++){
		for(int j=0; j<ny; j ++){
			for(int k=0; k<nz; k ++){
				int state0= states[i][j][k];

                if(-1==state0 || 3==state0){ // AB are considered B too
                    ncheck ++;
                    int nAnbr= 0;

				    for(int a=0; a<n1nbr; a ++){ // 1st neighbors
					    int x= pbc(i+(*(v1nbr+a))[0], nx);
					    int y= pbc(j+(*(v1nbr+a))[1], ny);
					    int z= pbc(k+(*(v1nbr+a))[2], nz);
                        int state1= states[x][y][z];

                        if(1==state1) nAnbr ++;
                    }
				    
                    for(int a=0; a<n2nbr; a ++){ // 1st neighbors
					    int x= pbc(i+(*(v2nbr+a))[0], nx);
					    int y= pbc(j+(*(v2nbr+a))[1], ny);
					    int z= pbc(k+(*(v2nbr+a))[2], nz);
                        int state2= states[x][y][z];

                        if(1==state2) nAnbr ++;
                    }

                    sro += 1.0 - nAnbr*1.0/(n1nbr+n2nbr)/cA;
                }
    }}}
    
    
    if(ncheck != (nAB+nB)) error(2, "(cal_sro) number inconsistent", 2, ncheck, nAB+nB);
    return sro/(nAB+nB);
}

