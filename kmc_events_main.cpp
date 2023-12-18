#include <cstdio>
#include <iostream>
#include <vector>
#include "kmc_global.h"
#include "kmc_events.h"

using namespace std;

double class_events::main(double &time_correction_factor_A, double &time_correction_factor_B)
{
	// a probability map will first generated by calculating all possible moves
	// then randomly picked the ACTUAL move based on the probability map

	// defect information
	vector <int>    etype; // type of the event: 0: ITL JUMP; 1: VCC JUMP; 7: F-P GENR; 8: VCC CRTN
	vector <double> rates; // transition rates
	vector <int>    ilist; // IDs in the lists
	vector <int>	inbr;  // the index of nbr
	
	// perform imaginary jumps and cal rates
	double irates= cal_ratesI(etype, rates, ilist, inbr); // WARNING: irates before vrates so the recb map can be generated
	double vrates= cal_ratesVsp(etype, rates, ilist, inbr);
    double crates= cvcc_rates;
	
	double sum_rates= vrates + irates + crates; // sum of all rates
    if(is_genr){ // the genr event
        etype.push_back(7); 
        rates.push_back(rate_genr);
        sum_rates += rate_genr;
    }

	
    if(abs(sum_rates)<1e-10) error(2, "(main) rate= 0, end up simulation");

	


    // check
    if(nA+nB+nV+nAA+nBB+nAB+nM != nx*ny*nz) error(2, "(jump) numbers of ltc points arent consistent, diff=", 1, nA+nB+nV+nAA+nBB+nAB+nM-nx*ny*nz);


    // perform the actual jump
	double ran= ran_generator();     //fh ran is a number between 0-1
        
    if(ran < crates/sum_rates)   //fh crates is the rate of forming a vacancy
	{ // vcc creation from srf
        double acc_cr= 0;          
        for(auto it= cvcc.begin(); it != cvcc.end(); it ++)  //fh loop through the surface atoms list
		{
            for(int a=0; a< (it->second).rates.size(); a ++)  //fh loop through the rates each surface atom possesses
			{ // loop tho creation paths of the atom
                double rate_a= (it->second).rates[a];
                if( (ran >= acc_cr) && (ran < (acc_cr + rate_a/sum_rates) ) )  //fh if random number is larger than acc_cr
					                                                                 //and random number less than accumulated rate and sum_rates
				{
                    create_vcc(it->first, (it->second).mltcp[a]);           //fh create vacancy
                    goto actionDONE;
                }
                        
                acc_cr += rate_a/sum_rates;
            }
        }
    }
    else{ // jump & genr
		//fh this section is for the case that surface reaction is not chosen
	    double acc_rate= crates/sum_rates; // accumulated rate
	    for(int i=0; i<rates.size(); i ++){
            if( (ran >= acc_rate) && (ran < (acc_rate + rates[i]/sum_rates) ) ){
                double sro0;
                switch(etype[i]){
                    case 0:  actual_jumpI(ilist[i], inbr[i]); recb_checki(ilist[i]); break;
					case 1:  
					{
						actual_jumpV(ilist[i], inbr[i]);
						//fh addition//////////////////////////////////////////////
						double vacancy_concentration = (1.0) / (nx*ny*nz);
						double vacancy_formation_energy_A;
						double vacancy_formation_energy_B;
						ecal_vacancy_formation(ilist[i], vacancy_formation_energy_A, vacancy_formation_energy_B);
						double equilibirum_vacancy_concentration_A = exp(-beta * vacancy_formation_energy_A);
						double equilibirum_vacancy_concentration_B = exp(-beta * vacancy_formation_energy_B);

						time_correction_factor_A = vacancy_concentration / equilibirum_vacancy_concentration_A;
						time_correction_factor_B = vacancy_concentration / equilibirum_vacancy_concentration_B;

						break;
					}
                    case 7:  if(iscaldsro) sro0= cal_sro(); 
                             genr(); N_genr ++;         //fh N_genr is for everytime an FP is generated
                             if(iscaldsro) acc_dsroG += cal_sro()-sro0;
                             break;
                    default: error(2, "(main) an unknown event type", 1, etype[i]);
                }

                goto actionDONE;
            }
		    
            acc_rate += rates[i]/sum_rates;
        }
    }
	
actionDONE:
    return 1.0/sum_rates;
}


//fh edit this to account for pbc
void class_events::actual_jumpV(int vid, int inbr){ // vcc id, neighbor id
    double sro0;
    if(iscaldsro) sro0= cal_sro();

    int xv= list_vcc[vid].x; // vcc position
	int yv= list_vcc[vid].y;
	int zv= list_vcc[vid].z;
    if(states[xv][yv][zv] != 0) error(2, "(actual_jumpV) the jumping vcc is not a vcc (type)", 1, states[xv][yv][zv]);

	//fh edit
	int crossed_x=0;
	int crossed_y=0;
	int crossed_z=0;
	int x= pbc_change(xv+v1nbr[inbr][0], nx, 0, crossed_x);
	int y= pbc_change(yv+v1nbr[inbr][1], ny, 1, crossed_y);
	int z= pbc_change(zv+v1nbr[inbr][2], nz, 2, crossed_z);
    if(states[x][y][z] != 1 && states[x][y][z] != -1) error(2, "(actual_jumpV) the jumping atom is not an atom (type)", 1, states[x][y][z]);
    
    int ltcp_v= xv*ny*nz+yv*nz+zv;
    //ltcp_array[ltcp][0]=1;
    ltcp_array[ltcp_v][1]=vacancy_correction_x;
    ltcp_array[ltcp_v][2]=vacancy_correction_y;
    ltcp_array[ltcp_v][3]=vacancy_correction_z;
	//fh edit
	int ltcp_a= x*ny*nz+y*nz+z;

	if (states[x][y][z] == -1)
	{
		//fh probably don't need this section
		solute_correction_x += crossed_x;
		solute_correction_y += crossed_y;
		solute_correction_z += crossed_z;

		
    	//ltcp_array[ltcp][0]=1;
   	 	ltcp_array[ltcp_a][1]+=crossed_x;
    	ltcp_array[ltcp_a][2]+=crossed_y;
    	ltcp_array[ltcp_a][3]+=crossed_z;

	}

	int typev=ltcp_array[ltcp_v][0];
	int corrxv=ltcp_array[ltcp_v][1];
    int corryv=ltcp_array[ltcp_v][2];
    int corrzv=ltcp_array[ltcp_v][3];

    ltcp_array[ltcp_v][0]=ltcp_array[ltcp_a][0];
    ltcp_array[ltcp_v][1]=ltcp_array[ltcp_a][1];
    ltcp_array[ltcp_v][2]=ltcp_array[ltcp_a][2];
    ltcp_array[ltcp_v][3]=ltcp_array[ltcp_a][3];
    	
    ltcp_array[ltcp_a][0]=typev;	
    ltcp_array[ltcp_a][1]=corrxv;
    ltcp_array[ltcp_a][2]=corryv;
    ltcp_array[ltcp_a][3]=corrzv;
    ////////////////////////////////////////////////////////
    
    if(states[x][y][z]==1)  Vja[0] ++; // track # of jumping atoms (see log file) 
    else                    Vja[1] ++;

	states[xv][yv][zv]= states[x][y][z];

    if(srf[x][y][z]){ // if jump into srf atom, becomes vacuum
        nV --;
        nM ++;

	    states[x][y][z]= 4;
        srf[x][y][z]= false;
        
        list_vcc.erase(list_vcc.begin()+vid);
    }
    else{
	    states[x][y][z]= 0;
    	list_vcc[vid].x= x;
    	list_vcc[vid].y= y;
    	list_vcc[vid].z= z;
    	if((x-xv)>nx/2) list_vcc[vid].ix --; if((x-xv)<-nx/2) list_vcc[vid].ix ++;
    	if((y-yv)>ny/2) list_vcc[vid].iy --; if((y-yv)<-ny/2) list_vcc[vid].iy ++;
    	if((z-zv)>nz/2) list_vcc[vid].iz --; if((z-zv)<-nz/2) list_vcc[vid].iz ++;
    
        if(list_vcc[vid].njump != -1) list_vcc[vid].njump ++;
        
        if(iscaldsro) acc_dsroV += cal_sro()-sro0; 
        recb_checkv(vid);
    }
    
    if(nM>0){ // update srf & creation
        srf_check(x, y, z);
        srf_check(xv, yv, zv);
        cvcc_rates += update_ratesC(x*ny*nz + y*nz + z);
        cvcc_rates += update_ratesC(xv*ny*nz+ yv*nz+ zv);
    }
}

void class_events::actual_jumpI(int iid, int inbr){
    int xi= list_itl[iid].x; // itl position
	int yi= list_itl[iid].y;
	int zi= list_itl[iid].z;
	if(2!=states[xi][yi][zi] && 3!=states[xi][yi][zi]) error(2, "(actual_jumpI) the jumping itl is not an itl (type)", 1, states[xi][yi][zi]);


	//fh edit
	int crossed_x = 0;
	int crossed_y = 0;
	int crossed_z = 0;
	
    int x, y, z;
    if(2==states[xi][yi][zi]){ // SIA jumps to 1st-nn

		x = pbc_change_itl(xi + v1nbr[inbr][0], nx, 0, crossed_x);
		y = pbc_change_itl(yi + v1nbr[inbr][1], ny, 1, crossed_y);
		z = pbc_change_itl(zi + v1nbr[inbr][2], nz, 2, crossed_z);
    }
    else{ // mixed itl jumps to 2nd-nn

		x = pbc_change_itl(xi + v2nbr[inbr][0], nx, 0, crossed_x);
		y = pbc_change_itl(yi + v2nbr[inbr][1], ny, 1, crossed_y);
		z = pbc_change_itl(zi + v2nbr[inbr][2], nz, 2, crossed_z);
    }

	//fh edit
	if (states[x][y][z] == -1)
	{
		solute_correction_x += crossed_x;
		solute_correction_y += crossed_y;
		solute_correction_z += crossed_z;
	}


    if(states[x][y][z] != 1 && states[x][y][z] != -1)  error(2, "(actual_jumpI) the jumping atom is not an atom (type)", 1, states[x][y][z]);
    
    if(2==states[xi][yi][zi])   Ija[0] ++;
    else                        Ija[1] ++;
	    
    if(states[x][y][z]==-1)
    {
        cout << "the AB is swapping with a B atom" << endl;
        exit(1);
    }
    states[x][y][z]= states[xi][yi][zi];
    states[xi][yi][zi]= 1; // change if BB is included

	
	list_itl[iid].x= x;
	list_itl[iid].y= y;
	list_itl[iid].z= z;
    list_itl[iid].dir=  inbr;

	
	if((x-xi)>nx/2) list_itl[iid].ix --; if((x-xi)<-nx/2) list_itl[iid].ix ++;
	if((y-yi)>ny/2) list_itl[iid].iy --; if((y-yi)<-ny/2) list_itl[iid].iy ++;
	if((z-zi)>nz/2) list_itl[iid].iz --; if((z-zi)<-nz/2) list_itl[iid].iz ++;

    if(nM>0){ // update srf & creation
        srf_check(x, y, z);
        srf_check(xi, yi, zi);
        cvcc_rates += update_ratesC(xi*ny*nz+ yi*nz+ zi);
        cvcc_rates += update_ratesC(x *ny*nz+ y *nz+ z);
    }

    if(trap_included && 3==states[x][y][z]) list_itl[iid].trapped= trap_check(x, y, z); // if trapping included, check if AB itl trapped
     
   
}
	
void class_events::create_vcc(int altcp, int mltcp)  //fh an atom goes from the surface into the vacuum while a vacuum site becomes a vacancy that moves
{ // vcc created from srf
	int xa= (int) (altcp/nz)/ny; // jumping atom
    int ya= (int) (altcp/nz)%ny; 
    int za= (int) altcp%nz;
	int xm= (int) (mltcp/nz)/ny; // to the vacuum site
    int ym= (int) (mltcp/nz)%ny; 
    int zm= (int) mltcp%nz;
	int ja= states[xa][ya][za];
    if(ja != 1 && ja != -1) error(2, "(create_vcc) an ja isnt 1 or -1", 1, ja);

	// Update states
	states[xm][ym][zm]= ja;
	states[xa][ya][za]= 0;

    // initialize the vcc in the list_vcc
	int vid= list_vcc.size();
	list_vcc.push_back(vcc());
	
	list_vcc[vid].x= xa;
	list_vcc[vid].y= ya;
	list_vcc[vid].z= za;
	list_vcc[vid].ix= 0;
	list_vcc[vid].iy= 0;
	list_vcc[vid].iz= 0;
    list_vcc[vid].njump= 0;

    nV ++;
    nM --;

    srf_check(xa, ya, za);
    srf_check(xm, ym, zm);
    cvcc_rates += update_ratesC(altcp);
    cvcc_rates += update_ratesC(mltcp);
    recb_checkv(vid);
}