#include <iostream>
#include <cmath>
#include <vector>
#include "kmc_global.h"
#include "kmc_events.h"
#include "kmc_par.h"

using namespace std;

//fh function that performs recombination
//fh i corresponds to interstitial, v corresponds to vacancy probably
void class_events::rules_recb(int ii, int xi, int yi, int zi, int iv, int xv, int yv, int zv) 
{ // execute the recombination: vcc or vacuum
    if(-1==ii) //fh: search for coordinates of i
	{ // if ii not know, search it
        for(int a=0; a<list_itl.size(); a ++) 
            if(list_itl[a].x==xi && list_itl[a].y==yi && list_itl[a].z==zi) ii= a;
        if(-1==ii) error(2, "(rules_recb) cant find the itl");
    }
    if(-1==iv && 0==states[xv][yv][zv])              //fh search for coordinates of v
	{
        for(int a=0; a<list_vcc.size(); a ++)
            if(list_vcc[a].x==xv && list_vcc[a].y==yv && list_vcc[a].z==zv) iv= a;
        if(-1==iv) error(2, "(rules_recb) cant find the vcc");
    }

    if(-1==states[xv][yv][zv])   //fh perform AA +B to A + AB event
	{ // AA+B->A+AB
        if(states[xi][yi][zi] != 2) error(2, "(rules_recb) a recb with B atom not AA itl (type)", 1, states[xi][yi][zi]);
        nAA --; nB --;
        states[xi][yi][zi]= 1;
        states[xv][yv][zv]= 3;
        nAB ++; nA ++;
        list_itl[ii].x= xv;
        list_itl[ii].y= yv;
        list_itl[ii].z= zv;
        recb_checki(ii);          //fh perform event and then check for recombination again
    }
    else
	{ // recb: I+V or I+M
        if(4==states[xv][yv][zv]) nM --;        //fh: if the state is a vacuum state
        else if(0==states[xv][yv][zv])          //fh: if the state is a vacancy
		{
            nV --;
            list_vcc.erase(list_vcc.begin()+iv);
        }
        else error(2, "(rules_recb) the vcc isnt 0 or 4", 1, states[xv][yv][zv]);

        double e1, e2;
        switch(states[xi][yi][zi])    //fh switch statement to see whether you have AA or AB
		{
		    case 2: // AA
			    nAA --; nA +=2;
	            states[xi][yi][zi]= 1;
	            states[xv][yv][zv]= 1;       //fh perform a swap
			    break;
            case 3: // AB
			    nAB --; nA ++; nB ++;
	            states[xi][yi][zi]= 1;
	            states[xv][yv][zv]=-1;
                e1= ecal_bond(xi, yi, zi, xv, yv, zv) + ecal_nonb(xi, yi, zi, xv, yv, zv);     //calculate energy
	            states[xi][yi][zi]=-1;
	            states[xv][yv][zv]= 1;
                e2= ecal_bond(xi, yi, zi, xv, yv, zv) + ecal_nonb(xi, yi, zi, xv, yv, zv);    //calculate energy
                if(e1<e2)
				{ // comparing energy to decide whether A or B jumps to V      
					//fh: choose A or B that corresponds to lowest energy
	                states[xi][yi][zi]= 1;
	                states[xv][yv][zv]=-1;
                }
                break;
		    default: error(2, "(rules_recb) an unknown itl type", 1, states[xi][yi][zi]);
	    }
			
	    list_itl.erase(list_itl.begin()+ii);           //erase interstitial from interstitial list
    
        if(nM !=0)                 //if you still have vacuum states
		{
            srf_check(xi, yi, zi);                         //check for recombination
            srf_check(xv, yv, zv);                            //check for recombination
            cvcc_rates += update_ratesC(xi*ny*nz+yi*nz+zi, true);       //update rates for vacancy formation
            cvcc_rates += update_ratesC(xv*ny*nz+yv*nz+zv, true);
        }
    }
}

bool class_events::recb_checki(int id)   //function to perform recombination of interstitial
{
    int i= list_itl[id].x;             //gather x,y,z coordinates of this interstitial
    int j= list_itl[id].y;
    int k= list_itl[id].z;
	//removed sink
	/*
    if(i==x_sink){ // check if at sink   // x_sink is the key parameter here
        sink(false, id);
        return true;
    }
	*/
    int stateI= states[i][j][k];      //create variable for interstitial state
    
    vector<vector<int>> list_recb;     //recombination list
    vector<vector<int>> list_AAtoAB;   //AA to AB list
    double minE= 999999999; // if multiple cases, choose a minE case
    int x, y, z;
	for(int a=0; a<n123nbr; a ++){ // check recb         //loop through 1st, 2nd, 3rd nearest neighbors of interstitial
        x= pbc(i+v123nbr[a][0], nx);    
		y= pbc(j+v123nbr[a][1], ny);
		z= pbc(k+v123nbr[a][2], nz);
        if(0==states[x][y][z] || 4==states[x][y][z]){     
            int stateV= states[x][y][z];                  //create variable for the "vacancy" state
            recb_check_ecal(true, minE, list_recb, i, j, k, stateI, x, y, z, stateV);  //calculate energy of event
            // minE, list_recb are references, change in the func

            states[i][j][k]= stateI;           //reseting the staes
            states[x][y][z]= stateV;           //reseting the states
        }
        if(a<n1nbr && ( -1==states[x][y][z] && 2==stateI ))     //if this is a first nearest neighbor 
			                                                                 //and you have a B neighbor and AA interstitial
			list_AAtoAB.push_back({x, y, z}); // AA+B->A+AB
    }

    if(list_recb.size() !=0)       //after looping through the three nearest neighbors, if recombination list is not empty
	{ //recb
        double sro0;                     //set variable for sro
        if(iscaldsro) sro0= cal_sro();    //calculate sro
        
        int ran= (int) ( ran_generator()*list_recb.size() );    //make a random number within the recombination list
        x= list_recb[ran].at(0);      //at is another way to index vectors
        y= list_recb[ran].at(1);
        z= list_recb[ran].at(2);
        rules_recb(id, i, j, k, -1, x, y, z);    //this seems to be the function that performs recombination
		// vccID is unknown, give -1
        
        if(iscaldsro) acc_dsroRi += cal_sro()-sro0;      //calculate the change in sro
        return true;
    }
    else if(list_AAtoAB.size() !=0)  //recombination with vacancy seems to take precedence over AA to AB transition
	{ // AA+B->A+AB
        double sro0;                             //create variable for sro
        if(iscaldsro) sro0= cal_sro();           //calculate change in sro
        
        int ran= (int) ( ran_generator()*list_AAtoAB.size() );        //choose a random event from the event list
        x= list_AAtoAB[ran].at(0);                //index
        y= list_AAtoAB[ran].at(1);
        z= list_AAtoAB[ran].at(2);
        rules_recb(id, i, j, k, -1, x, y, z);      //this function seems to perform recombination
		// vccID is unknown, give -1
        
        if(iscaldsro) acc_dsroRi += cal_sro()-sro0;         //calculate change in sro
        return true;                                  
    }
    else return false;                            //return false if you do not perform recombination
}

bool class_events::recb_checkv(int id)  //vacancy recombination
{
    int i= list_vcc[id].x;        //gather coordinates of the vacancy
    int j= list_vcc[id].y;
    int k= list_vcc[id].z;

	/*       //removed sink
    if(i==x_sink){ // check if at sink  //x_sink is the key parameter to determine if you want to make it into sink
        sink(true, id);
		cout << "sink recombination occurred" << endl; 
        return true;
    }
	*/
    int stateV= states[i][j][k];   //create variable to capture vacancy state
    
    vector<vector<int>> list_recb;   //create recombination list
    double minE= 999999999;          //minimum energy
    int x, y, z;
	for(int a=0; a<n123nbr; a ++)     //loop through the first 3 neighbors of this vacancy
	{ // check recb
        x= pbc(i+v123nbr[a][0], nx);    //gather coordinates of the neighbor
		y= pbc(j+v123nbr[a][1], ny);
		z= pbc(k+v123nbr[a][2], nz);

        if((2==states[x][y][z] || 3==states[x][y][z])) //if you have an AA or an AB neighbor
		{
            int stateI= states[x][y][z];                 //create variable to describe the interstitial neighbor
            recb_check_ecal(false, minE, list_recb, x, y, z, stateI, i, j, k, stateV);    //check the energy of this event
            // minE, list_recb are references, change in the func

            states[x][y][z]= stateI;                   //reset the variables to original state
            states[i][j][k]= stateV;
        }
    }

    if(list_recb.size() !=0)  //if after looping through the 3 nearest neighbors, you have recombination events
	{ //recb
        double sro0;                   //calculate sro
        if(iscaldsro) sro0= cal_sro();   
        
        int ran= (int) ( ran_generator()*list_recb.size() );    //randomly select one event
        x= list_recb[ran].at(0);
        y= list_recb[ran].at(1);
        z= list_recb[ran].at(2);
        rules_recb(-1, x, y, z, id, i, j, k);       //use this function to perform recombination
		// itlID is unknown, give -1
        
        if(iscaldsro) acc_dsroRv += cal_sro()-sro0;        //calculate the change in sro from this event
        return true;
    }
    else return false;                             //return false if you did not recombine your vacancy
}

void class_events::recb_check_ecal(bool isitl, double& minE, vector<vector<int>>& list_recb, int xi, int yi, int zi, int stateI, int xv, int yv, int zv, int stateV){
    double e0, ediff;              //create variable to describe initial energy and change in energy
    switch(stateI){ // perform image recb to cal ediff
        case 2:
            e0= ecal_bond(xi, yi, zi, xv, yv, zv);
            states[xi][yi][zi]= 1; 
            states[xv][yv][zv]= 1;
            ediff= ecal_bond(xi, yi, zi, xv, yv, zv) - e0;
            if(abs(ediff-minE)<1e-7){
                if(isitl)   list_recb.push_back({xv, yv, zv});
                else        list_recb.push_back({xi, yi, zi});
            }
            else if((ediff-minE)<0){
                list_recb.clear();
                if(isitl)   list_recb.push_back({xv, yv, zv});
                else        list_recb.push_back({xi, yi, zi});
                minE= ediff;
            }
                    
            break;
        case 3:
            e0= ecal_bond(xi, yi, zi, xv, yv, zv) + ecal_nonb(xi, yi, zi, xv, yv, zv);
            states[xi][yi][zi]=  1; 
            states[xv][yv][zv]= -1;
            ediff= ecal_bond(xi, yi, zi, xv, yv, zv) + ecal_nonb(xi, yi, zi, xv, yv, zv) - e0;
            if(abs(ediff-minE)<1e-7){
                if(isitl)   list_recb.push_back({xv, yv, zv});
                else        list_recb.push_back({xi, yi, zi});
            }
            else if((ediff-minE)<0){
                list_recb.clear();
                if(isitl)   list_recb.push_back({xv, yv, zv});
                else        list_recb.push_back({xi, yi, zi});
                minE= ediff;
            }
            
            states[xi][yi][zi]= -1; 
            states[xv][yv][zv]=  1;
            ediff= ecal_bond(xi, yi, zi, xv, yv, zv) + ecal_nonb(xi, yi, zi, xv, yv, zv) - e0;
            if(abs(ediff-minE)<1e-7){
                if(isitl)   list_recb.push_back({xv, yv, zv});
                else        list_recb.push_back({xi, yi, zi});
            }
            else if((ediff-minE)<0){
                list_recb.clear();
                if(isitl)   list_recb.push_back({xv, yv, zv});
                else        list_recb.push_back({xi, yi, zi});
                minE= ediff;
            }
                    
            break;
    }
}

void class_events::srf_check(int i, int j, int k)
{  
	// when vacuum changed, check if srf array changes
    for(int a= -1; a<n1nbr; a ++)    //loop through neighbors of the inputted specie
	{ // here we check if some surface atoms become non-surface atoms (hv no bond with vacuum)
        int x, y, z;
        if(-1==a){x= i; y= j; z= k;} // check itself
        else{
            x= pbc(i+v1nbr[a][0], nx);
            y= pbc(j+v1nbr[a][1], ny);
            z= pbc(k+v1nbr[a][2], nz);
        }
        
        if(states[x][y][z] != 1 && states[x][y][z] != -1) //if the neighbor is not an A or a B atom, set to not be surface atom
		{
            srf[x][y][z]= false;
            continue;
        }

        bool is_srf= false;
        for(int b=0; b<n1nbr; b ++)  //loop through first neighbors of the atom
		{
            int d= pbc(x+v1nbr[b][0], nx);
            int e= pbc(y+v1nbr[b][1], ny);
            int f= pbc(z+v1nbr[b][2], nz);
            if(4==states[d][e][f])   //if you detect a vacuum state, this atom is a surface atom
			{
                is_srf= true;
                break;
            }
        }

        srf[x][y][z]= is_srf; 
    }
}

void class_events::sink(bool isvcc, int index){ // execute the sink
	int x, y, z;
	double ran;
	
	if(isvcc)  //if this is a vacancy
	{
		nV --;                    //reduce vacancy number by 1
		x= list_vcc[index].x;     //gather x,y,z coordinates of the vacancy
		y= list_vcc[index].y;
		z= list_vcc[index].z;
		list_vcc.erase(list_vcc.begin()+index);  //erase the vacancy from vacancy list

        ran= ran_generator();                    //generate random number
		if(list_sink.size() != 0) // if atoms in sink, use them  .if the sink list is not zero
		{ 
            int i= (int) ( ran*list_sink.size() );             //choose a random sink in the sink list
			states[x][y][z]= list_sink[i];                     //set states to that sink state
            if(list_sink[i] != 1 && list_sink[i] != -1) error(2, "a atom from sink isnt atom", 1, list_sink[i]);
			list_sink.erase(list_sink.begin()+i);
		}
		//else if sink list is blank
		else states[x][y][z]= (ran<par_compA) ? 1:-1;    

		if(1==states[x][y][z])  nA ++;                //if the location is A atom, increase number of A atoms
		else                    nB ++;                //if the location is B atom, increase number of B atoms
	}
	else{     //if this is an interstitial
        x= list_itl[index].x;        //grab coordinates of the interstitial
        y= list_itl[index].y;
        z= list_itl[index].z;
		list_itl.erase(list_itl.begin()+index);    //erase the interstitial
		
		switch(states[x][y][z]){        //switch statement to see if AA or AB
			case  2:
				nAA --;
				states[x][y][z]= 1; nA ++;  //make it into an A atom
				list_sink.push_back(1);     //add to the sink list
				break;
			case  3:
				nAB --;
				states[x][y][z]= (ran_generator()>0.5) ? 1:-1;        //randomly choose to make it an A or B atom
                if(1==states[x][y][z]){ nA ++; list_sink.push_back(-1);}   //assign the type to the sink list
                else                  { nB ++; list_sink.push_back( 1);}
				break;
			default: error(2, "(sink) an unknown type", 1, states[x][y][z]);
		}
	}
}

bool class_events::trap_check(int i, int j, int k)                    //check if AB interstitial is trapped
{ // check if AB itl trapped
    if(states[i][j][k] != 3) error(2, "(trap_check) input not AB itl", 1, states[i][j][k]);
    if(! trap_included)      error(2, "(trap_check) trap_included not turn on");
    if(temp>1000)            error(2, "(trap_check) T >1000k, should not trapped", 1, temp);

    for(int a= 0; a<n1nbr; a ++) //loop through first neighbors of AB interstitial
	{
        int x= pbc(i+v1nbr[a][0], nx);
        int y= pbc(j+v1nbr[a][1], ny);
        int z= pbc(k+v1nbr[a][2], nz);
        
        if(3==states[x][y][z] || -1==states[x][y][z]) return true; // if near by B or AB, trapped
		                                                        //this implies that AB being around B or AB traps it
    }

    return false; // not trapped
}


