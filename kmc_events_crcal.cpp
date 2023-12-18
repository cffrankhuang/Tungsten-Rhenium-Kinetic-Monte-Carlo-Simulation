#include <iostream>
#include <vector>
#include "kmc_global.h"
#include "kmc_events.h"
#include <unordered_map>

using namespace std;

double class_events::update_ratesC(int ltcp_in, bool is_recb){ // remove rate info of them and rebuilt it
	//frank: this function updates rates of recombination at surface
	//frank: ltcp_in is some kind of index for the lattice
	//frank: this function return rate change from surface hop event

	double rate_change= 0;

    int i0= (int) (ltcp_in/nz)/ny;    //fh custom algorithms to find x,y,z indices of the vacancies
    int j0= (int) (ltcp_in/nz)%ny;
    int k0= (int)  ltcp_in%nz;

	for(int a= -1; a<n1nbr; a ++) //fh loop through first nearest neighbors of this specie
	{ // loop through itself & nbrs
        int i, j, k, ltcp;
        if(-1==a){ i= i0; j= j0; k= k0; ltcp= ltcp_in;}     //fh i,j,k are indices of the neighbor
        else{
	        i= pbc(i0+v1nbr[a][0], nx);
	        j= pbc(j0+v1nbr[a][1], ny);
	        k= pbc(k0+v1nbr[a][2], nz);
    	    ltcp= i*ny*nz + j*nz + k;
        }

        if(cvcc.find(ltcp) != cvcc.end())   //fh look through the cvcc vector to find of the custom ltcp index corresponds to anything
		{
            if(! is_recb && cvcc[ltcp].step == timestep) continue; // if has been modified, dont do again

            for(int b=0; b<cvcc[ltcp].rates.size(); b ++)   //fh loop through rates list in the cvcc
															//fh I'm assuming the rates list is for the possible jumps of the 
				rate_change -= cvcc[ltcp].rates[b];
            cvcc.erase(ltcp); // re-calculate it
        }
    
        if(srf[i][j][k])   //fh make sure to go check how srf array is initiated
		{
            if(states[i][j][k] != 1 && states[i][j][k] != -1) error(2, "(update_ratesC) a srf atom not A nor B", 1, states[i][j][k]);
            cvcc[ltcp].step= timestep;                //fh set the timestep to the iteration we're at
            
            double mu= (1==states[i][j][k]) ? muvA:muvB;  //fh set mu to to whichever atom
            double em= (1==states[i][j][k]) ? emvA:emvB;   //fh set em to whichver atom

            double rate_temp= 0; // temporary rate
	        bool isne= false;    // is neg energy
            int AM= 0;           // count of A-M bonds
            for(int b= 0; b<n1nbr; b ++) //fh loop through the first neighbor atoms of this surface site
			{
	            int x= pbc(i+v1nbr[b][0], nx);
	            int y= pbc(j+v1nbr[b][1], ny);
	            int z= pbc(k+v1nbr[b][2], nz);
    	        int ltcp2= x*ny*nz + y*nz + z;    //fh get a custom index of the first neighbor

                if(4==states[x][y][z])    //fh if the first neighbor is a surface "atom"
				{
                    AM ++;

                    states[x][y][z]= 0; // considering as V+atom->atom+V, not M+atom->atom+V   //fh consider surface as vacancy
                    double e0= ecal_bond(x, y, z, i, j, k) + ecal_nonb(x, y, z, i, j, k);      //fh calculate energy before atom and "vacancy" swap
				    states[x][y][z]= states[i][j][k];
                    states[i][j][k]= 0;                       //fh switch their states
		            double ediff= ecal_bond(x, y, z, i, j, k) + ecal_nonb(x, y, z, i, j, k) - e0;     //fhcalculate energy diff
				    states[i][j][k]= states[x][y][z];
				    states[x][y][z]= 4;
                    
                    ediff += em; // ediff becomes ediff+em
                    if(ediff<0){
                        isne= true;
                        ediff= 0;
                    }
                    
				    cvcc[ltcp].mltcp.push_back(ltcp2);                  //fh in the custom index list, push back coordinate of neighbor
                    cvcc[ltcp].rates.push_back(mu * exp(-beta*ediff));   //fh in the rates list, push back rate
				
				    rate_temp += cvcc[ltcp].rates.back();
                }
		    }

            if(0==AM) error(2, "(cvcc) an srf atom does not bond with vacuum");
            else if(AM >= N_NOCVCC) cvcc.erase(ltcp); // either extrusion or intrusion
            else{
                rate_change += rate_temp;
                if(isne) cout << "*** neg e in cvcc ***";
            }
        }
    }

	return rate_change;
}

double class_events::init_ratesC(){
	double sum_rate= 0;

	for(int i=0; i < nx; i ++){
	    for(int j=0; j < ny; j ++){
	        for(int k=0; k < nz; k ++){
		        if(srf[i][j][k]) sum_rate += update_ratesC(i*ny*nz+j*nz+k);
    }}}

	return sum_rate;
}
