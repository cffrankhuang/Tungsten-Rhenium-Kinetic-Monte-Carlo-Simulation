#ifndef KMC_EVENTS_INCLUDED
#define KMC_EVENTS_INCLUDED
#include <vector>
#include <iostream>
#include <cmath>
#include <unordered_map>

using namespace std;

class class_events{
	public:
		class_events(){
			cout << "##Generation parameters (rate_genr) " << rate_genr << " (damage/s)" << endl;
			cout << "##Recombination parameters: 3rd nearest-neighbor distance (FIXED in SURFACE simulations) " << endl;
		
            cvcc_rates= init_ratesC();
            N_nediff= 0;
        }
		
        // variables
        double cvcc_rates;  // total vcc creation rate  //frank: rate of vacancy formation from surface
        int N_nediff;       // times that neg ediff

        // functions
		double main(double &time_correction_factor_A, double &time_correction_factor_B);
        double ecal_range(int xlo=0, int xhi=nx-1, int ylo=0, int yhi=ny-1, int zlo=0, int zhi=nz-1, bool is_corr= false);
	
	private:
        // ********************  VARIABLES  ******************** //
        // vcc creation at surface //
#define N_NOCVCC 3 // in update cvcc, not put into list if # of bonds are less
        struct cvcc_info{
            long long int step; // modified step
            vector <int> mltcp;
            vector <double> rates;
        };
        unordered_map <int, struct cvcc_info> cvcc;
        // ********************  VARIABLES  ******************** //
        
		///// functions of energy calculation /////
		double ecal_bond(int x1, int y1, int z1, int x2, int y2, int z2) const; 
		void ecal_vacancy_formation(int vacancy_index, double &vacancy_formation_A, double &vacancy_formation_B) const;    //fh function
        double ecal_nonb(int x1, int y1, int z1, int x2, int y2, int z2) const; // cal non-broken A-B bonds 
        double ecal_sp(int stateA1, int inbr, int i, int j, int k) const;       // calculate saddle-point e (vcc)

		///// functions of events /////
		void actual_jumpI(int iid, int jatom);
        void actual_jumpV(int vid, int jatom);
		void genr();
        void create_vcc  (int altcp, int mltcp);
		
        ///// functions for rate calculations) /////
		double cal_ratesI  (vector <int> &etype, vector <double> &rates, vector <int> &ilist, vector <int> &jatom);
//		double cal_ratesV  (vector <int> &etype, vector <double> &rates, vector <int> &ilist, vector <int> &jatom); // replaced by sp
        double cal_ratesVsp(vector <int> &etype, vector <double> &rates, vector <int> &ilist, vector <int> &jatom);
        double update_ratesC(int ltcp_in, bool is_recb= false); // update rates of vcc creation at srf
        double init_ratesC();
		
        ///// functions for recombination /////
        void rules_recb(int ii, int xi, int yi, int zi, int iv, int xv, int yv, int zv);
        bool recb_checki(int id);
        bool recb_checkv(int id);
        void recb_check_ecal(bool isitl, double& minE, vector<vector<int>>& list_recb, int xi, int yi, int zi, int stateI, int xv, int yv, int zv, int stateV);
        void srf_check(int i, int j, int k);
        void sink(bool isvcc, int index); // execute the sink
        bool trap_check(int i, int j, int k);
};

#endif // KMC_EVENTS_INCLUDED
        


