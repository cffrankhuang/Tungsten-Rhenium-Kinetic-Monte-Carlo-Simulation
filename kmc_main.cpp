#include <iostream>
#include <fstream>
#include <ctime>
#include "kmc_par.h"
#include "kmc_global.h"
#include "kmc_initial.h"
#include "kmc_events.h"
using namespace std;

int main(int nArg, char *Arg[]){
	int t0cpu= time(0);

	long long int ts_bg;   //background timestep
	double time_bg;       //background time
	double time_bg_corrected_A;  //background corrected time; 
	double time_bg_corrected_B;  //background corrected time; 

	cout << "########## Initializing System... ##########" << endl;
	//class_initial init(ts_bg, time_bg, nArg, Arg[1]);     //fh: is restarting now
	class_initial init(ts_bg, time_bg, time_bg_corrected_A,time_bg_corrected_B, nArg, Arg[1], Arg[2]);     //fh: new initial function
	
	cout << "\n########## Initializing Events ... ##########" << endl;
	class_events events; 

	cout << "\n########## The Simulation Begins !! ##########" << endl;
	timestep=  ts_bg; totaltime= time_bg;
	totaltime_vacancy_adjust_A = time_bg_corrected_A;              //fh addition
	totaltime_vacancy_adjust_B = time_bg_corrected_B;              //fh addition

    double init_sro= cal_sro();

    if(0==timestep){ fprintf(out_engy, "%lld %e %f\n", timestep, totaltime, events.ecal_range()); 
                     fprintf(out_sro,  "%lld %e %e %e %e %e %e %e\n", timestep, totaltime, init_sro, 0.0, 0.0, 0.0, 0.0, 0.0); }
    cout << "TIMESTEP() TIME(s) GENR()	NA() NB()	NV() NAA() NAB() NBB()	AJUMP_V% AJUMP_I%";
	printf("\n%lld %.10e %d     %d %d     %d %d %d %d     %f %f", timestep, totaltime, 0, nA, nB, nV, nAA, nAB, nBB, 0.0, 0.0);
    int N_0def= 0;    //fh: defines the event when a generation occurs when there is no other
    int N_conf= 0;

    

	while((totaltime<= time_bg+par_time) && (timestep != ts_bg+par_step)){
		frenkel_pairs_generated = 0;   //fh: addition

        	timestep ++;

		double time_correction_vacancy_A=1.0;
		double time_correction_vacancy_B = 1.0;
		
		// CALCULATIONS
		double dt= events.main(time_correction_vacancy_A, time_correction_vacancy_B); // Defect jumps
        	totaltime += dt;
		totaltime_vacancy_adjust_A += dt * time_correction_vacancy_A;
		totaltime_vacancy_adjust_B += dt * time_correction_vacancy_B;

		
		
		// OUTPUT DATA
		if(0==timestep%step_log)
		{
            cout << endl;

			printf("%lld %.10e %d     %d %d     %d %d %d %d     %f %f", timestep, totaltime, N_genr, nA, nB, nV, nAA, nAB, nBB, 
                    ((double) Vja[0])/(Vja[0]+Vja[1]), ((double) Ija[0])/(Ija[0]+Ija[1]));
			if(N_0def != 0)
			{
				cout << "  *** 0-defect genr: " << N_0def << " ***"; 
                N_0def= 0;
			}
            if(events.N_nediff != 0)  //times you have negative ediff
			{
				cout << "  *** neg e: " << events.N_nediff << " ***"; 
                events.N_nediff= 0;
            }
			


			//fh might want to implement feature to print out amount of FP pairs generated --> N_genr


            fflush(stdout);
		}

		if(0==timestep%step_conf || totaltime>=(N_conf+1)*time_conf){
			if(0==timestep%step_conf) write_conf(1);
            else                      write_conf(2);
            cout << "   <Output conf files at: " << timestep << ">";
            N_conf= (int) (totaltime / time_conf);
            
    		//fprintf(correction_factors_file, "%d %d %d %d %d\n", type, ltcp, i, j, k);
			fprintf(correction_factors_file, "%lld %d %d %d %d %d %d %d %d %d\n", timestep, vacancy_correction_x, vacancy_correction_y, vacancy_correction_z, solute_correction_x, solute_correction_y, solute_correction_z, itl_correction_x, itl_correction_y, itl_correction_z);
			fflush(correction_factors_file);


        }

		if(0==timestep%step_out){
			fprintf(out_engy, "%lld %e %f\n", timestep, totaltime, events.ecal_range()); fflush(out_engy);
            double sro= cal_sro();
            double acc_dsroI= sro - init_sro - acc_dsroV - acc_dsroRi - acc_dsroRv - acc_dsroG;
			fprintf(out_sro,  "%lld %e %e %e %e %e %e %e\n", timestep, totaltime, sro, acc_dsroI, acc_dsroV, acc_dsroRi, acc_dsroRv, acc_dsroG); fflush(out_sro);
            
            if(par_isOUTrestart) write_conf(3);
        }
		
		if(0==timestep%step_his){
 			//write_hissol(); // write sol first to construct list_srf
			//write_hisdef(); // and then output his_srf here
			
			

			if (par_typeD == 0)         //fh if you only have vacancy
			{
				write_hisdis(); //fh: output file for r2 purposes
				//write_hisdis_correction();           //fh addition 
			}
			else if (par_typeD == 3)
			{
				write_hisdis_itl();      //fh addition temporary
			}
		}
	}

	// finalizing
	if(timestep%step_log != 0) printf("\n%lld %f %d	%d %d	%d %d %d %d ", timestep, totaltime, N_genr, nA, nB, nV, nAA, nAB, nBB);
	write_conf(1); cout << "<Output conf files at: " << timestep << ">";
    cout << "\nAccumulative sro change from vcc: " << acc_dsroV << endl;

    FILE* correction_factors_file2 = fopen("correction_factors2.txt", "a");
    //fprintf(correction_factors_file, "%d %d %d %d %d\n", type, ltcp, i, j, k);
	fprintf(correction_factors_file2, "%lld %d %d %d %d %d %d %d %d %d\n", timestep, vacancy_correction_x, vacancy_correction_y, vacancy_correction_z, solute_correction_x, solute_correction_y, solute_correction_z, itl_correction_x, itl_correction_y, itl_correction_z);
	fflush(correction_factors_file2);
	
	
	

	int tfcpu= time(0);
	cout << "**** The simulation is done! Total CPU time: " << tfcpu - t0cpu << " secs ****" << endl;

	return 0;
}
