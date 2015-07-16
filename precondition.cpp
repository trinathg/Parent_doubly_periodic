#include "headers.h"
#include "init_2.h"
#include "declarations.h"

//Function that takes in a residual vector and preconditions with MG method

vec precond(vector<mg_grid>& level, vec x, vec rhs, pbcs& pbc)
{
	int no_cycles = 2; 
	
	mg_pre_initialize(level,x,rhs);  //Most important function of mapping the shortened vector to the field with pinned points included
	
	/****Only V-cycle preconditioning*****/ 
	
	for(int no_v=0;no_v<no_cycles;no_v++)
	{
		v_cycle(level,0,pbc); 
	}
	
	/****FMG preconditioning*****/
	
	//pre_fmg(level,pbc); 		
	
	return mg_post_precond(level);
}

/*********************************************************/ 

void mg_pre_initialize(vector<mg_grid>& level,vec & X, vec &rhs)
{
	int ind, ind_m, i_m, j_m; 
	int lev=0; 			
	int sp = pow(2,lev);
	
	int st_inx=0, st_iny=0;  	//True indices of the computational domain
	int en_inx = nx_per, en_iny = ny_per;   //True indices of the computational domain
		
	int nx_sol = (nx_per/pow(2,lev)), ny_sol = (ny_per/pow(2,lev)); //No of spacings of the coarsest grid
	int tot_p_sol = (nx_sol+1)*(ny_sol+1); //No.of points that are solved for directly. The top most and the farthest right line of points are removed due to periodicity
		 
	int str_m = nx_sol+1; //Stride for the coarsest grid walk (also reduces by 1 for the periodic case)

	for(int j=0;j<=ny_per;j=j+sp)
	{
		for(int i=0;i<=nx_per;i=i+sp)
		{
			i_m = (i/sp); 
			j_m = (j/sp); 
			
			ind = i + j*str_x; 				
			ind_m = i_m + j_m*str_m; 
						
			//cout<<"i = "<<i<<"\nj = "<<j<<"i_m = "<<i_m<<"\nj_m = "<<j_m<<"\n";
					
			if(ind_m<=tot_p_sol-2)
			{
				//level[lev].phi_s[ind] = X(ind_m);
				level[lev].phi_s[ind] = 0.0;
				level[lev].rhs[ind] = rhs(ind_m); 
			}
			else
			{
				level[lev].phi_s[ind] = 0.0; 	
				level[lev].rhs[ind] = 0.0; 		
			}
		}
	}		
	
	mg_bcs_neu(level,lev);				
}

/***************************************************************************/ 

vec mg_post_precond(vector<mg_grid>& level)
{	
	int ind, ind_m, i_m,j_m;
	 
	int lev=0; 				
	int sp = pow(2,lev);
	
	int st_inx=0, st_iny=0;  	//True indices of the computational domain
	int en_inx = nx_per, en_iny = ny_per;   //True indices of the computational domain
		
	int nx_sol = (nx_per/pow(2,lev)), ny_sol = (ny_per/pow(2,lev)); //No of spacings of the coarsest grid
	
	int tot_p_sol = (nx_sol+1)*(ny_sol+1); //No.of points that are solved for directly. The top and the farthest right points are removed due to periodicity
		 
	int str_m = nx_sol+1; //Stride for the coarsest grid walk (also reduces by 1 for the periodic case)

	vec pcond_vec(tot_p_sol); 
	
	for(int j=0;j<=ny_sol;j=j+sp)
	{
		for(int i=0;i<=nx_sol;i=i+sp)
		{
			ind = i + j*str_x; 							
			
			i_m = (i/sp);
			j_m = (j/sp); 
			
			ind_m = i_m + j_m*str_m; 
			
			pcond_vec(ind_m) = level[lev].phi_s[ind];	//Untrimmed preconditioned vector
		}
	}
	
	return pcond_vec.submat(0,0,tot_p_sol-2,0); 
} 

/***************************************************************************/
//This function adjusts the preconditioned residual such that the last point in the 
//periodic computable domain is zero because that value is globally taken to be zero
//even when it is pinned in the CG method 

void adjust_precond(vector<mg_grid>& level,vec & zi)
{
	int str_m = nx_per+1; 
	int ind; 
	
	for(int j=0;j<ny_per;j++)
	{
		for(int i=0;i<nx_per;i++)
		{
			ind = i + j*str_x; 
						
			zi(i+j*str_m) = zi(i+j*str_m) - level[0].phi_s[nx_per+ny_per*str_x]; 							
		}	
	}
}

/***************************************************************************/
