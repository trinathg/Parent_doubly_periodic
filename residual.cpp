#include "headers.h"
#include "init_2.h"
#include "declarations.h"

double div_calc(vector<vertex>& node, vector<fval>& fvar)
{
	double max_div=0.0; 
	int ind; 
	
	tdma1x(fvar,0);
	tdma1y(fvar,1); 
	
	for(int j=0;j<=ny_per;j++)
	{
		for(int i=0;i<=nx_per;i++)
		{
			ind = i + j*str_x; 		
			
			fvar[ind].div = fvar[ind].ux[0] + fvar[ind].uy[1]; 
			
		 	if( abs(fvar[ind].div)>max_div ) max_div = abs(fvar[ind].div); 		
		}
	}
	
	return max_div; 
}

/*********************************************************************************************/
/****The L2 norm of the residuals from time step to time step for all the basic variables******/
/****The file pointer is passed to the function from the solver function*/

void global_residuals(vector<fval>& fvar, ofstream& fp_res, int ite)
{
  int ind; 
  double time = ite*dt; //Computing the time at which these are computed
  double u_res=0.0,v_res=0.0, p_res=0.0; 
  
  for(int j=0;j<=ny_per;j++)
  {
    for(int i=0;i<=nx_per;i++)
    {
	ind = i + j*str_x;   
	
	u_res = u_res + (fvar[ind].u[0] - fvar[ind].u0[0])*(fvar[ind].u[0] - fvar[ind].u0[0]); 
	v_res = v_res + (fvar[ind].u[1] - fvar[ind].u0[1])*(fvar[ind].u[1] - fvar[ind].u0[1]); 
	p_res = p_res + (fvar[ind].u[2] - fvar[ind].u0[2])*(fvar[ind].u[2] - fvar[ind].u0[2]);		
    }  
  }
  
  fp_res<<time<<"		"<<u_res<<"	"<<v_res<<"	"<<p_res<<"\n"; 
}

/*********************************************************************************************/
