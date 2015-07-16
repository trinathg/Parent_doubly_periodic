#include "headers.h"
#include "init_2.h"
#include "declarations.h"

void make_mesh()
{

	int ind;  // The global index
	double start_time, end_time;
	int precond_count=0, restart_count; 

	vector<vertex> node (tot_p);
	vector<fval> fvar (tot_p); 
	vector<mg_grid> level(mg_levels);
	pbcs pbc; 

	for(int j=0;j<=ny;j++) 		//Grid generation 
	{
		for(int i=0;i<=nx;i++)
		{
			ind = i + j*str_x; 
			node[ind].x[0] = i*dx;
			node[ind].x[1] = j*dy;   
		}
	} 
	
	omp_set_num_threads(num_threads); 
  	
	
	if(restart==0)
        {
	  initialize(node,fvar);
	  mg_initialize(level,node);			
          restart_count=0;
        }
        else
        {
          read_restart(fvar,restart_count);
          mg_restart(fvar,level);
        }
        
        mg_coeff();
        mg_read_coeff(level);         	
		
	//mg_poisson_solver(level,fvar, pbc);

	mg_conjugate_gradient(level,0,200,0,0);
				
	mg_bcs_neu(level,0,pbc);
	mg_final(level,fvar,0); 
	
	write_restart(fvar,0); 
	write_to_file(node,fvar); 
}		

/***************************************************/

double get_wall_time()
{
    struct timeval time;

    if(gettimeofday(&time,NULL))
    {
        return 0; //Handle error
    }

    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

/***************************************************/
