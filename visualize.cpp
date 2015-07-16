#include "headers.h"
#include "init_2.h"
#include "declarations.h"

///////////////////Function overloaded of Writing to a file////////////////////
//Writing data at frequency specified by file_freq in the "init" file

void write_to_file(vector<vertex>& node, vector<fval>& fvar, int t)
{
	ofstream out_put, vel_put0, vel_put1; //0->x-velocity along y-direction at geometric center and 1->y-velocity along x-direction at geometric center  
	ofstream vort1, vort2, flux; 
	int p,p1,p2;
	double udy=0.0, vdx=0.0; 

        tdma1y(fvar,0);
	tdma1x(fvar,1);  

	string file("data_2D");
	string file0("x-vel");
	string file1("y-vel");
        string file2("vort-mov");  
  	string file3("vort-cen");
	string file4("mass_flux_gc");  

	stringstream tag; 

	tag<<t;

	file = file + "_" + tag.str() + ".vtk";   
	file0 = file0 + "_" + tag.str() + ".dat"; 
	file1 = file1 + "_" + tag.str() + ".dat";
	file2 = file2 + "_" + tag.str() + ".dat";
	file3 = file3 + "_" + tag.str() + ".dat";
	file4 = file4 + "_" + tag.str() + ".dat";

	out_put.open(file.c_str());

	vel_put0.open(file0.c_str());
	vel_put1.open(file1.c_str());
        vort1.open(file2.c_str()); 
        vort2.open(file3.c_str());
	flux.open(file4.c_str());  

	out_put<<"# vtk DataFile Version 2.0"<<"\n"; 
	out_put<<"2DGrid"<<"\n";
	out_put<<"ASCII"<<"\n";
	out_put<<"DATASET"<<" "<<"STRUCTURED_GRID"<<"\n";
	out_put<<"DIMENSIONS "<<nx+1<<" "<<ny+1<<" 1"<<"\n";
	out_put<<"POINTS"<<" "<<tot_p<<" float"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<node[p].x[0]<<" "<<node[p].x[1]<<" 0.0"<<"\n"; 
		}
	}

	out_put<<"POINT_DATA "<<tot_p<<"\n";
	out_put<<"SCALARS"<<" Pressure"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" Pressure_lp"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].u[2]<<"\n"; 
		}
	}
	out_put<<"SCALARS"<<" divergence"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" div_lp"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].div<<"\n"; 
		}
	}
	out_put<<"SCALARS"<<" residual"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" res_lp"<<"\n";

	for(int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].res<<"\n"; 
		}
	}
	
	out_put<<"SCALARS"<<" F"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" F_lp"<<"\n";

	for(int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].F<<"\n"; 
		}
	}

	out_put<<"SCALARS"<<" Vorticity"<<" double"<<"\n";
        out_put<<"LOOKUP_TABLE"<<" Vorticity_lp"<<"\n";

        for(int j=0;j<=ny;j++)
        {
                for (int i=0;i<=nx;i++)
                {
                        p = i + j*str_x;
                        out_put<<fvar[p].uy[0]-fvar[p].ux[1]<<"\n";
                }
        }


	out_put<<"VECTORS"<<" Velocity"<<" double"<<"\n"; 

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].u[0]<<" "<<fvar[p].u[1]<<" "<<"0.0"<<"\n";
		}
	}

	//Writing the x-velocity in the y-direction at the geometric center of the domain 
        
	for(int j=0;j<=ny;j++)
	{
		vel_put0<<node[nx/2 + j*str_x].x[1]<<"	"<<fvar[nx/2 + j*str_x].u[0]<<"\n"; 
	}

	//Writing the y-velocity in the x-direction at the geometric center of the domain 

	for(int i=0;i<=nx;i++)
	{
		vel_put1<<node[i+(ny/2)*str_x].x[0]<<"	"<<fvar[i+(ny/2)*str_x].u[1]<<"\n"; 

	}

	//Writing vorticity along the moving plate

	for(int i=1;i<nx;i++)
	{
		p = i + ny*str_x;                
		vort1<<node[p].x[0]<<"	"<<fvar[p].uy[0] - fvar[p].ux[1]<<"\n";
	}

        //Writing vorticity at the geometric center along a vertical line

	for(int j=0;j<=ny;j++)
        {
		p = nx/2 + j*str_x; 
                vort2<<node[p].x[1]<<"  "<<fvar[p].uy[0] - fvar[p].ux[1]<<"\n";
        }

	//Writing mass fluxes 

	for(int j=1;j<=ny;j++)
	{
		p1 = (nx/2) + (j-1)*str_x; 
		p2 = (nx/2) + j*str_x; 

		udy = udy + (fvar[p1].u[0] + fvar[p2].u[0])*0.5*dy; 
	}

	for(int i=1;i<=nx;i++)
	{
		p = i + (ny/2)*str_x; 		
		
		vdx = vdx + (fvar[p-1].u[1] + fvar[p].u[1])*0.5*dx; 
	}	

	flux<<"udy= "<<udy<<"\n"; 
	flux<<"vdx= "<<vdx<<"\n"; 


	out_put.close(); 
	vel_put0.close();
	vel_put1.close();
	vort1.close();
	vort2.close(); 
	flux.close(); 
}

///////////////////Writing to a file////////////////////

void write_to_file(vector<vertex>& node, vector<fval>& fvar)
{

	ofstream out_put;  
	int p; 
	
	vec error_diff = zeros<vec>(tot_p); 

	out_put.open("doubly_periodic.vtk");

	out_put<<"# vtk DataFile Version 2.0"<<"\n"; 
	out_put<<"2DGrid"<<"\n";
	out_put<<"ASCII"<<"\n";
	out_put<<"DATASET"<<" "<<"STRUCTURED_GRID"<<"\n";
	out_put<<"DIMENSIONS "<<nx+1<<" "<<ny+1<<" 1"<<"\n";
	out_put<<"POINTS"<<" "<<tot_p<<" float"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<node[p].x[0]<<" "<<node[p].x[1]<<" 0.0"<<"\n"; 
		}
	}

	out_put<<"POINT_DATA "<<tot_p<<"\n";
	out_put<<"SCALARS"<<" Pressure"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" Pressure_lp"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			//out_put<<i<<"	"<<j<<"	   "<<fvar[p].u[2]<<"\n"; 
			out_put<<fvar[p].u[2]<<"\n"; 
		}
	}
	
	out_put<<"SCALARS"<<" rhs"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" rhs"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].rhs<<"\n"; 
		}
	}
	
	out_put<<"SCALARS"<<" exact_sol"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" exact_sol"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<sin(2.*M_PI*node[p].x[0])*sin(2.*M_PI*node[p].x[1]) <<"\n"; 
		}
	}
	
	out_put<<"SCALARS"<<" error_diff"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" res_lp"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			
			out_put<<fvar[p].u[2] - ( sin(2.*M_PI*node[p].x[0])*sin(2.*M_PI*node[p].x[1]) ) + sin(2.*M_PI*(nx-1)*dx)*sin(2.*M_PI*(ny-1)*dy) <<"\n";					
		}
	}
	
	/*cout<<"The error norm is "<<norm(error_diff,2)<<"\n";*/ 
	
	out_put<<"SCALARS"<<" residual"<<" double"<<"\n"; 
	out_put<<"LOOKUP_TABLE"<<" res_lp"<<"\n";

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 

			out_put<<fvar[p].res<<"\n"; 
		}
	}

	out_put<<"VECTORS"<<" Velocity"<<" double"<<"\n"; 

	for (int j=0;j<=ny;j++)
	{ 
		for (int i=0;i<=nx;i++)
		{
			p = i + j*str_x; 
			out_put<<fvar[p].u[0]<<" "<<fvar[p].u[1]<<" "<<"0.0"<<"\n";
		}
	}

	out_put.close(); 

}
