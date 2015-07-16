#include <iostream> 
#include <cmath> 
#include <vector>
#include <fstream> 

using namespace std;

void write_order_field(vector<double>&, vector<double>&, vector<double>&, int, int); 
 

int main()
{
	int nx_h = 256, ny_h = 256;         //Spacings in the finest grid 
	int nx_2h = nx_h/2, ny_2h = ny_h/2; //Spacings in the 2h grid 
	int nx_4h = nx_h/4, ny_4h = ny_h/4; //Spacings in the 4h grid 

	int tot_hp = (nx_h+1)*(ny_h+1), tot_2hp = (nx_2h+1)*(ny_2h+1), tot_4hp = (nx_4h+1)*(ny_4h+1); 
	int str_h_x = nx_h + 1, str_2h_x = nx_2h + 1, str_4h_x = nx_4h + 1;

	int ind, p, q, ind_h, ind_2h, ind_4h; 

	vector<double> u_h(tot_hp), v_h(tot_hp), p_h(tot_hp); 
	vector<double> u_2h(tot_2hp), v_2h(tot_2hp), p_2h(tot_2hp); 
	vector<double> u_4h(tot_4hp), v_4h(tot_4hp), p_4h(tot_4hp);  
	vector<double> order_u(tot_4hp), order_v(tot_4hp), order_p(tot_4hp); 	

 	ifstream h_out, h2_out, h4_out; 
	ofstream order_out; 

	h_out.open("phi_h.dat");
	h2_out.open("phi_2h.dat");
	h4_out.open("phi_4h.dat");
	order_out.open("order_4h.dat");  

        /*Reading data from h-file*/
	
	for(int j=0;j<=ny_h;j++)
	{
		for(int i=0;i<=nx_h;i++)
		{
			ind = i + j*str_h_x; 
			h_out >> p >> q >> u_h[ind] >> v_h[ind] >> p_h[ind]; 
		}	       	
	}
	 
        /*Reading data from 2h-file*/
	
	for(int j=0;j<=ny_2h;j++)
        {
                for(int i=0;i<=nx_2h;i++)
                {
                        ind = i + j*str_2h_x; 
                        h2_out>>p>>q>>u_2h[ind]>>v_2h[ind]>>p_2h[ind]; 
                } 
  
        }
	
        /*Reading data from 4h-file*/ 

        for(int j=0;j<=ny_4h;j++)
        {
                for(int i=0;i<=nx_4h;i++)
                {
                        ind = i + j*str_4h_x; 
                        h4_out>>p>>q>>u_4h[ind]>>v_4h[ind]>>p_4h[ind]; 
                }  
        }

	/*******Computing order on the 4h grid***********/ 
	/**This is because it has all points in common**/ 

	for(int j=0;j<=ny_4h;j++)
	{
		for(int i=0;i<=nx_4h;i++)
		{
			ind_4h = i + j*str_4h_x; 
			ind_2h = 2*i + 2*j*str_2h_x; 
			ind_h =  4*i + 4*j*str_h_x; 
			
			//order_u[ind_4h] = log10( abs( (u_2h[ind_2h] - u_4h[ind_4h] )/(u_h[ind_h] - u_2h[ind_2h])  ) )/ log10(2.0);
			//order_v[ind_4h] = log10( abs( (v_2h[ind_2h] - v_4h[ind_4h] )/(v_h[ind_h] - v_2h[ind_2h]) )  )/ log10(2.0);
			order_p[ind_4h] = log10( abs( (p_2h[ind_2h] - p_4h[ind_4h] )/(p_h[ind_h] - p_2h[ind_2h]) ) )/ log10(2.0);

			order_out<<i<<"		"<<j<<"		"<<order_u[ind_4h]<<"		"<<order_v[ind_4h]<<"		"<<order_p[ind_4h]<<"\n"; 
			  			
		}
	}
	
	/***Writing the order field to a file***/

	write_order_field(order_u, order_v, order_p, nx_4h, ny_4h); 	    
        
	h_out.close();
	h2_out.close();
	h4_out.close(); 
	order_out.close(); 
}

/*********************************************************************/ 
/*Function writing the order-field to a file for visualization*/

void write_order_field(vector<double>& order_u, vector<double>& order_v, vector<double>& order_p, int nx_4h, int ny_4h)
{
	int str_4h_x = nx_4h + 1;
	int p, tot_p; 

	tot_p = (nx_4h+1)*(ny_4h+1); 
 
	ofstream order; 
	order.open("order_field.vtk");  	

	order<<"# vtk DataFile Version 2.0"<<"\n";
        order<<"2DGrid"<<"\n";
        order<<"ASCII"<<"\n";
        order<<"DATASET"<<" "<<"STRUCTURED_GRID"<<"\n";
        order<<"DIMENSIONS "<<nx_4h+1<<" "<<ny_4h+1<<" 1"<<"\n";
        order<<"POINTS"<<" "<<tot_p<<" float"<<"\n";

	for(int j=0;j<=ny_4h;j++)
        {
                for (int i=0;i<=nx_4h;i++)
                {
                        p = i + j*str_4h_x;
                        order<<i*(1./nx_4h)<<" "<<j*(1./ny_4h)<<" 0.0"<<"\n";
                }
        }


        order<<"POINT_DATA "<<tot_p<<"\n";
        order<<"SCALARS"<<" u_order"<<" double"<<"\n";
        order<<"LOOKUP_TABLE"<<" u_order"<<"\n";

        for (int j=0;j<=ny_4h;j++)
        {
                for (int i=0;i<=nx_4h;i++)
                {
                        p = i + j*str_4h_x;
                        order<<order_u[p]<<"\n";
                }
        }


        order<<"SCALARS"<<" v_order"<<" double"<<"\n";
        order<<"LOOKUP_TABLE"<<" v_order"<<"\n";

        for (int j=0;j<=ny_4h;j++)
        {
                for (int i=0;i<=nx_4h;i++)
                {
                        p = i + j*str_4h_x;
                        order<<order_v[p]<<"\n";
                }
        }
	

	order<<"SCALARS"<<" p_order"<<" double"<<"\n";
        order<<"LOOKUP_TABLE"<<" p_order"<<"\n";

        for (int j=0;j<=ny_4h;j++)
        {
                for (int i=0;i<=nx_4h;i++)
                {
                        p = i + j*str_4h_x;
                        order<<order_p[p]<<"\n";
                }
        }
}
