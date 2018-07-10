#include "solve_system_SS_adaptive.h"
#include "read_matrices.h"

using namespace dealii;
void develop_system(system_data &system_matrices,const int &M,const int &neqn_M,const int &nbc_M,const double &Kn);

// now we specify the iniital and the boundary conditions
template<int dim>
class
ic_bc:public ic_bc_base<dim>
{
	public:
		ic_bc() {;};
		virtual double ic(const Point<dim> &p,const int &id);
		virtual void exact_solution(const Point<dim> &p,Vector<double> &value,const double &t);
		virtual void bc_inhomo(const Sparse_Matrix &B,const unsigned int &bc_id,
								Vector<double> &value,const double &t);
};


void 
set_square_bid(parallel::distributed::Triangulation<2> &triangulation)
{
	typename parallel::distributed::Triangulation<2>::active_cell_iterator
											 cell = triangulation.begin_active(),
											 endc = triangulation.end();

    const double left_edge = 0;
    const double right_edge = 1;

    for(; cell != endc ; cell++)
    	if(cell->is_locally_owned())
    	{
                for (unsigned int face = 0 ; face < GeometryInfo<2>::faces_per_cell ; face++)
              
                  if (cell->face(face)->at_boundary())
                  { 
                    double x_cord = cell->face(face)->center()(0);
                    double y_cord = cell->face(face)->center()(1);

                    // this boundary ID is coherent (-1) with the matlab code

                    // Following are the boundary ids:
                    // Right Wall = 0
                    // Top wall = 1
                    // left wall = 2
                    // bottom wall = 3
                    
                   // right
                    if (x_cord == right_edge)
                      cell->face(face)->set_boundary_id(0);

                    // top edge
                    if (y_cord == right_edge)
                      cell->face(face)->set_boundary_id(1);

                    // left edge
                    if (x_cord == left_edge)
                      cell->face(face)->set_boundary_id(2);

                    // bottom 
                    if (y_cord == left_edge)
                      cell->face(face)->set_boundary_id(3);
                   }
    	}
}

void count_bc_id(parallel::distributed::Triangulation<2> &triangulation)
{
	  typename parallel::distributed::Triangulation<2>::active_cell_iterator cell = triangulation.begin_active(),
	  													 endc = triangulation.end();


	  
	  int id0 = 0,id1 = 0, id2 = 0, id3 = 0; 

	  for(; cell != endc ; cell++)   
	  	if (cell->is_locally_owned())
	  {

	  	for(unsigned int face = 0 ; face < GeometryInfo<2>::faces_per_cell ; face++)
	  		if (cell->face(face)->at_boundary())
	  			switch (cell->face(face)->boundary_id())
	  			{
	  				case 0:
	  				{
	  					id0++;
	  					break;
	  				}

	  				case 1:
	  				{
	  					id1++;
	  					break;
	  				}

	  				case 2:
	  				{
	  					id2++;
	  					break;
	  				}

	  				case 3:
	  				{
	  					id3++;
	  					break;
	  				}
	  			}
	  		}
	  		
	  

	  	  	std::cout << "id0 " << id0 
	  			<< " id1 " << id1 
	  			<< " id2 " << id2 
	  			<< " id3 " << id3 << std::endl;

	  		fflush(stdout);


}


int main(int argc, char *argv[])
{
      using namespace dealii;
      
      const unsigned int num_threads = atoi(argv[1]);
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, num_threads);

      const int dim = 2;
      const int poly_degree = 0;

      // the output folder name
      std::string foldername = "2x3v_moments_HC";

     
     // store the number of equations for a given M (starts at M = 3)
     const int neqn_M[18] = {13,22,34,50,70,95,125,161,203,252,308,372,444,525,615,715,825,946};
     const int nbc_M[18] = {5,8,14,20,30,40,55,70,91,112,140,168,204,240,285,330,385,440};
     const double Kn = 0.1;

     std::vector<int> M(1);
     M[0] = 3;
     std::vector<system_data> system_matrices(M.size());

     Assert(*std::max_element(M.begin(),M.end())<=20,ExcNotImplemented());

     for(unsigned int i = 0 ; i < M.size(); i++)
     {
     	develop_system(system_matrices[i],M[i],neqn_M[M[i]-3],nbc_M[M[i]-3],Kn);
     	system_matrices[i].bc_inhomo_time = true;
     }

      // create a rectangular mesh 
      parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
      Point<dim> p1;
      Point<dim> p2;
      std::vector<unsigned int> repetitions(dim);

      const double left_edge = 0;
      const double right_edge = 1;

      // corners of the diagonal
      p1(0) = left_edge;
      p1(1) = left_edge;

      p2(0) = right_edge;
      p2(1) = right_edge;

     
      repetitions[0] = atoi(argv[2]);
      repetitions[1] = 1;

      //The diagonal of the rectangle is the line joining p1 and p2
      GridGenerator::subdivided_hyper_rectangle(triangulation,repetitions,p1,p2);
      set_square_bid(triangulation);

     
      triangulation.signals.post_refinement.connect(std_cxx11::bind (&set_square_bid,
                      								std_cxx11::ref(triangulation)));

      ic_bc<dim> initial_boundary;	

	  Solve_System_SS_adaptive<dim> solve_system(system_matrices,
	  								 triangulation,
	  								 poly_degree,
	  								 &initial_boundary,
	  								 foldername);

	 solve_system.run_time_loop(triangulation);

	 // develop the filename for output
	 // std::string filename = "2x3v_moments_HC/M" + std::to_string(M)
	 // 						 + "/result" + std::to_string(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
	 // 						 + "_Kn_" + "0p1" + ".txt";


	 // solve_system.create_output(filename);
 
     
}


void develop_system(system_data &system_matrices,const int &M,const int &neqn_M,
												 const int &nbc_M,const double &Kn)
{
	std::cout << "developing systems: " << std::endl;
	// we first initialise all the matrices
	system_matrices.Ax.resize(neqn_M,neqn_M);
	system_matrices.Ay.resize(neqn_M,neqn_M);
	system_matrices.P.resize(neqn_M,neqn_M);

	system_matrices.B.resize(4);
	system_matrices.penalty_B.resize(4);
	system_matrices.penalty.resize(4);

	for(unsigned int id = 0 ; id < 4 ; id ++)
	{
		system_matrices.B[id].resize(nbc_M,neqn_M);
		system_matrices.penalty_B[id].resize(neqn_M,neqn_M);
		system_matrices.penalty[id].resize(neqn_M,nbc_M);
	}

	std::vector<triplet> Row_Col_Value;
	std::string filename = "3v_Moments/Ax/Ax" + std::to_string(M) + ".txt";
	build_triplet(Row_Col_Value,filename);
	build_matrix_from_triplet(system_matrices.Ax,Row_Col_Value);

	filename = "3v_Moments/P/P" + std::to_string(M) + ".txt";
	build_triplet(Row_Col_Value,filename);
	build_matrix_from_triplet(system_matrices.P,Row_Col_Value);
	system_matrices.P = system_matrices.P/Kn;
	system_matrices.P.makeCompressed();

	filename = "3v_Moments/Bwall/Bwall" + std::to_string(M) + ".txt";
	build_triplet(Row_Col_Value,filename);
	build_matrix_from_triplet(system_matrices.B[0],Row_Col_Value);

	filename = "3v_Moments/Bwall/penalty_wall" + std::to_string(M) + ".txt";
	build_triplet(Row_Col_Value,filename);
	build_matrix_from_triplet(system_matrices.penalty[0],Row_Col_Value);	
	
	std::vector<Sparse_Matrix> rotator(4);

	for(unsigned int id = 0 ; id < 4 ; id++ )
	{
		rotator[id].resize(neqn_M,neqn_M);
		filename = "3v_Moments/Rotator/rotator" + std::to_string(M)
					 + '_' + std::to_string(id+1) + ".txt";
					 
		build_triplet(Row_Col_Value,filename);
		build_matrix_from_triplet(rotator[id],Row_Col_Value);		
	}

	system_matrices.Ax.makeCompressed();

	for(unsigned int i = 0 ; i < 4 ; i++)
	{
		system_matrices.B[i] = system_matrices.B[0] * rotator[i];
		system_matrices.penalty[i] = rotator[i].transpose() * system_matrices.penalty[0];

		if (i == 1 || i == 3) // only changes along the x-direction
			system_matrices.penalty[i] = rotator[i].transpose() * system_matrices.penalty[0] * 0;

		system_matrices.penalty_B[i] = system_matrices.penalty[i] * system_matrices.B[i];

		system_matrices.B[i].makeCompressed();
		system_matrices.penalty[i].makeCompressed();
		system_matrices.penalty_B[i].makeCompressed();
	}




}


template<int dim>
double 
ic_bc<dim>::ic(const Point<dim> &p,const int &id)
{
	const double x = p[0];
	const double y = p[1];

	const double density = 0;	
	// const double density = exp(-pow((x-0.5),2)*100);	

	switch (id)
	{
		case 0:
		{
			return(density);
			break;
		}
		default:
		{
			return(0);
			break;
		}
	}
	
	return 0;
}

template<int dim>
void 
ic_bc<dim>::exact_solution(const Point<dim> &p,Vector<double> &value,const double &t)
{
	const double x = p[0];
	const double y = p[1];

	value(0) = 0;

}

template<int dim>
void 
ic_bc<dim>::bc_inhomo(const Sparse_Matrix &B,const unsigned int &bc_id,
								Vector<double> &value,const double &t)
{

	const int num_bc = B.rows();
	double thetaW;
	value.reinit(num_bc);
	value = 0;

	if (t <= 1)
        thetaW = exp(-1/(1-pow((t-1),2))) * exp(1);
    else
        thetaW = 1;

    
	switch (bc_id)
	{
		case 1:
		case 3:
		{
			value = 0;
			break;
		}
		case 0:
		{

			double bc_normal = 1.0;
			for (unsigned int m = 0 ; m < B.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(B,m); n ; ++n)
                    	if (n.col() == 3 || n.col() == 5 || n.col() == 6)
                    		value(n.row()) += bc_normal * thetaW * n.value()/sqrt(2.0);
                    	
			break;
		}

		case 2:
		{
			double bc_normal = -1.0;
			for (unsigned int m = 0 ; m < B.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(B,m); n ; ++n)
                    	if (n.col() == 3 || n.col() == 5 || n.col() == 6)
                    		value(n.row()) += bc_normal * thetaW * n.value()/sqrt(2.0);
                    	
			break;
		}

		default:
		{
			Assert(1 == 0, ExcMessage("should not have reached"));
			break;
		}
	}

}




