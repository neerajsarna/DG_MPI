#include "solve_system.h"

using namespace dealii;
void develop_system(system_data &system_matrices);

// now we specify the iniital and the boundary conditions
template<int dim>
class
ic_bc:public ic_bc_base<dim>
{
	public:
		ic_bc() {;};
		virtual double ic(const Point<dim> &p,const int &id);
		virtual void exact_solution(const Point<dim> &p,Vector<double> &value,const double &t);
};

int main(int argc, char *argv[])
{
      using namespace dealii;
      
      const unsigned int num_threads = atoi(argv[1]);
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, num_threads);

      const int dim = 2;
      const int poly_degree = 0;
      std::string foldername = "2D_advection_gaussian";

     system_data system_matrices;
     develop_system(system_matrices);
	  // develop mesh
     for(int i = 1 ; i < 30 ; i++ ) 
     {
      int repetitions = 10 * i;
      parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
      GridGenerator::subdivided_hyper_cube(triangulation,repetitions);


      ic_bc<dim> initial_boundary;

	  
	  Solve_System<dim> solve_system(system_matrices,
	  								 triangulation,
	  								 poly_degree,
	  								 &initial_boundary,
	  								 foldername);
 
	  solve_system.run_time_loop();
     }
}


void develop_system(system_data &system_matrices)
{

	// we can use the same sparsity pattern for all the test cases in this case
	system_matrices.Ax.resize(1,1);
	system_matrices.Ay.resize(1,1);
	system_matrices.P.resize(1,1);

	system_matrices.Ax.coeffRef(0,0) = 1;
	system_matrices.Ay.coeffRef(0,0) = 1;
	system_matrices.P.coeffRef(0,0) = 1;

	system_matrices.Ax.makeCompressed();
	system_matrices.Ay.makeCompressed();
	system_matrices.P.makeCompressed();
	
	//allocation of boundary matrices
	system_matrices.B.resize(4);
	system_matrices.penalty.resize(4);
	system_matrices.penalty_B.resize(4);


	for (unsigned int i = 0 ; i < 4 ; i ++)
	{
		system_matrices.B[i].resize(1,1);
		system_matrices.penalty[i].resize(1,1);
		system_matrices.penalty_B[i].resize(1,1);
	}
	//boundary at x = 1 
	unsigned int bc_id = 1;

	system_matrices.B[bc_id-1].coeffRef(0,0) = 0;
	system_matrices.penalty[bc_id-1].coeffRef(0,0) = 0;
	system_matrices.penalty_B[bc_id-1].coeffRef(0,0) = 0;

	// top boundary  y = 1
	bc_id = 2;

	system_matrices.B[bc_id-1].coeffRef(0,0) = 0;
	system_matrices.penalty[bc_id-1].coeffRef(0,0) = 0;
	system_matrices.penalty_B[bc_id-1].coeffRef(0,0) = 0;

	// left boundary x = 0
	bc_id = 3;

	system_matrices.B[bc_id-1].coeffRef(0,0) = 1;
	system_matrices.penalty[bc_id-1].coeffRef(0,0) = -1;
	system_matrices.penalty_B[bc_id-1].coeffRef(0,0) = -1;

	// boundary at y = 0, bottom boundary
	bc_id = 4;

	system_matrices.B[bc_id-1].coeffRef(0,0) = 1;
	system_matrices.penalty[bc_id-1].coeffRef(0,0) = -1;
	system_matrices.penalty_B[bc_id-1].coeffRef(0,0) = -1;

	// loop over all the boundaries
	for (unsigned int i = 0 ; i < 4 ; i ++)
	{
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

	const double density = exp(-pow((x-0.5),2)*100) * exp(-pow((y-0.5),2)*100);

	switch (id)
	{
		case 0:
			return(density);
		default:
			Assert(1 == 0 ,ExcMessage("should not have reached here"));
	}
	
	return 0;
}

template<int dim>
void 
ic_bc<dim>::exact_solution(const Point<dim> &p,Vector<double> &value,const double &t)
{
	const double x = p[0];
	const double y = p[1];

	value(0) = exp(-pow((x-t-0.5),2)*100) * exp(-pow((y-t-0.5),2)*100);

}
