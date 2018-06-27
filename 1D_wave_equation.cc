#include "solve_system.h"
#include "read_matrices.h"

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


int main(int argc, char *argv[])
{
      using namespace dealii;
      
      const unsigned int num_threads = atoi(argv[1]);
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, num_threads);

      const int dim = 2;
      const int poly_degree = 0;

      // the output folder name
      std::string foldername = "2x3v_moments_gaussian_inflow";

     system_data system_matrices;
     develop_system(system_matrices);
     system_matrices.bc_inhomo_time = false;

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
      repetitions[1] = 3;

      //The diagonal of the rectangle is the line joining p1 and p2
      GridGenerator::subdivided_hyper_rectangle(triangulation,repetitions,p1,p2);
      set_square_bid(triangulation);

      ic_bc<dim> initial_boundary;	  
	  Solve_System<dim> solve_system(system_matrices,
	  								 triangulation,
	  								 poly_degree,
	  								 &initial_boundary,
	  								 foldername,
	  								 1.0/repetitions[0]);
 
	  solve_system.run_time_loop();
     
}


void develop_system(system_data &system_matrices)
{
	system_matrices.Ax.resize(2,2);
	system_matrices.Ay.resize(2,2);
	system_matrices.P.resize(2,2);

	system_matrices.Ax.coeffRef(0,1) = 1;
	system_matrices.Ax.coeffRef(1,0) = 1;
	system_matrices.Ax.makeCompressed();

	system_matrices.B.resize(4);
	system_matrices.penalty.resize(4);
	system_matrices.penalty_B.resize(4);

	for (unsigned int id = 0 ; id < 4 ; id ++)
	{
		system_matrices.B[id].resize(1,2);
		system_matrices.penalty[id].resize(2,1);
		system_matrices.penalty_B[id].resize(2,2);
	}

		system_matrices.B[0].coeffRef(0,0) = 1;
		system_matrices.B[0].coeffRef(0,1) = -1;
		system_matrices.penalty[0].coeffRef(1,0) = 1;
		system_matrices.penalty_B[0] = system_matrices.penalty[0] * system_matrices.B[0];


		system_matrices.B[2].coeffRef(0,0) = 1;
		system_matrices.B[2].coeffRef(0,1) = 1;
		system_matrices.penalty[2].coeffRef(1,0) = -1;
		system_matrices.penalty_B[2] = system_matrices.penalty[2] * system_matrices.B[2];

		// std::cout << "Ax " << std::endl;
		// std::cout << system_matrices.Ax << std::endl;
 
 	// 	std::cout << "Ay " << std::endl;
		// std::cout << system_matrices.Ay << std::endl;

		// std::cout << "B0 " << std::endl;
		// std::cout << system_matrices.B[0] << std::endl;

		// std::cout << "B1 " << std::endl;
		// std::cout << system_matrices.B[1] << std::endl;

		// std::cout << "B2 " << std::endl;
		// std::cout << system_matrices.B[2] << std::endl;

		// std::cout << "B3 " << std::endl;
		// std::cout << system_matrices.B[3] << std::endl;

		// std::cout << "penalty0 " << std::endl;
		// std::cout << system_matrices.penalty[0] << std::endl;

		// std::cout << "penalty1 " << std::endl;
		// std::cout << system_matrices.penalty[1] << std::endl;

		// std::cout << "penalty2 " << std::endl;
		// std::cout << system_matrices.penalty[2] << std::endl;

		// std::cout << "penalty3 " << std::endl;
		// std::cout << system_matrices.penalty[3] << std::endl;


		// std::cout << "penalty_B0 " << std::endl;
		// std::cout << system_matrices.penalty_B[0] << std::endl;

		// std::cout << "penalty_B1 " << std::endl;
		// std::cout << system_matrices.penalty_B[1] << std::endl;

		// std::cout << "penalty_B2 " << std::endl;
		// std::cout << system_matrices.penalty_B[2] << std::endl;

		// std::cout << "penalty_B3 " << std::endl;
		// std::cout << system_matrices.penalty_B[3] << std::endl;		
}


template<int dim>
double 
ic_bc<dim>::ic(const Point<dim> &p,const int &id)
{
	const double x = p[0];
	const double y = p[1];

	//const double density = exp(-pow((x-0.5),2)*100) * exp(-pow((y-0.5),2)*100);	
	const double density = exp(-pow((x-0.5),2)*100);	

	switch (id)
	{
		case 0:
			return(density);
		case 1:
			return(0);
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
	value.reinit(num_bc);

	value = 0;

}




