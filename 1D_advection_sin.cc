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

                    // Following are the boundary ids:
                    // Left Wall = 0
                    // Bottom Wall = 1
                    // Right Wall = 2
                    // Top Wall = 3
                    
                    
                    // left edge
                    if (x_cord == left_edge)
                      cell->face(face)->set_boundary_id(0);

                    // this is the bottom wall
                    if (y_cord == left_edge)
                      cell->face(face)->set_boundary_id(1);

                   // right
                    if (x_cord == right_edge)
                      cell->face(face)->set_boundary_id(2);

                    // top edge, This is the second wall
                    if (y_cord == right_edge)
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
      std::string foldername = "1D_advection_sin";

     system_data system_matrices;
     develop_system(system_matrices);
     system_matrices.bc_inhomo_time = true;

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

      // for (unsigned int i = 29 ; i < 30 ; i = i + 2)
      // {
      repetitions[0] = atoi(argv[2]);
      //repetitions[0] = 15 * i;
      repetitions[1] = 2;

            //The diagonal of the rectangle is the line joining p1 and p2
      triangulation.clear();
      GridGenerator::subdivided_hyper_rectangle(triangulation,repetitions,p1,p2);
      set_square_bid(triangulation);

      ic_bc<dim> initial_boundary;
      
      ConditionalOStream pcout(std::cout,
          						(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0));
  	  TimerOutput   computing_timer(MPI_COMM_WORLD,
                    				pcout,
                    				TimerOutput::summary,
                    				TimerOutput::wall_times);	  

  	  computing_timer.enter_subsection("solving");
	  Solve_System<dim> solve_system(system_matrices,
	  								 triangulation,
	  								 poly_degree,
	  								 &initial_boundary,
	  								 foldername,
	  								 1.0/repetitions[0]);
 
	  solve_system.run_time_loop();

	  computing_timer.leave_subsection();
      pcout << std::endl;
	//}
     
}


void develop_system(system_data &system_matrices)
{

	// we can use the same sparsity pattern for all the test cases in this case
	system_matrices.Ax.resize(1,1);
	system_matrices.Ay.resize(1,1);
	system_matrices.P.resize(1,1);

	// zero in Ay since we have 1D advection
	system_matrices.Ax.coeffRef(0,0) = 1;
	system_matrices.Ay.coeffRef(0,0) = 0;
	system_matrices.P.coeffRef(0,0) = 0;

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

	// Following are the boundary ids:
                    // Left Wall = 0
                    // Bottom Wall = 1
                    // Right Wall = 2
                    // Top Wall = 3
	//left boundary at x = 0 
	unsigned int bc_id = 0;

	system_matrices.B[bc_id].coeffRef(0,0) = 1;
	system_matrices.penalty[bc_id].coeffRef(0,0) = -1;
	system_matrices.penalty_B[bc_id].coeffRef(0,0) = -1;

	// bottom boundary  y = 0
	bc_id = 1;

	system_matrices.B[bc_id].coeffRef(0,0) = 0;
	system_matrices.penalty[bc_id].coeffRef(0,0) = 0;
	system_matrices.penalty_B[bc_id].coeffRef(0,0) = 0;

	// right boundary x = 1
	bc_id = 2;

	system_matrices.B[bc_id].coeffRef(0,0) = 0;
	system_matrices.penalty[bc_id].coeffRef(0,0) = 0;
	system_matrices.penalty_B[bc_id].coeffRef(0,0) = 0;

	// top boundary at y = 1
	bc_id = 3;

	system_matrices.B[bc_id].coeffRef(0,0) = 0;
	system_matrices.penalty[bc_id].coeffRef(0,0) = 0;
	system_matrices.penalty_B[bc_id].coeffRef(0,0) = 0;

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

	//const double density = exp(-pow((x-0.5),2)*100);
	const double density = sin(M_PI * x);
	

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

	//value(0) = exp(-pow((x-t-0.5),2)*100);
	value(0) = sin(M_PI*(x-t));

}

template<int dim>
void 
ic_bc<dim>::bc_inhomo(const Sparse_Matrix &B,const unsigned int &bc_id,
								Vector<double> &value,const double &t)
{
	const int num_bc = B.rows();
	value.reinit(num_bc);

	value = 0;

	switch(bc_id)
	{
		case 0:
			// only inhomogeneity for the first element
			value(0) = -sin(M_PI * t);
			//value(0) = 0;
	}
}




