#include "include_dealii.h"
#include "basic_data_structure.h"
#include "ic_bc_base.h"


using namespace dealii;
template<int dim>
class 
Solve_System
{
	public:
		Solve_System( system_data &system_mat,
						parallel::distributed::Triangulation<dim> &triangulation,
						const int poly_degree,
						ic_bc_base<dim> *ic_bc);
		~Solve_System();
		
		MPI_Comm mpi_comm;

		DoFHandler<dim> dof_handler;
		FE_DGQ<dim> fe_basic;
		FESystem<dim> fe;

		IndexSet locally_relevant_dofs;
		IndexSet locally_owned_dofs;

		LA::MPI::Vector locally_relevant_solution;
		LA::MPI::Vector locally_owned_solution;


		ic_bc_base<dim> *initial_boundary;

		void prescribe_initial_conditions();
		const unsigned int n_eqn;
		system_data system_matrices;
		
		void run_time_loop();

		const double CFL = 0.5;
		double dt;
		double t_end = 0.3;


};