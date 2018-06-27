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
						ic_bc_base<dim> *ic_bc,
						std::string &foldername,
						const double min_h);
		~Solve_System();
		
		MPI_Comm mpi_comm;

		DoFHandler<dim> dof_handler;
		FE_DGQ<dim> fe_basic;
		FESystem<dim> fe;

		IndexSet locally_relevant_dofs;
		IndexSet locally_owned_dofs;
		IndexSet locally_owned_cells; // helpful while computing error

		LA::MPI::Vector locally_relevant_solution;
		LA::MPI::Vector locally_owned_solution;
		LA::MPI::Vector error_per_cell;   


		ic_bc_base<dim> *initial_boundary;

		void prescribe_initial_conditions();
		const unsigned int n_eqn;
		system_data system_matrices;
		
		void run_time_loop();

		const double CFL = 1.0;
		double dt;
		double t_end = 1.0;
		double max_speed;

		// data structure for RK time stepping
		struct RK
		{
			Vector<double> weights;	// weights of the RK stages
			Vector<double> t_temp;  // intermidiate time steps
			Vector<double> dt_temp; // delta_t of different RK stages
			std::vector<LA::MPI::Vector> k_RK; // intermidate stages
			LA::MPI::Vector locally_owned_solution_temp; // temporary variable
		};

		RK RK_data;

		std::string output_foldername;

		void create_output();
		const unsigned int this_mpi_process;
		void compute_error();
		void create_IndexSet_triangulation();
		

		DeclException1 (ExcFileNotOpen, std::string,
                        << "Could not open: " << arg1);

		struct PerCellError
		{
			Vector<double> solution_value;
			Vector<double> exact_solution;
			double error_value;
			double cell_index;
			double volume;

			PerCellError(const unsigned int &n_eqn)
			:
			solution_value(n_eqn),
			exact_solution(n_eqn)
			{}
		};

		
		struct PerCellErrorScratch
		{};

		void compute_error_per_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
									 PerCellErrorScratch &scratch,
									 PerCellError &data);

		void copy_error_to_global(const PerCellError &data);
		ConditionalOStream pout;

    	//TimerOutput computing_timer;

};
