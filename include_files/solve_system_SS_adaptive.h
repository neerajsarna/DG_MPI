
#include "include_dealii.h"
#include "basic_data_structure.h"
#include "ic_bc_base.h"


// SS referes to the steady state
using namespace dealii;
template<int dim>
class 
Solve_System_SS_adaptive
{
	public:
		Solve_System_SS_adaptive(std::vector<system_data> &system_mat,
						parallel::distributed::Triangulation<dim> &triangulation,
					 	const int poly_degree,
						ic_bc_base<dim> *ic_bc,
						std::string &foldername);
		
		~Solve_System_SS_adaptive();

		MPI_Comm mpi_comm;
		hp::DoFHandler<dim> dof_handler;
      	FE_DGQ<dim> fe_basic;
      	hp::FECollection<dim> fe;

		IndexSet locally_relevant_dofs;
		IndexSet locally_owned_dofs;
		
		LA::MPI::Vector locally_relevant_solution;
		LA::MPI::Vector locally_owned_solution;
		LA::MPI::Vector locally_owned_residual;
		LA::MPI::Vector locally_relevant_solution_temp;
		LA::MPI::Vector error_per_cell;  
		LA::MPI::Vector edge_per_cell;

		double min_h(parallel::distributed::Triangulation<dim> &triangulation);
		

		ic_bc_base<dim> *initial_boundary;

		void prescribe_initial_conditions();
		void distribute_dofs();
		std::vector<unsigned int> n_eqn;
		std::vector<system_data> system_matrices;
	
		const double CFL = 1.0;
		double dt;
		double t_end = 0.3;
		std::vector<double> max_speed;
		double current_max_speed;

		std::string output_foldername;

		void create_output(const std::string &filename);
		const unsigned int this_mpi_process;
		void create_IndexSet_triangulation(IndexSet &locally_owned_cells);
		std::vector<double> compute_max_speed();		

		DeclException1 (ExcFileNotOpen, std::string,
                        << "Could not open: " << arg1);

		struct PerCellIC
		{
			Vector<double> local_contri;
			std::vector<types::global_dof_index> local_dof_indices;
			unsigned int dofs_per_cell;
		};

		struct PerCellICScratch
		{};

		struct PerCellAssemble
		{
			std::vector<types::global_dof_index> local_dof_indices;
			Vector<double> local_contri;
			unsigned int dofs_per_cell;
			unsigned int this_fe_index;
		};

		struct PerCellAssembleScratch
		{
			    PerCellAssembleScratch(const hp::FECollection<dim> &fe,
                	        		   const hp::QCollection<dim-1> &   quadrature); // initialisation

    			PerCellAssembleScratch(const PerCellAssembleScratch &scratch);	// copy constructor

    			hp::FEFaceValues<dim> fe_v_face;    			
		};

		void assemble_per_cell(const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
                                      PerCellAssembleScratch &scratch,
                                      PerCellAssemble &data,
                                      const std::vector<std::vector<Vector<double>>> &component_to_system,
                                      const double &t,
                                      const std::vector<std::vector<Vector<double>>> &g);

		void assemble_to_global(const PerCellAssemble &data,const std::vector<Vector<double>> &component);


		ConditionalOStream pout;

		void compute_ic_per_cell(const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
                                  PerCellICScratch &scratch,PerCellIC &data,
                                  const std::vector<Vector<double>> &component);

		void copy_ic_to_global(const PerCellIC &data);


    	double residual_ss;

    	TimerOutput computing_timer;

    	void run_time_loop(parallel::distributed::Triangulation<dim> &triangulation);

    	void develop_neqn();

    	void construct_fe_collection(); 

    	void construct_block_structure(std::vector<int> &block_structure,
                                       const std::vector<unsigned int> &n_eqn);

    	void allocate_fe_index(const unsigned int present_cycle);

    	unsigned int current_max_fe_index();

    	int solve_steady_state(parallel::distributed::Triangulation<dim> &triangulation,double &t);

		Sparse_Matrix construct_An_effective(const Sparse_Matrix &An_cell,const Sparse_Matrix &An_neighbor);

		std::vector<Vector<double>> return_component();
		std::vector<std::vector<Vector<double>>> return_component_to_system();

		// void refine_and_interpolate(parallel::distributed::Triangulation<dim> &triangulation);

		// void compute_error();

		// void compute_error_per_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
		// 							 PerCellErrorScratch &scratch,
		// 							 PerCellError &data);

		// struct PerCellError
		// {
		// 	Vector<double> solution_value;
		// 	Vector<double> exact_solution;
		// 	double error_value;
		// 	double cell_index;
		// 	double volume;

		// 	PerCellError(const unsigned int &n_eqn)
		// 	:
		// 	solution_value(n_eqn),
		// 	exact_solution(n_eqn)
		// 	{}
		// };

		// struct PerCellErrorScratch
		// {};

		// void copy_error_to_global(const PerCellError &data);

};

// list of changes 
// 5. allocate fe index
