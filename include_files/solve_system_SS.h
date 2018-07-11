#include "include_dealii.h"
#include "basic_data_structure.h"
#include "ic_bc_base.h"


// SS referes to the steady state
using namespace dealii;
template<int dim>
class 
Solve_System_SS
{
	public:
		Solve_System_SS(system_data &system_mat,
						parallel::distributed::Triangulation<dim> &triangulation,
					 	const int poly_degree,
						ic_bc_base<dim> *ic_bc,
						std::string &foldername);
		
		~Solve_System_SS();

		MPI_Comm mpi_comm;
		DoFHandler<dim> dof_handler;
      	FE_DGQ<dim> fe_basic;
      	FESystem<dim> fe;

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
		unsigned int n_eqn;
		system_data system_matrices;
	
		const double CFL = 0.1;
		double dt;
		double t_end = 0.3;
		double max_speed;

		std::string output_foldername;

		void create_output(const std::string &filename);
		const unsigned int this_mpi_process;
		void compute_error();
		void create_IndexSet_triangulation(IndexSet &locally_owned_cells);
		double compute_max_speed();		

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
		};

		struct PerCellAssembleScratch
		{
			    PerCellAssembleScratch(const FiniteElement<dim> &fe,
                	        		   const Quadrature<dim-1> &   quadrature); // initialisation

    			PerCellAssembleScratch(const PerCellAssembleScratch &scratch);	// copy constructor

    			FEFaceValues<dim> fe_v_face;    			
		};

		void assemble_per_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                      PerCellAssembleScratch &scratch,
                                      PerCellAssemble &data,
                                      const Vector<double> &component,
                                      const std::vector<Vector<double>> &component_to_system,
                                      const double &t,
                                      const std::vector<Vector<double>> &g);

		void assemble_to_global(const PerCellAssemble &data,const Vector<double> &component);

		void compute_error_per_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
									 PerCellErrorScratch &scratch,
									 PerCellError &data);

		void copy_error_to_global(const PerCellError &data);
		ConditionalOStream pout;

		void compute_ic_per_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                  PerCellICScratch &scratch,PerCellIC &data,const Vector<double> &component);

		void copy_ic_to_global(const PerCellIC &data);


    	double residual_ss;

    	TimerOutput computing_timer;

    	void run_time_loop(parallel::distributed::Triangulation<dim> &triangulation);

};
