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
						std::string &foldername);
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

		void assemble_to_global(const PerCellAssemble &data);

		void compute_error_per_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
									 PerCellErrorScratch &scratch,
									 PerCellError &data);

		void copy_error_to_global(const PerCellError &data);
		ConditionalOStream pout;

		void compute_ic_per_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                  PerCellICScratch &scratch,PerCellIC &data,const Vector<double> &component);

		void copy_ic_to_global(const PerCellIC &data);

		double min_h(const parallel::distributed::Triangulation<dim> &triangulation);


    	//TimerOutput computing_timer;

};
