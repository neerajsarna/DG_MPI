
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
						Triangulation<dim> &triangulation,
					 	const int poly_degree,
						ic_bc_base<dim> *ic_bc,
						const unsigned int &maximum_neqn,
						const unsigned int &dim_problem);
		
		~Solve_System_SS_adaptive();

		const unsigned int dim_problem;
		MPI_Comm mpi_comm;
		DoFHandler<dim> dof_handler;
      	FE_DGQ<dim> fe_basic;
      	FESystem<dim> fe;

		IndexSet locally_relevant_dofs;
		IndexSet locally_owned_dofs;
		
		Vector<double> locally_relevant_solution;
		Vector<double> locally_owned_solution;
		Vector<double> locally_owned_residual;

		double min_h(Triangulation<dim> &triangulation);
		

		ic_bc_base<dim> *initial_boundary;

		void prescribe_initial_conditions();
		void distribute_dofs();
		void allocate_memory();
		std::vector<unsigned int> n_eqn;
		const unsigned int max_neqn;
		std::vector<system_data> system_matrices;

		const double CFL = 1.0;
		double dt;
		double t_end = 0.3;
		std::vector<double> max_speed;
		double current_max_speed;
		unsigned int current_max_index;


		void create_output(const std::string &filename);
		const unsigned int this_mpi_process;
		std::vector<double> compute_max_speed();		

		DeclException1 (ExcFileNotOpen, std::string,
                        << "Could not open: " << arg1);

		DeclException1 (ExcOrderingChanged, double,
                        << "distance between cells : " << arg1);

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
			unsigned int this_neqn;
		};

		struct PerCellAssembleScratch
		{
			    PerCellAssembleScratch(const FESystem<dim> &fe,
                	        		   const QGauss<dim-1> &   quadrature); // initialisation

    			PerCellAssembleScratch(const PerCellAssembleScratch &scratch);	// copy constructor

    			FEFaceValues<dim> fe_v_face;    			
		};



		void assemble_per_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                      PerCellAssembleScratch &scratch,
                                      PerCellAssemble &data,
                                      const std::vector<Vector<double>> &component_to_system,
                                      const double &t,
                                      const std::vector<Vector<double>> &g,
                                      const std::vector<Vector<double>> &force_vector);



		void integrate_face(Vector<double> &result,
			const typename DoFHandler<dim>::cell_iterator &neighbor,
			const std::vector<system_data> &system_matrices,
			const std::vector<Vector<double>> &component_to_system,
			const unsigned int &this_fe_index,
			const double &face_length,
			const double &volume,
			const double &nx,
			const double &ny,
			const Sparse_Matrix &An_cell,
			const std::vector<types::global_dof_index> &local_dof_indices);



		void assemble_to_global(const PerCellAssemble &data,
								const std::vector<Vector<double>> &component_to_system);


		ConditionalOStream pout;

		void compute_ic_per_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                  PerCellICScratch &scratch,PerCellIC &data,
                                  const Vector<double> &component);

		void copy_ic_to_global(const PerCellIC &data);


    	double residual_ss;


    	void run_time_loop(Triangulation<dim> &triangulation,
    					   const unsigned int &cycle,
    					   double &t,
    					   const std::vector<Vector<double>> &force_vector);

    	
    	void develop_neqn();

    	unsigned int current_max_fe_index();

    	int solve_steady_state(Triangulation<dim> &triangulation,double &t,
    							const std::vector<Vector<double>> &force_vector);

		Sparse_Matrix construct_An_effective(const Sparse_Matrix &An_cell,const Sparse_Matrix &An_neighbor);

		Vector<double> return_component();
		std::vector<Vector<double>> return_component_to_system();

		void compute_error();

		struct PerCellError
		{
			Vector<double> solution_value;
			Vector<double> exact_solution;
			double error_value;
			double cell_index;
			double volume;
		};

		double discretization_error;

		struct PerCellErrorScratch
		{};

		void compute_error_per_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
									 PerCellErrorScratch &scratch,
									 PerCellError &data);

		 void copy_error_to_global(const PerCellError &data);

};

// list of changes 
// 5. allocate fe index
