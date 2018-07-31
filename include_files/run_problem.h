#include "solve_system_SS_adaptive.h"

template<int dim>
class
run_problem
{
	public:
		run_problem(std::vector<system_data> &system_mat,	  // system data
					std::vector<system_data> &system_mat_error,	  // system data
				    std::vector<system_data> &system_mat_adj, // adjoint data
					Triangulation<dim> &triangulation, // triangulation
					const int poly_degree,
					ic_bc_base<dim> *ic_bc_primal,
					ic_bc_base<dim> *ic_bc_adjoint,
					const std::string &foldername);


		void compute_error_velocity(const hp::DoFHandler<dim> &dof_handler_primal,
						   const std::vector<system_data> &system_matrices,
						   const std::vector<unsigned int> &n_eqn_primal,
							const std::vector<unsigned int> &n_eqn_adjoint,
						   ic_bc_base<dim> *ic_bc,
						   const std::vector<Vector<double>> &adjoint_solution,
						   const std::vector<Vector<double>> &primal_solution,
						   const std::vector<Point<dim>> &cell_index_primal,
                           const std::vector<Point<dim>> &cell_index_adjoint);

 		void compute_error_h(const Triangulation<dim> &triangulation,
                             const hp::DoFHandler<dim> &dof_handler_adjoint,
                             const Vector<double> &value_adjoint,
                             const std::vector<Vector<double>> &primal_solution,
                             const std::vector<system_data> &system_matrices,
                             const std::vector<unsigned int> &n_eqn_adjoint,
                             ic_bc_base<dim> *ic_bc);

		// x^T A y
		double xAy(const Vector<double> &x,
							 const Sparse_Matrix &A,
							 const Vector<double> &y);

		 void error_face(const typename hp::DoFHandler<dim>::cell_iterator &neighbor,
                              const unsigned int &this_fe_index,
                              const std::vector<system_data> &system_matrices,
                              const Vector<double> &adjoint_value,
                              const Vector<double> &solution_value,
                              const Vector<double> &neighbor_value,
                              const double &face_length,
                              const double &nx,
                              const double &ny,
                              double &result);

		 void error_face_h(const typename DoFHandler<dim>::cell_iterator &neighbor,  // iterator of the neighbor
                                const std::vector<Vector<double>> &primal_solution,       // primal solution
                                const Vector<double> &primal_value,                       // value of the current cell
                                const unsigned int &this_fe_index,                        // fe index of the present cell
                                const std::vector<system_data> &system_matrices,          // matrices
                                const std::vector<types::global_dof_index> &local_dof_indices, // local dof indices for solution construction
                                const std::vector<Vector<double>> &component_to_system,         // component to system index of the adjoint
                                const Vector<double> &adjoint_solution,                         // the adjoint solution
                                const FEValuesBase<dim> &fe_v_face,                       
                                const Vector<int> &cell_fe_index,
                                double &result);

		DeclException1 (ExcOrderingChanged, double,
                        << "distance between cells : " << arg1);

		Vector<double> error_per_cell;

		double t;		// time of the computation


		struct PerCellError
		{
			
			double local_contri;
			unsigned int index;

		};

		struct PerCellErrorScratch
		{
			    PerCellErrorScratch(const hp::FECollection<dim> &fe,
                	        		const hp::QCollection<dim-1> &   quadrature); // initialisation

    			PerCellErrorScratch(const PerCellErrorScratch &scratch);	// copy constructor

    			hp::FEFaceValues<dim> fe_v_face;    			
		};

		struct PerCellErrorScratch_h
		{
			    PerCellErrorScratch_h(const FiniteElement<dim> &fe,
			    					            const QGauss<dim> &   quadrature_int,
                	        		  const QGauss<dim-1> &   quadrature); // initialisation

    			PerCellErrorScratch_h(const PerCellErrorScratch_h &scratch);	// copy constructor

    			FEFaceValues<dim> fe_v_face;    			
    			FESubfaceValues<dim> fe_v_subface;    			
    			FEValues<dim> fe_v;    			
		};

		void compute_error_per_cell_velocity(const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
                                      	PerCellErrorScratch &scratch,
                                      	PerCellError &data,
                                        const std::vector<std::vector<Vector<double>>> &g,
                                        const std::vector<system_data> &system_matrices,
                                        const std::vector<unsigned int> &n_eqn_primal,
                                        const std::vector<unsigned int> &n_eqn_adjoint,
                                        const std::vector<Vector<double>> &adjoint_solution,
						   				const std::vector<Vector<double>> &primal_solution,
						   				const std::vector<Point<dim>> &cell_index_primal,
                                        const std::vector<Point<dim>> &cell_index_adjoint);

		void compute_error_per_cell_h(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                        PerCellErrorScratch_h &scratch,
                                        PerCellError &data,
                                        const std::vector<std::vector<Vector<double>>> &g,
                                        const std::vector<system_data> &system_matrices,
                                        const Vector<double> &adjoint_solution,
                                        const std::vector<Vector<double>> &primal_solution,
                                        const std::vector<Vector<double>> &component_to_system_adjoint,
                                        const Vector<int> &cell_fe_index);


		void assemble_to_global(const PerCellError &data);

		Sparse_Matrix construct_An_effective(const Sparse_Matrix &An_cell,const Sparse_Matrix &An_neighbor);

		void write_error(const std::string &filename,const Triangulation<dim> &triangulation);
    void write_grid(const std::string &filename,const Triangulation<dim> &triangulation);
		void write_fe_index(const std::string &filename,
                            const hp::DoFHandler<dim> &dof_handler,
                            const std::vector<unsigned int> &n_eqn);


		void print_fe_index(const hp::DoFHandler<dim> &dof_handler,const std::vector<unsigned int> &n_eqn);
		
    	void construct_fe_collection(const FiniteElement<dim> &fe_basic,
                                                       const std::vector<unsigned int> &n_eqn,
                                                       hp::FECollection<dim> &fe); 


    	void get_value_at_quad(const std::vector<types::global_dof_index> &local_dof_indices,
                                    const std::vector<Vector<double>> &component_to_system,
                                    Vector<double> &adjoint_value,  // value to be filled
                                    const unsigned int &dofs_per_component,
                                    const Vector<double> &adjoint_solution,
                                    const FEValuesBase<dim> &fe_v,
                                    const unsigned int &q);

    void construct_block_structure(std::vector<int> &block_structure,
                                       const std::vector<unsigned int> &n_eqn);

    void   extrapolate(const hp::DoFHandler<dim> &dofIn,
                        const Vector<double> &InVec,
                        const DoFHandler<dim>  &dofOut,
                        Vector<double> &OutVec,
                        const Triangulation<dim> &triangulation,
                        Vector<int> &cell_fe_index);


    void interpolate_to_highest_fe_index(const hp::DoFHandler<dim> &dofIn,
                                         const Vector<double> &InVec,
                                         const DoFHandler<dim>  &dofOut,
                                         Vector<double> &OutVec,
                                         Vector<int> &cell_fe_index);

};
