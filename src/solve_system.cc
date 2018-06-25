
#include "solve_system.h"


template<int dim>
Solve_System<dim>::Solve_System(system_data &system_mat,
								parallel::distributed::Triangulation<dim> &triangulation,
								const int poly_degree,
								ic_bc_base<dim> *ic_bc)
:
mpi_comm(MPI_COMM_WORLD),
dof_handler(triangulation),
fe_basic(poly_degree),
fe(fe_basic,system_mat.Ax.rows()),
initial_boundary(ic_bc),
n_eqn(fe.n_components()),
system_matrices(system_mat)
{
	
	  // we store data of the system ( flux matrices and boundary matrices)

	  dof_handler.distribute_dofs(fe);
	  
	  locally_owned_dofs = dof_handler.locally_owned_dofs ();
      DoFTools::extract_locally_relevant_dofs (dof_handler,
                          locally_relevant_dofs) ;
    

      locally_owned_solution.reinit(locally_owned_dofs,mpi_comm);
      locally_relevant_solution.reinit(locally_owned_dofs,locally_relevant_dofs,mpi_comm);

      prescribe_initial_conditions();
      locally_relevant_solution = locally_owned_solution;

      const double min_h = GridTools::minimal_cell_diameter(triangulation);

      // we need to fix this in case of moments
      const double max_speed = 1;

      // an approximation to delta_t
      dt = CFL * min_h/max_speed;

}


// prescribe the initial conditions
template<int dim>
void
Solve_System<dim>::prescribe_initial_conditions()
{

	typename DoFHandler<dim>::active_cell_iterator cell =
                                                   dof_handler.begin_active(),
                                                 endc = dof_handler.end();

	const unsigned int dofs_per_cell = fe.dofs_per_cell;
    Vector<double> component(dofs_per_cell);
    Vector<double> value(dofs_per_cell);


    for (unsigned int i = 0 ; i < dofs_per_cell ; i++)
      component[i] = fe.system_to_component_index(i).first;

    for (; cell != endc; ++cell)
    	if (cell->is_locally_owned())
      	{
      		const Point<dim> location = cell->center();

      		for (unsigned int dof = 0 ; dof < dofs_per_cell ; dof++)
      		{
      			value[dof] = initial_boundary->ic(location,component[dof]);
      		}

      		cell->distribute_local_to_global(value,locally_owned_solution);

      	}
}


template<int dim>
void 
Solve_System<dim>::run_time_loop()
{
    const QGauss<dim-1> face_quadrature(1);

		typename DoFHandler<dim>::cell_iterator neighbor;
		typename DoFHandler<dim>::active_cell_iterator endc = dof_handler.end();

    const UpdateFlags face_update_flags =  update_normal_vectors;

    FEFaceValues<dim> fe_v_face(fe, face_quadrature, face_update_flags);

		const unsigned int dofs_per_cell = fe.dofs_per_cell;
		const unsigned int dofs_per_component = dofs_per_cell/n_eqn;

		// only implemented a finite volume scheme
		Assert(dofs_per_component == 1,ExcNotImplemented());
		Assert(dofs_per_cell == n_eqn,ExcNotImplemented());
  		
  		Vector<double> component(dofs_per_cell);
    
    	for (unsigned int i = 0 ; i < dofs_per_cell ; i++)
      		component[i] = fe.system_to_component_index(i).first;

      	std::vector<Vector<double>> component_to_system(n_eqn,
                                                          Vector<double>(dofs_per_component));

      for (unsigned int i = 0 ; i < n_eqn ; i ++)
        for (unsigned int j = 0 ; j < dofs_per_component ; j ++)
          	component_to_system[i](j) = fe.component_to_system_index(i,j); 

        std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      	std::vector<types::global_dof_index> local_dof_indices_neighbor(dofs_per_cell);

      	Vector<double> neighbor_values(dofs_per_cell);

      	double t = 0;
      	while (t < t_end)
      	{
      		if (t+dt > t_end)
      			dt = t_end-t;

          // we need to initialise the iterator again and again
          typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active();

      		for (; cell != endc; ++cell)
      			if (cell->is_locally_owned())
      			{
            // no hanging nodes check
      				Assert(!cell->has_children(), ExcInternalError());
      				
      				cell->get_dof_indices(local_dof_indices);
      				const double volume = cell->measure();

          // loop over the faces, we assume no hanging nodes 
      				for(unsigned int face  = 0; face< GeometryInfo<dim>::faces_per_cell; face++ )
      				{
                fe_v_face.reinit(cell,face);
                
      		   		// normal to the face assuming cartesian grid
      				 	Tensor<1,dim> normal_vec = fe_v_face.normal_vector(0);

      					const double nx = normal_vec[0];
      					const double ny = normal_vec[1];

					     // construct An
      				 	Sparse_Matrix An = system_matrices.Ax * nx + system_matrices.Ay * ny;

      					const typename DoFHandler<dim>::face_iterator face_itr = cell->face(face);

      					if (face_itr->at_boundary())
      					{
      						for (unsigned int m = 0 ; m < An.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(An,m); n ; ++n)
                  {
							// 0 because of finite volume
      							unsigned int dof_sol = local_dof_indices[component_to_system[n.row()](0)];

							// the solution id which meets An
      							unsigned int dof_sol_col = local_dof_indices[component_to_system[n.col()](0)];

							// explicit euler update
      							locally_owned_solution(dof_sol) = locally_relevant_solution(dof_sol) - 
      							                                  dt * n.value() 
                                                      * locally_relevant_solution(dof_sol_col)/volume;
      						}
      					}
      					else
      					{
                  neighbor = cell->neighbor(face);

      						 neighbor->get_dof_indices(local_dof_indices_neighbor);

							     Assert(!neighbor->has_children(), ExcInternalError());

                  for (unsigned int m = 0 ; m < An.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(An,m); n ; ++n)
                  {
								// 0 because of finite volume
      							unsigned int dof_sol = local_dof_indices[component_to_system[n.row()](0)];

								// the solution point which meets An
      							unsigned int dof_sol_col = local_dof_indices[component_to_system[n.col()](0)];
      							unsigned int dof_neighbor_col = local_dof_indices_neighbor[component_to_system[n.col()](0)];

					 		// explicit euler update
      							locally_owned_solution(dof_sol) = locally_relevant_solution(dof_sol) - 
      							dt * n.value() * (locally_relevant_solution(dof_sol_col)
      								+ locally_relevant_solution(dof_neighbor_col))/volume;
      						}
       					}


       				}

      			} // end of loop over cells

      			locally_relevant_solution = locally_owned_solution;
      			t += dt;
      			std::cout << "time: " << t << std::endl;

      		}	

}


template<int dim>
Solve_System<dim>::~Solve_System()
{
	dof_handler.clear();
}

// explicit initiation to avoid linker error
template class Solve_System<2>;


