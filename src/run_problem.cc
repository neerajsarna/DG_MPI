#include "run_problem.h"

template<int dim>
run_problem<dim>::run_problem(std::vector<system_data> &system_mat_primal,	  // system data
                              std::vector<system_data> &system_mat_error_comp,    // system data to compute error
				  			              std::vector<system_data> &system_mat_adjoint, // adjoint data
							                Triangulation<dim> &triangulation, // triangulation
							                 const int poly_degree,
							                 ic_bc_base<dim> *ic_bc_primal,
					                     ic_bc_base<dim> *ic_bc_adjoint,
                               const std::string &foldername,
                               const unsigned int &max_equations_primal,
                               const unsigned int &max_equations_adjoint,
                               const unsigned int &dim_problem)
:
dummy_dof_handler_grid(triangulation),
dummy_dof_handler_velocity(triangulation),
dummy_fe_grid(FE_DGQ<dim>(1),max_equations_adjoint),
dummy_fe_velocity(FE_DGQ<dim>(0),1),
dim_problem(dim_problem)
{
     AssertThrow(max_equations_primal <= max_equations_adjoint,ExcInternalError()); 
     const unsigned int refinement_type_grid = 0;
     const unsigned int refinement_type_velocity = 1;
     std::string filename;

     switch (refinement_type_grid)
     {
      case 0:
      {
        std::cout << "UNIFORM GRID REFINEMENT...." << std::endl;
        break;
      }
      case 1:
      {
        std::cout << "ADAPTIVE GRID REFINEMENT...." << std::endl;
        break;
      }
      default:
      {
        break;
      }
     }

     dummy_dof_handler_grid.distribute_dofs(dummy_fe_grid);
     dummy_dof_handler_velocity.distribute_dofs(dummy_fe_velocity);
     store_user_index.reinit(dummy_dof_handler_velocity.n_dofs());
     store_user_index = 0;

     fill_user_index_from_index_vector();

      TimerOutput timer (std::cout, TimerOutput::summary,
                     TimerOutput::wall_times);


     error_per_cell_velocity.reinit(triangulation.n_active_cells());
     error_per_cell_velocity = 0; // error per cell in the velocity space

     error_per_cell_grid.reinit(triangulation.n_active_cells());
     error_per_cell_grid_target.reinit(triangulation.n_active_cells());
     error_per_cell_grid = 0;     // error per cell in the discretization
     error_per_cell_grid_target = 0;  // error per cell in the target functional. 

     std::cout << "preparing primal: " << std::endl;
     
	   Solve_System_SS_adaptive<dim> solve_primal(system_mat_primal,
	  								 			  triangulation,
	  								 			  poly_degree,
	  								 			  ic_bc_primal,
                            max_equations_adjoint,
                            dim_problem);

      solve_primal.distribute_dofs();
      solve_primal.allocate_memory();
      solve_primal.prescribe_initial_conditions();

     std::cout << "preparing adjoint: " << std::endl;

	   Solve_System_SS_adaptive<dim> solve_adjoint(system_mat_adjoint,
	  								 			  triangulation,
	  								 			  poly_degree,
	  								 			  ic_bc_adjoint,
                            max_equations_adjoint,
                            dim_problem);

     solve_adjoint.distribute_dofs();
     solve_adjoint.allocate_memory();
     solve_adjoint.prescribe_initial_conditions();

	   //const unsigned int refine_cycles = 8;
     unsigned int cycle = 0;
	   t = 0;						// we solve for the steady state so set t only initially
	   std::vector<Vector<double>> component_to_system = solve_primal.return_component_to_system(); 
     std::vector<Vector<double>> component_to_system_adjoint = solve_adjoint.return_component_to_system(); 
	   std::vector<Vector<double>> temp;
     const unsigned int max_dofs = 6 * 320;



     while(compute_active_dofs(triangulation,solve_primal.n_eqn) <= max_dofs)
	   {

        std::cout << "refinement cycle: " << cycle <<  std::endl;
	   		std::cout << "solving primal: " << std::endl;
        timer.enter_subsection("solve primal");
	   		solve_primal.run_time_loop(triangulation,cycle,t,temp);
        timer.leave_subsection();
   
        // std::string filename = foldername + std::string("/result_cycle") 
        //                       + std::to_string(cycle)
					   //                  + std::string(".txt");

        // solve_primal.create_output(filename);

        if(refinement_type_grid == 1)  // write the grid when doing adaptive refinement
        {
            filename = foldername + "/grid" + std::to_string(cycle);
            write_grid(filename,triangulation);
        }
            

	   		// if(cycle != refine_cycles-1)
	   		// {
          if(refinement_type_grid == 1 || refinement_type_velocity == 1)
          {
	   			  std::cout << "solving adjoint: " << std::endl;
            timer.enter_subsection("solve adjoint");
	   			  solve_adjoint.run_time_loop(triangulation,cycle,t,temp);
            
            timer.leave_subsection();

            filename = foldername + std::string("/resultAdj_cycle") + std::to_string(cycle)
                      + std::string(".txt");
          
            solve_adjoint.create_output(filename);

            filename = foldername + std::string("/result_cycle") + std::to_string(cycle)
                      + std::string(".txt");
          
            solve_primal.create_output(filename);
          }

          //solve_primal.compute_error(); 
          compute_error_in_target(triangulation,
                                  ic_bc_primal,
                                  ic_bc_adjoint,
                                  solve_primal.dof_handler,
                                  solve_primal.locally_owned_solution,
                                  system_mat_adjoint,
                                  t);


          if(refinement_type_velocity == 1)
          {

          timer.enter_subsection("compute velocity error");
          compute_error_velocity(solve_primal.locally_owned_solution,
                        solve_primal.dof_handler,
                        solve_adjoint.locally_owned_solution,
                        solve_adjoint.dof_handler,
                        system_mat_error_comp,
                        ic_bc_primal);
          timer.leave_subsection();
          }

          if(refinement_type_grid == 0)
          {
          timer.enter_subsection("compute discretization error");
          compute_error_grid(solve_primal.locally_owned_solution,
                        solve_primal.dof_handler,
                        solve_adjoint.locally_owned_solution,
                        solve_adjoint.dof_handler,
                        system_mat_primal,
                        ic_bc_primal,
                        triangulation);
          timer.leave_subsection();
          }

          std::cout << "Error in target..." << 
                      error_per_cell_grid_target.mean_value() * error_per_cell_grid_target.size() << std::endl;
          std::cout << "Predicted error in target from grid..." << 
                      error_per_cell_grid.mean_value() * error_per_cell_grid.size() << std::endl;
          std::cout << "Predicted error in target from velocity..." << 
                      error_per_cell_velocity.mean_value() * error_per_cell_velocity.size() << std::endl;
          std::cout << "Total predicted error..." << 
                      error_per_cell_velocity.mean_value() * error_per_cell_velocity.size()
                      + error_per_cell_grid.mean_value() * error_per_cell_grid.size() << std::endl;
          std::cout << "efficiency index grid: " << 
                    error_per_cell_grid_target.mean_value() * error_per_cell_grid_target.size()/(error_per_cell_grid.mean_value() * error_per_cell_grid.size()) << std::endl;
          
          print_fe_index(solve_primal.dof_handler,solve_primal.n_eqn);

          filename = foldername + std::string("/error_observed_cycle") + std::to_string(cycle)
                     + std::string(".txt");
          write_error(filename,triangulation,error_per_cell_grid_target);
          filename = foldername + std::string("/error_grid_cycle") + std::to_string(cycle)
                     + std::string(".txt");
          write_error(filename,triangulation,error_per_cell_grid);
          filename = foldername + std::string("/error_velocity_cycle") + std::to_string(cycle)
                     + std::string(".txt");
          write_error(filename,triangulation,error_per_cell_velocity);
                    filename = foldername + std::string("/fe_index_cycle") + std::to_string(cycle)
                     + std::string(".txt");
          write_fe_index(filename,solve_primal.dof_handler,solve_primal.n_eqn);

          develop_convergence_table(solve_primal.discretization_error,
                                    solve_adjoint.discretization_error,
                                    solve_primal.min_h(triangulation),
                                    triangulation,
                                    solve_primal.n_eqn);

          timer.enter_subsection("perform adaptivity");
          
          update_index_vector(triangulation,refinement_type_velocity,system_mat_primal.size());  // update the user indices based upon the error
          update_grid_refine_flags(triangulation,refinement_type_grid);  
          perform_grid_refinement_and_sol_transfer(triangulation,solve_primal,solve_adjoint); // perform grid refinement and solution transfer
          fill_user_index_from_index_vector();  // update the user indices 

          timer.leave_subsection();

          error_per_cell_velocity.reinit(triangulation.n_active_cells());
          error_per_cell_grid.reinit(triangulation.n_active_cells());
          error_per_cell_grid_target.reinit(triangulation.n_active_cells());

	   		//} // end of if condition

        cycle++;
      
	   } // end of loop over cycles

     print_convergence_table(foldername);
     dummy_dof_handler_grid.clear();
     dummy_dof_handler_velocity.clear();

}


template<int dim>
void 
run_problem<dim>::write_error(const std::string &filename,
                              const Triangulation<dim> &triangulation,
                              const Vector<double> &error_per_cell)
{
	FILE *fp;
	fp = fopen(filename.c_str(),"w+");

	typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(),
                                                    endc = triangulation.end();

	for(; cell != endc ; cell++)
	{
		 for (unsigned int space = 0 ; space < dim ; space ++)
        		fprintf(fp, "%0.16f\t",cell->center()(space));

         fprintf(fp, "%0.16f\n",error_per_cell(cell->active_cell_index()));
	}

	fclose(fp);
}

template<int dim>
void 
run_problem<dim>::write_grid(const std::string &filename,const Triangulation<dim> &triangulation)
{
      std::ofstream out (filename.c_str());
      GridOut grid_out;
      grid_out.write_eps (triangulation, out);
}

template<int dim>
void
run_problem<dim>::write_fe_index(const std::string &filename,
                                 const DoFHandler<dim> &dof_handler,
                                 const std::vector<unsigned int> &n_eqn)
{
  FILE *fp;
  fp = fopen(filename.c_str(),"w+");

  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                     endc = dof_handler.end();

  for(; cell != endc ; cell++)
  {
    for(unsigned int i = 0 ; i < dim; i++)
      fprintf(fp, "%f\t",cell->center()(i));

    fprintf(fp, "%d\t%d\n",cell->user_index(),n_eqn[cell->user_index()]);

  }
         
  

  fclose(fp);
}

template<int dim>
void
run_problem<dim>::print_fe_index(const DoFHandler<dim> &dof_handler,
                                 const std::vector<unsigned int> &n_eqn)
{
  std::cout << "FE index data....." << std::endl;
  typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
  std::vector<unsigned int> fe_indices(n_eqn.size());

  for(; cell!=endc ; cell++)
    fe_indices[cell->user_index()]++;
  

  for(unsigned int i = 0 ; i < fe_indices.size(); i++)
    std::cout << "FE index: " << i << 
                " appeared: " << fe_indices[i] << 
                " equations: " << n_eqn[i] << std::endl;


}

template<int dim>
unsigned int 
run_problem<dim>::compute_active_dofs(const Triangulation<dim> &triangulation,
                                           const std::vector<unsigned int> &n_eqn)
{
  typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();
  unsigned int active_dofs = 0;

  for(; cell != endc ; cell++)
    active_dofs += n_eqn[cell->user_index()];

  return(active_dofs);
}
template<int dim>
void
run_problem<dim>::develop_convergence_table(const double &error_primal,
                                           const double &error_adjoint,
                                           const double &min_h,
                                           const Triangulation<dim> &triangulation,
                                           const std::vector<unsigned int> &n_eqn)
{
  const unsigned int num_dofs = compute_active_dofs(triangulation,n_eqn);

  std::string col1 = "primal_error";
  std::string col2 = "adjoint_error";
  std::string col3 = "min_h";
  std::string col4 = "dofs_primal";

  std::string col5 = "target_error";
  std::string col6 = "predicted_error_grid";
  std::string col7 = "predicted_error_velocity";
  std::string col8 = "corrected_error";

  

  const double error_observed = error_per_cell_grid_target.mean_value() * error_per_cell_grid_target.size();
  const double error_predicted_grid = error_per_cell_grid.mean_value() * error_per_cell_grid.size();
  

  const double error_predicted_velocity = error_per_cell_velocity.mean_value() * error_per_cell_velocity.size();
  const double error_corrected = error_observed - error_predicted_grid - error_predicted_velocity;

  convergence_table.add_value(col1,error_primal);
  convergence_table.add_value(col2,error_adjoint);
  convergence_table.add_value(col3,min_h);
  convergence_table.add_value(col4,num_dofs);
  convergence_table.add_value(col5,fabs(error_observed));
  convergence_table.add_value(col6,fabs(error_predicted_grid));
  convergence_table.add_value(col7,fabs(error_predicted_velocity));
  convergence_table.add_value(col8,fabs(error_corrected));

  convergence_table.set_scientific(col1,true);
  convergence_table.set_scientific(col2,true);
  convergence_table.set_scientific(col3,true);
  convergence_table.set_scientific(col4,true);
  convergence_table.set_scientific(col5,true);
  convergence_table.set_scientific(col6,true);
  convergence_table.set_scientific(col7,true);
  convergence_table.set_scientific(col8,true);
}

template<int dim>
void
run_problem<dim>::print_convergence_table(const std::string &foldername)
{ 
       std::ofstream output_convergence(foldername + std::string("/validate_grid_error/convergence_table_uniformM8.txt"));

      convergence_table.evaluate_convergence_rates("primal_error",
                                                  "dofs_primal",
                                                  ConvergenceTable::reduction_rate_log2,dim_problem);

      convergence_table.evaluate_convergence_rates("adjoint_error",
                                                  "dofs_primal",
                                                  ConvergenceTable::reduction_rate_log2,dim_problem);

      convergence_table.evaluate_convergence_rates("target_error",
                                                  "dofs_primal",
                                                  ConvergenceTable::reduction_rate_log2,dim_problem);

      convergence_table.evaluate_convergence_rates("predicted_error_grid",
                                                  "dofs_primal",
                                                  ConvergenceTable::reduction_rate_log2,dim_problem);

      convergence_table.evaluate_convergence_rates("predicted_error_velocity",
                                                  "dofs_primal",
                                                  ConvergenceTable::reduction_rate_log2,dim_problem);

      convergence_table.evaluate_convergence_rates("corrected_error",
                                                  "dofs_primal",
                                                  ConvergenceTable::reduction_rate_log2,dim_problem);
      convergence_table.write_text(output_convergence);

}

template<int dim>
void 
run_problem<dim>::compute_error(const Vector<double> &primal_solution,
                                const Vector<double> &adjoint_solution,
                                const DoFHandler<dim> &dof_handler_adjoint,
                                const std::vector<system_data> &system_matrices,
                                ic_bc_base<dim> *ic_bc_primal,
                                Vector<double> &error_vector,
                                const unsigned int &quad_points)
{
        std::vector<Vector<double>> component_to_system(dof_handler_adjoint.get_fe().n_components());
        const unsigned int dofs_per_comp = dof_handler_adjoint.get_fe().dofs_per_cell/dof_handler_adjoint.get_fe().n_components();     // considering a finite volume scheme

        for (unsigned int k = 0 ; k < component_to_system.size() ;k ++)
        {
          component_to_system[k].reinit(dofs_per_comp);
              for (unsigned int j = 0 ; j < dofs_per_comp ; j++)
                component_to_system[k](j) = dof_handler_adjoint.get_fe().component_to_system_index(k,j);          
        }



        // one points is enough for constant functions
        QGauss<dim> quadrature_basic(quad_points);
        QGauss<dim-1> face_quadrature_basic(quad_points);

        PerCellError per_cell_error;
        PerCellErrorScratch per_cell_error_scratch(dof_handler_adjoint.get_fe(),
                                                 quadrature_basic,
                                                 face_quadrature_basic);


        typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_adjoint.begin_active(),
                                                       endc = dof_handler_adjoint.end();

        std::vector<Vector<double>> g(4);

        for (unsigned int id = 0 ; id < 4 ; id++) // develop the boundary for the biggest inhomogeneity anyhow
            ic_bc_primal->bc_inhomo(system_matrices[system_matrices.size()-1].B[id],id,g[id],t);


          WorkStream::run(cell,
           endc,
           std::bind(&run_problem<dim>::compute_error_per_cell,
             this,
             std::placeholders::_1,
             std::placeholders::_2,
             std::placeholders::_3,
             std::cref(g),
             std::cref(system_matrices),
             std::cref(adjoint_solution),
             std::cref(primal_solution),
             std::cref(component_to_system),
             ic_bc_primal),
           std::bind(&run_problem<dim>::assemble_to_global,
             this,
             std::placeholders::_1,
             std::ref(error_vector)),
           per_cell_error_scratch,
           per_cell_error);

}

template<int dim>
void
run_problem<dim>::compute_error_velocity(const Vector<double> &primal_solution,
                                const DoFHandler<dim> &dof_handler_primal,
                                Vector<double> &adjoint_solution,
                                const DoFHandler<dim> &dof_handler_adjoint,
                                const std::vector<system_data> &system_matrices,
                                ic_bc_base<dim> *ic_bc_primal)
{
        Vector<double> temp(dof_handler_adjoint.n_dofs());

        // velocity space error computation
        // first interpolate
        FETools::interpolate(dof_handler_primal,
                             primal_solution,
                             dof_handler_adjoint,
                             temp);

        compute_error(temp,
                      adjoint_solution,
                      dof_handler_adjoint,
                      system_matrices,
                      ic_bc_primal,
                      error_per_cell_velocity,
                      1);
}


template<int dim>
void
run_problem<dim>::compute_error_grid(const Vector<double> &primal_solution,
                                const DoFHandler<dim> &dof_handler_primal,
                                Vector<double> &adjoint_solution,
                                const DoFHandler<dim> &dof_handler_adjoint,
                                const std::vector<system_data> &system_matrices,
                                ic_bc_base<dim> *ic_bc_primal,
                                const Triangulation<dim> &triangulation)
{

        Assert(dof_handler_primal.get_fe().n_components() == dummy_dof_handler_grid.get_fe().n_components(),ExcInternalError());
        Assert(dof_handler_primal.get_fe().n_components() == dof_handler_adjoint.get_fe().n_components(),ExcInternalError());

        Vector<double> temp_primal(dummy_dof_handler_grid.n_dofs());
        Vector<double> temp_adjoint(dummy_dof_handler_grid.n_dofs());

        // interpolate the primal solution onto the higher dimensional space. But only interpolate, not extrapolate
        // interpolation provides us with a constant function anyhow.
        FETools::interpolate(dof_handler_primal,
                            primal_solution,
                            dummy_dof_handler_grid,
                            temp_primal);

        // extrapolate the adjiont solution to P1 space
        FETools::extrapolate(dof_handler_adjoint,
                             adjoint_solution,
                             dummy_dof_handler_grid,
                             temp_adjoint);

        // for(unsigned int i = 0 ; i < temp_adjoint.size() ; i++)
        //   temp_adjoint(i) = 1;

        compute_error(temp_primal,
                      temp_adjoint,
                      dummy_dof_handler_grid,
                      system_matrices,
                      ic_bc_primal,error_per_cell_grid,
                      2);
}


template<int dim>
void
run_problem<dim>::assemble_to_global(const PerCellError &data,Vector<double> &input)
 {
        input(data.active_index) = data.local_contri;
 }

 template<int dim>
 void 
 run_problem<dim>::compute_error_per_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                        PerCellErrorScratch &scratch,
                                        PerCellError &data,
                                        const std::vector<Vector<double>> &g,
                                        const std::vector<system_data> &system_matrices,
                                        const Vector<double> &adjoint_solution,
                                        const Vector<double> &primal_solution,
                                        const std::vector<Vector<double>> &component_to_system,
                                        ic_bc_base<dim> *ic_bc_primal)
 {
              std::vector<Vector<double>> temp_g = g;


              const unsigned int num_comp = cell->get_fe().n_components();
              const unsigned int this_fe_index = cell->user_index(); // fe index of the primal solution
              const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
              const unsigned int dofs_per_component = dofs_per_cell/num_comp;
              const unsigned int index = cell->active_cell_index();
              const double volume = cell->measure();

              std::vector<types::global_dof_index> local_dof_indices(cell->get_fe().dofs_per_cell);
              cell->get_dof_indices(local_dof_indices);

              data.local_contri = 0;
              data.active_index = index;

              scratch.fe_v.reinit(cell);

              const std::vector<double> &Jacobians_interior = scratch.fe_v.get_JxW_values();
              const unsigned int total_ngp = Jacobians_interior.size();

             Vector<double> primal_value(num_comp);
             Vector<double> adjoint_value(num_comp);    // a vector to store the value of the adjoint solution 
             Vector<double> temp_force_vec;


              for(unsigned int q = 0 ; q < total_ngp ; q++)
              {

                  get_value_at_quad(local_dof_indices,
                                    component_to_system,
                                    adjoint_value,
                                    dofs_per_component,
                                    adjoint_solution,
                                    scratch.fe_v,
                                    q);
                  
                  get_value_at_quad(local_dof_indices,
                                    component_to_system,
                                    primal_value,
                                    dofs_per_component,
                                    primal_solution,
                                    scratch.fe_v,
                                    q);


                  data.local_contri += xAy(adjoint_value,system_matrices[this_fe_index].P,primal_value)
                                        * Jacobians_interior[q];

                  if(system_matrices[this_fe_index].have_force)
                  {

                      Vector<double> force_value(num_comp);
                      Vector<double> temp(1);

                      ic_bc_primal->force(force_value,temp,
                                              scratch.fe_v.quadrature_point(q),t);

                      
                      data.local_contri += force_value * adjoint_value * Jacobians_interior[q];                
                  }
                  
              }


              for(unsigned int face = 0 ; face < GeometryInfo<dim>::faces_per_cell ; face++)
              {
                scratch.fe_v_face.reinit(cell,face);  // reinitialise for the present face
                const std::vector<double> &Jacobians_face = scratch.fe_v_face.get_JxW_values();
                const unsigned int total_ngp_face = Jacobians_face.size();
                
                // normal to the face assuming cartesian grid
                Tensor<1,dim> normal_vec = scratch.fe_v_face.normal_vector(0);
                Vector<double> temp_normal_vec(2);      // 2 is the current maximum dimension

                for(unsigned int space = 0 ; space < dim; space ++)
                  temp_normal_vec(space) = normal_vec[space];

                const double nx = temp_normal_vec[0];
                const double ny = temp_normal_vec[1];

                const typename Triangulation<dim>::face_iterator face_itr = cell->face(face);

                // construct An of the current cell
                Sparse_Matrix An_cell = system_matrices[this_fe_index].Ax * nx
                                      + system_matrices[this_fe_index].Ay * ny;

                if(face_itr->at_boundary())
                {
                  const unsigned int bc_id = face_itr->boundary_id();

                  if(bc_id == 1 && dim_problem == 1)
                    continue;

                  if(bc_id == 3 && dim_problem == 1)
                    continue;

                  for(unsigned int q = 0 ; q < total_ngp_face ; q++)
                  {

                  get_value_at_quad(local_dof_indices,
                                    component_to_system,
                                    adjoint_value,
                                    dofs_per_component,
                                    adjoint_solution,
                                    scratch.fe_v_face,
                                    q);

                  get_value_at_quad(local_dof_indices,
                                    component_to_system,
                                    primal_value,
                                    dofs_per_component,
                                    primal_solution,
                                    scratch.fe_v_face,
                                    q);

                    // integral for SigmaB
                    data.local_contri += xAy(adjoint_value,system_matrices[this_fe_index].penalty_B[bc_id],
                                          primal_value)
                                        * Jacobians_face[q];

                    // integral for Sigmag
                    data.local_contri -= xAy(adjoint_value,system_matrices[this_fe_index].penalty[bc_id],
                                          temp_g[bc_id])
                                        * Jacobians_face[q];
                    
                    }


                  } // end of if
                  else
                  {
                   typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face);

                   if(!neighbor->has_children()) // having no children
                   {
                    scratch.fe_v_neighbor.reinit(neighbor);
                    error_face(neighbor,
                               cell,
                                primal_solution,
                                adjoint_solution,
                                system_matrices,
                                component_to_system,
                                scratch.fe_v_face,
                                scratch.fe_v,
                                scratch.fe_v_neighbor,
                                data.local_contri);
                     
                  }//end over if neighbor children
                  else
                  {
                    Assert(neighbor->has_children() && dim_problem != 1,ExcInternalError());
                    for(unsigned int subface = 0 ; subface < face_itr->n_children() ; subface ++) // loop over the subfaces of the present cell
                     {
                        Assert(subface < 2,ExcInternalError());
                        const typename DoFHandler<dim>::active_cell_iterator neighbor_child 
                                          = cell->neighbor_child_on_subface(face,subface);

                        scratch.fe_v_subface.reinit(cell,face,subface);
                        scratch.fe_v_neighbor.reinit(neighbor_child);

                        error_face(neighbor_child,
                                   cell,
                                    primal_solution,
                                    adjoint_solution,
                                     system_matrices,
                                     component_to_system,
                                     scratch.fe_v_subface,
                                     scratch.fe_v,
                                     scratch.fe_v_neighbor,
                                     data.local_contri);
                     } // end of loop over the subfaces                      

                     

                  } // 
                }//end over else of interior edges
              } // end over loops on faces
 }

 template<int dim>
 void
 run_problem<dim>::error_face(const typename DoFHandler<dim>::cell_iterator &neighbor,  // iterator of the neighbor
                                const typename DoFHandler<dim>::cell_iterator &cell,  // iterator of the cell
                                const Vector<double> &primal_solution,       // primal solution
                                const Vector<double> &adjoint_solution,                         // the adjoint solution
                                const std::vector<system_data> &system_matrices,          // matrices
                                const std::vector<Vector<double>> &component_to_system,         // component to system index of the adjoint
                                const FEValuesBase<dim> &fe_v_face,  
                                const FEValuesBase<dim> &fe_v,                     
                                const FEValuesBase<dim> &fe_v_neighbor, 
                                double &result)
 {
                    const unsigned int num_comp = cell->get_fe().n_components();
                    const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
                    const unsigned int dofs_per_component = dofs_per_cell/num_comp;

                    const unsigned int neighbor_fe_index = neighbor->user_index(); 
                    const unsigned int this_fe_index = cell->user_index(); 

                    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
                    std::vector<types::global_dof_index> local_dof_indices_neighbor(dofs_per_cell);

                    cell->get_dof_indices(local_dof_indices);
                    neighbor->get_dof_indices(local_dof_indices_neighbor);

                    Vector<double> adjoint_value(num_comp);
                    Vector<double> primal_value(num_comp);
                    Vector<double> primal_value_neighbor(num_comp);

                    // reinitialise the finite element object over the cell
                    get_value_at_quad(local_dof_indices,
                                    component_to_system,
                                    primal_value,
                                    dofs_per_component,
                                    primal_solution,
                                    fe_v,
                                    0);   
                    // since we are dealing with a finite volume primal solution, the value remains the same at all the quadrature 
                    // points
                    get_value_at_quad(local_dof_indices_neighbor,
                                    component_to_system,
                                    primal_value_neighbor,
                                    dofs_per_component,
                                    primal_solution,
                                    fe_v_neighbor,
                                    0);   

                    // VectorTools::point_value(cell->get_dof_handler(),primal_solution,cell->center(),primal_value);
                    // VectorTools::point_value(cell->get_dof_handler(),primal_solution,neighbor->center(),primal_value_neighbor);
                   
                    Tensor<1,dim> normal_vec = fe_v_face.normal_vector(0);

                    Vector<double> temp_normal_vec(2);      // 2 is the current maximum dimension

                    for(unsigned int space = 0 ; space < dim; space ++)
                      temp_normal_vec(space) = normal_vec[space];

                    const double nx = temp_normal_vec[0];
                    const double ny = temp_normal_vec[1];

                    const std::vector<double> &Jacobians_face = fe_v_face.get_JxW_values();
                    const unsigned int total_ngp_face = Jacobians_face.size();

                    // take the bigger one during adaptivity
                   const unsigned int effective_index = std::max(this_fe_index,neighbor_fe_index);

                   Sparse_Matrix An_effective = system_matrices[effective_index].Ax * nx
                                              + system_matrices[effective_index].Ay * ny;

                   Sparse_Matrix Amod = system_matrices[effective_index].Ax_mod * fabs(nx)
                                        + system_matrices[effective_index].Ay_mod * fabs(ny);


                   // over the cell
                   for(unsigned int q = 0 ; q < total_ngp_face ; q++) 
                   {
                      get_value_at_quad(local_dof_indices,
                                    component_to_system,
                                    adjoint_value,
                                    dofs_per_component,
                                    adjoint_solution,
                                    fe_v_face,
                                    q);


                      result += xAy(adjoint_value,An_effective,primal_value)
                                           * Jacobians_face[q]/2; 

                      result -= xAy(adjoint_value,Amod,primal_value)
                                           * Jacobians_face[q]/2;

                      result -= xAy(adjoint_value,An_effective,primal_value_neighbor)
                                           * Jacobians_face[q]/2;

                      result += xAy(adjoint_value,Amod,primal_value_neighbor)
                                           * Jacobians_face[q]/2;

                   }  // end of loop over quad points  
 }

 template<int dim>
void
run_problem<dim>::get_value_at_quad(const std::vector<types::global_dof_index> &local_dof_indices,
                                    const std::vector<Vector<double>> &component_to_system,
                                    Vector<double> &adjoint_value,  // value to be filled
                                    const unsigned int &dofs_per_component,
                                    const Vector<double> &adjoint_solution,
                                    const FEValuesBase<dim> &fe_v,
                                    const unsigned int &q)
{
                  adjoint_value = 0;
                  unsigned int dof_test;

                  for(unsigned int i = 0 ; i < adjoint_value.size(); i++) // loop over the number of equations
                    for(unsigned int j =  0 ; j < dofs_per_component; j++)  // loop over the dofs per component
                    {

                      dof_test = local_dof_indices[component_to_system[i](j)]; // dof value 
                      adjoint_value(i) += adjoint_solution(dof_test) * fe_v.shape_value(j,q); // value of the adjoint
                    }
}

// computes x^T A y
template<int dim>
double
run_problem<dim>::xAy(const Vector<double> &x,
                                const Sparse_Matrix &A,
                                const Vector<double> &y)
{
  const unsigned int size_x = x.size();
  const unsigned int size_y = y.size();

  const unsigned int rows_A = A.rows();
  const unsigned int cols_A = A.cols();

  Assert(size_x != 0,ExcNotInitialized());
  Assert(size_y != 0,ExcNotInitialized());
  Assert(A.rows() != 0,ExcNotInitialized());
  Assert(A.cols() != 0,ExcNotInitialized());

  const unsigned int effective_rows = std::min(rows_A,size_x);
  const unsigned int effective_cols = std::min(cols_A,size_y);
  
  double result = 0;
  Sparse_Matrix temp = A.block(0,0,effective_rows,effective_cols);

  // x^T A y
  for (unsigned int m = 0 ; m < temp.outerSize() ; m++)
          for (Sparse_Matrix::InnerIterator n(temp,m); n ; ++n)
            result += x(n.row()) * n.value() * y(n.col());


  return(result);

}

template<int dim>
void
run_problem<dim>::fill_index_vector_from_user_index()
{
  typename DoFHandler<dim>::active_cell_iterator cell = dummy_dof_handler_velocity.begin_active(),
                                                  endc = dummy_dof_handler_velocity.end();


  for(; cell != endc ; cell++)
  {
    std::vector<types::global_dof_index> local_dof_indices(1);
    cell->get_dof_indices(local_dof_indices);
    store_user_index(local_dof_indices[0]) = cell->user_index();
  }
}

template<int dim>
void
run_problem<dim>::fill_user_index_from_index_vector()
{
   typename DoFHandler<dim>::active_cell_iterator cell = dummy_dof_handler_velocity.begin_active(),
                                                  endc = dummy_dof_handler_velocity.end();


  for(; cell != endc ; cell++)
  {
    std::vector<types::global_dof_index> local_dof_indices(1);
    cell->get_dof_indices(local_dof_indices);
    cell->set_user_index(store_user_index(local_dof_indices[0]));
  } 
}

template<int dim>
void
run_problem<dim>::update_index_vector(Triangulation<dim> &triangulation,
                                      const unsigned int &refinement_type,
                                      const unsigned int &num_systems)
{

  const double frac_refine = 0.3;
  const double frac_coarsen = 0.0;

  switch (refinement_type)
  {
    case 0: // uniform refinement
    {
      typename DoFHandler<dim>::active_cell_iterator cell = dummy_dof_handler_velocity.begin_active(),
                                                  endc = dummy_dof_handler_velocity.end();

      for(; cell != endc ; cell++)
      {
         std::vector<types::global_dof_index> local_dof_indices(1);
         cell->get_dof_indices(local_dof_indices);
         if(store_user_index(local_dof_indices[0]) < num_systems-1) // only increase if we have enough components
          store_user_index(local_dof_indices[0]) += 1; // increase by one
      }
      break;
    }
    case 1:
    {
      GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                      modulus_Vec(error_per_cell_velocity),
                                                      frac_refine,frac_coarsen); 

      typename DoFHandler<dim>::active_cell_iterator cell = dummy_dof_handler_velocity.begin_active(),
                                                  endc = dummy_dof_handler_velocity.end();
      for(; cell != endc ; cell++)
        if(cell->refine_flag_set())
          {
                std::vector<types::global_dof_index> local_dof_indices(1);
                cell->get_dof_indices(local_dof_indices);
                
                if(store_user_index(local_dof_indices[0]) < num_systems-1)
                  store_user_index(local_dof_indices[0]) += 1; // increase by one     

                cell->clear_refine_flag();
          } 
      break;     
    }
    default:
    {
      AssertThrow(1 == 0 ,ExcMessage("should not have reached here"));
      break;
    }
  }
  
}

template<int dim>
void
run_problem<dim>::update_grid_refine_flags(Triangulation<dim> &triangulation,
                                           const unsigned int &refinement_type)
{
  const double frac_refine = 0.3;
  const double frac_coarsen = 0.0;

  switch (refinement_type)
  {
    case 0: // uniform refinement
    {
      typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(),
                                                     endc = triangulation.end();

      for(; cell != endc ; cell++)
      {
        if(dim_problem != 1)
            cell->set_refine_flag();
        else
            cell->set_refine_flag(RefinementCase<dim>::cut_axis(0));
      }

      break;
    }
    case 1:
    {
      GridRefinement::refine_and_coarsen_fixed_number(triangulation,modulus_Vec(error_per_cell_grid),
                                                      frac_refine,frac_coarsen); 

      
    // refine_and_coarsen_cancellation(triangulation,
    //                                  error_per_cell_grid,
    //                                  frac_refine);

      if(dim_problem == 1)
      {
      typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(),
                                                     endc = triangulation.end();

      for(; cell != endc ; cell++)
        if(cell->refine_flag_set())
            cell->set_refine_flag(RefinementCase<dim>::cut_axis(0));
      }

      break;
    }
  }
}

// the following routine takes care of the inter element cancellation
template<int dim>
void 
run_problem<dim>::refine_and_coarsen_cancellation(Triangulation<dim> &triangulation,
                                     Vector<double> &error,
                                     const double &refine_frac)
{
  

  Vector<double> error_sorted = error;  // the array which stores the sorted error
  Vector<int> index_cell(triangulation.n_active_cells());  // array which store the index of the error, error(index_cell) = error_sorted
  const double total_error = error_sorted.mean_value() * error_sorted.size();

  for(unsigned int i = 0 ; i < index_cell.size() ; i++)
      index_cell(i) = i;

  typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active();
  bool neg_error = total_error > 0 ? false : true;  
  // check whether the present error is negative or positive. If negative then sort the error array in 
  // increasing order. If error positive then arrange in decreasing order. 

    unsigned int stop_location = 0;
    double sum = 0;



  if(neg_error)
  {
    auto compare = std::bind (&run_problem<dim>::compare_lesser,this,std::placeholders::_1,std::placeholders::_2,error);   
   
    std::sort(error_sorted.begin(),error_sorted.end(),std::less<double>());
   
   AssertThrow(error_sorted.size() == index_cell.size(), ExcInternalError());
    std::sort(index_cell.begin(),index_cell.end(),compare);
  

   
    
    while(fabs(sum) <= refine_frac * fabs(total_error) && error_sorted(stop_location) < 0)
    {
      sum += error_sorted(stop_location);
   
      stop_location++;
      
    }
  }
  else
  {
    
    auto compare = std::bind (&run_problem<dim>::compare_greater,this,std::placeholders::_1,std::placeholders::_2,error);   
   
    std::sort(error_sorted.begin(),error_sorted.end(),std::greater<double>());
   
    std::sort(index_cell.begin(),index_cell.end(),compare);

   
    
    while(fabs(sum) <= refine_frac * fabs(total_error) && error_sorted(stop_location) > 0)
    {
      
      sum += error_sorted(stop_location);
   
      stop_location++;
      
    }   
  }



    AssertThrow(stop_location >0 && stop_location < error.size(), ExcMessage("location cannot be negative"));

    for(unsigned int loc = 0 ; loc <= stop_location ; loc ++)
    {
      typename Triangulation<dim>::active_cell_iterator temp = triangulation.begin_active();
      std::advance(temp,index_cell(loc));  // go to the cell which has to be refined
      Assert(temp->active(),ExcInternalError());
      temp->set_refine_flag();
    }

}


template<int dim>
void
run_problem<dim>::perform_grid_refinement_and_sol_transfer(Triangulation<dim> &triangulation,
                                                           Solve_System_SS_adaptive<dim> &solve_primal,
                                                           Solve_System_SS_adaptive<dim> &solve_adjoint)
{
         triangulation.prepare_coarsening_and_refinement();

         SolutionTransfer<dim,Vector<double>,DoFHandler<dim>> solution_transfer_primal(solve_primal.dof_handler);
         SolutionTransfer<dim,Vector<double>,DoFHandler<dim>> solution_transfer_adjoint(solve_adjoint.dof_handler);
         SolutionTransfer<dim,Vector<double>,DoFHandler<dim>> solution_transfer_index_vector(dummy_dof_handler_velocity);

         solution_transfer_primal.prepare_for_coarsening_and_refinement(solve_primal.locally_owned_solution);
         solution_transfer_adjoint.prepare_for_coarsening_and_refinement(solve_adjoint.locally_owned_solution);
         solution_transfer_index_vector.prepare_for_coarsening_and_refinement(store_user_index);

         triangulation.execute_coarsening_and_refinement();

         solve_primal.distribute_dofs();
         solve_adjoint.distribute_dofs();
         dummy_dof_handler_velocity.distribute_dofs(dummy_fe_velocity);
         dummy_dof_handler_grid.distribute_dofs(dummy_fe_grid);

         Vector<double> tmp(solve_primal.dof_handler.n_dofs());
         solution_transfer_primal.interpolate(solve_primal.locally_owned_solution, tmp);
         solve_primal.locally_owned_solution = tmp;

         tmp.reinit(solve_adjoint.dof_handler.n_dofs());
         solution_transfer_adjoint.interpolate(solve_adjoint.locally_owned_solution, tmp);
         solve_adjoint.locally_owned_solution = tmp;

         tmp.reinit(dummy_dof_handler_velocity.n_dofs());
         solution_transfer_index_vector.interpolate(store_user_index, tmp);
         store_user_index = tmp;
}

template<int dim>
Vector<double>
run_problem<dim>::modulus_Vec(const Vector<double> &inVec)
{
  Vector<double> outVec(inVec.size());

  for(unsigned int i = 0 ; i < inVec.size() ; i++)
    outVec(i) =fabs(inVec(i));

  return(outVec);
}

// for problems with an exact solution, we compute the error in the target functional
template<int dim>
void 
run_problem<dim>::compute_error_in_target(const Triangulation<dim> &triangulation,
                                     ic_bc_base<dim> *ic_bc_primal,
                                     ic_bc_base<dim> *ic_bc_adjoint,
                                     const DoFHandler<dim> &dof_handler_primal,
                                     const Vector<double> &primal_solution,
                                     const std::vector<system_data> &system_matrices_adjoint,
                                     const double &t)
{
  typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();
  Vector<double> temp(1);

  std::vector<Vector<double>> jbc_Omega(4);


  for (unsigned int id = 0 ; id < 4 ; id++) // develop the boundary for the biggest inhomogeneity anyhow
            ic_bc_adjoint->bc_inhomo(system_matrices_adjoint[system_matrices_adjoint.size()-1].B[id],
                                     id,jbc_Omega[id],t);

  std::vector<int> bc_id_primal(4); // the id of primal which is the id of adjoint
  // bc_id_primal[0] = 2;  // adjoint boundary at x = 1 is the primal boundary at x = 0 (reversal in advection direction)
  // bc_id_primal[1] = 3;
  // bc_id_primal[2] = 0;
  // bc_id_primal[3] = 1;

  bc_id_primal[0] = 0;  // adjoint boundary at x = 1 is the primal boundary at x = 0 (reversal in advection direction)
  bc_id_primal[1] = 1;
  bc_id_primal[2] = 2;
  bc_id_primal[3] = 3;

  for(; cell != endc ; cell++)
  {
    const double volume = cell->measure();
    Vector<double> solution_value(dummy_fe_grid.n_components());
    Vector<double> exact_solution_value(dummy_fe_grid.n_components());
    Vector<double> jOmega(dummy_fe_grid.n_components());


    VectorTools::point_value(dof_handler_primal,primal_solution,cell->center(),solution_value); // compute exact solution
    ic_bc_primal->exact_solution(cell->center(),exact_solution_value,t);
    ic_bc_adjoint->force(jOmega,temp,cell->center(),t);

    for(unsigned int i = 0 ; i < jOmega.size(); i++)
      error_per_cell_grid_target(cell->active_cell_index()) += jOmega(i)
                                                           *(exact_solution_value(i)-solution_value(i)) * volume;

    for(unsigned int face = 0 ; face < GeometryInfo<dim>::faces_per_cell ; face++)
      if(cell->face(face)->at_boundary())
    {
      ic_bc_primal->exact_solution(cell->face(face)->center(),exact_solution_value,t);
      const double face_length = cell->face(face)->measure();
      const unsigned int bc_id = bc_id_primal[cell->face(face)->boundary_id()];
      for(unsigned int i = 0 ; i < jbc_Omega[bc_id].size(); i++)
      {
        error_per_cell_grid_target(cell->active_cell_index()) += jbc_Omega[bc_id](i)
                                                           *(exact_solution_value(i)-solution_value(i)) * face_length;
      }
    }

  }
}

template<int dim>
bool 
run_problem<dim>::compare_greater(int i,int j,const Vector<double> &error)
{
  return(error[i]>error[j]);
}

template<int dim>
bool 
run_problem<dim>::compare_lesser(int i,int j,const Vector<double> &error)
{
  return(error[i]<error[j]);
}


template<int dim>
run_problem<dim>::PerCellErrorScratch::PerCellErrorScratch(const FiniteElement<dim> &fe,
                                     const QGauss<dim> &   quadrature_int,
                                    const QGauss<dim-1> &   quadrature_face)
:
fe_v(fe,quadrature_int,update_values|update_JxW_values|update_quadrature_points),
fe_v_neighbor(fe,quadrature_int,update_values|update_JxW_values|update_quadrature_points),
fe_v_face(fe,quadrature_face,update_values|update_normal_vectors|update_JxW_values),
fe_v_subface(fe,quadrature_face,update_values|update_normal_vectors|update_JxW_values)
{;}

template<int dim>
run_problem<dim>::PerCellErrorScratch::PerCellErrorScratch(const PerCellErrorScratch &scratch)
:
fe_v(scratch.fe_v.get_fe(),scratch.fe_v.get_quadrature(),update_values|update_JxW_values|update_quadrature_points),
fe_v_neighbor(scratch.fe_v_neighbor.get_fe(),scratch.fe_v_neighbor.get_quadrature(),update_values|update_JxW_values|update_quadrature_points),
fe_v_face(scratch.fe_v_face.get_fe(),scratch.fe_v_face.get_quadrature(),update_values|update_normal_vectors|update_JxW_values),
fe_v_subface(scratch.fe_v_subface.get_fe(),scratch.fe_v_subface.get_quadrature(),update_values|update_normal_vectors|update_JxW_values)
{;}


template class run_problem<2>;
