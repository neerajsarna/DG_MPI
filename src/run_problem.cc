#include "run_problem.h"

template<int dim>
run_problem<dim>::run_problem(std::vector<system_data> &system_mat_primal,	  // system data
                              std::vector<system_data> &system_mat_error_comp,    // system data to compute error
				  			              std::vector<system_data> &system_mat_adjoint, // adjoint data
							                Triangulation<dim> &triangulation, // triangulation
							                 const int poly_degree,
							                 ic_bc_base<dim> *ic_bc_primal,
					                     ic_bc_base<dim> *ic_bc_adjoint,
                               const std::string &foldername)
{
      TimerOutput timer (std::cout, TimerOutput::summary,
                     TimerOutput::wall_times);

     error_per_cell.reinit(triangulation.n_active_cells());
     error_per_cell = 0;

	   Solve_System_SS_adaptive<dim> solve_primal(system_mat_primal,
	  								 			  triangulation,
	  								 			  poly_degree,
	  								 			  ic_bc_primal);

      solve_primal.allocate_fe_index(0,error_per_cell,triangulation);
      solve_primal.distribute_dofs();
      solve_primal.allocate_memory();
      solve_primal.prescribe_initial_conditions();

     
	   Solve_System_SS_adaptive<dim> solve_adjoint(system_mat_adjoint,
	  								 			  triangulation,
	  								 			  poly_degree,
	  								 			  ic_bc_adjoint);

     solve_adjoint.allocate_fe_index(0,error_per_cell,triangulation);
     solve_adjoint.distribute_dofs();
     solve_adjoint.allocate_memory();
     solve_adjoint.prescribe_initial_conditions();

	   const unsigned int refine_cycles = 3;
	   t = 0;						// we solve for the steady state so set t only initially
	   std::vector<std::vector<Vector<double>>> component_to_system = solve_primal.return_component_to_system(); 
     std::vector<std::vector<Vector<double>>> component_to_system_adjoint = solve_adjoint.return_component_to_system(); 
	   std::vector<Vector<double>> temp;

	   for(unsigned int cycle = 0 ; cycle < refine_cycles ; cycle++)
	   {

        std::cout << "refinement cycle: " << cycle <<  std::endl;
	   		std::cout << "solving primal: " << std::endl;
        timer.enter_subsection("solve primal");
	   		solve_primal.run_time_loop(triangulation,cycle,refine_cycles,t,temp);
        timer.leave_subsection();
   
        //print_fe_index(solve_primal.dof_handler,solve_primal.n_eqn);

        // std::string filename = foldername + std::string("/fe_index_cycle") + std::to_string(cycle)
        //                                   + "_Kn_" + "0p1" + std::string(".txt");

        // write_fe_index(filename,solve_primal.dof_handler,solve_primal.n_eqn);

        std::string filename = foldername + std::string("/result_cycle") 
                              + std::to_string(cycle)
					                    + std::string(".txt");

        solve_primal.create_output(filename);

        filename = "grid" + std::to_string(cycle);
        write_grid(filename,triangulation);


	   		if(cycle != refine_cycles-1)
	   		{
	   			std::cout << "solving adjoint: " << std::endl;
          timer.enter_subsection("solve adjoint");
	   			solve_adjoint.run_time_loop(triangulation,cycle,refine_cycles,t,solve_primal.cellwise_sol);
          timer.leave_subsection();

          solve_primal.compute_error(); 
          solve_adjoint.compute_error();
          develop_convergence_table(solve_primal.discretization_error,
                                    solve_adjoint.discretization_error,
                                    solve_primal.min_h(triangulation));


          filename = foldername + std::string("/resultAdj_cycle") + std::to_string(cycle)
                 + std::string(".txt");
          
          solve_adjoint.create_output(filename);

          // filename = foldername + std::string("/result_adjoint_cycle") + std::to_string(cycle)
          //                                 + "_Kn_" + "0p1" + std::string(".txt");

          // solve_adjoint.create_output(filename);

       //    timer.enter_subsection("compute velocity error");
	   			// compute_error_velocity(solve_adjoint.dof_handler,
	   			// 			                 system_mat_error_comp,
	   			// 		  	  solve_primal.n_eqn,
	   			// 		  	  solve_adjoint.n_eqn,
	   			// 		  	  solve_primal.initial_boundary,
	   			// 		  	  solve_adjoint.cellwise_sol,
       //              solve_primal.cellwise_sol,
       //              solve_primal.cell_index_center,
       //              solve_adjoint.cell_index_center);

       //    timer.leave_subsection();
       //    solve_primal.prepare_velocity_cycle(cycle,error_per_cell);
       //    solve_adjoint.prepare_velocity_cycle(cycle,error_per_cell);

          // now we conduct the grid adaptivity
          timer.enter_subsection("compute error discretization");
          std::cout << "computing discretization error " << std::endl;
          // compute_error_h(triangulation,
          //                 solve_adjoint.dof_handler,
          //                 solve_adjoint.locally_owned_solution,
          //                 solve_primal.cellwise_sol,
          //                 system_mat_primal,
          //                 solve_adjoint.n_eqn,
          //                 ic_bc_primal,
          //                 solve_primal.cell_index_center);

          // std::cout << "performing refinement.." << std::endl;
          // fflush(stdout);

          filename = foldername + "/error_cycle"  + std::to_string(cycle) + std::string(".txt");
          write_error(filename,triangulation);


         // GridRefinement::refine_and_coarsen_fixed_number (triangulation,
         //                                                  error_per_cell,
         //                                                  0.3, 0.00);

         typename Triangulation<dim>::active_cell_iterator cell =triangulation.begin_active(),
                                                           endc = triangulation.end();

         for(; cell != endc ; cell++)
                cell->set_refine_flag();  // global refinement


         triangulation.prepare_coarsening_and_refinement();

         SolutionTransfer<dim,Vector<double>,hp::DoFHandler<dim>> solution_transfer_primal(solve_primal.dof_handler);
         SolutionTransfer<dim,Vector<double>,hp::DoFHandler<dim>> solution_transfer_adjoint(solve_adjoint.dof_handler);

         solution_transfer_primal.prepare_for_coarsening_and_refinement(solve_primal.locally_owned_solution);
         solution_transfer_adjoint.prepare_for_coarsening_and_refinement(solve_adjoint.locally_owned_solution);

         triangulation.execute_coarsening_and_refinement();

         solve_primal.distribute_dofs();
         solve_adjoint.distribute_dofs();

         Vector<double> tmp(solve_primal.dof_handler.n_dofs());
         solution_transfer_primal.interpolate(solve_primal.locally_owned_solution, tmp);
         solve_primal.locally_owned_solution = tmp;

         tmp.reinit(solve_adjoint.dof_handler.n_dofs());
         solution_transfer_adjoint.interpolate(solve_adjoint.locally_owned_solution, tmp);
         solve_adjoint.locally_owned_solution = tmp;

         error_per_cell.reinit(triangulation.n_active_cells());
         timer.leave_subsection();

	   		} // end of if condition
      
	   } // end of loop over cycles

     print_convergence_table();


}

// template<int dim>
// void 
// run_problem<dim>::compute_error_velocity(const hp::DoFHandler<dim> &dof_handler,
// 										 const std::vector<system_data> &system_matrices,
// 										 const std::vector<unsigned int> &n_eqn_primal,
// 										 const std::vector<unsigned int> &n_eqn_adjoint,
// 										 ic_bc_base<dim> *ic_bc,
//                      const std::vector<Vector<double>> &adjoint_solution,
//                      const std::vector<Vector<double>> &primal_solution,
//                      const std::vector<Point<dim>> &cell_index_primal,
//                      const std::vector<Point<dim>> &cell_index_adjoint)
// {
// 	const typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
// 													 	                               endc = dof_handler.end();												   
// 	error_per_cell = 0;

// 	QGauss<dim-1> face_quadrature_basic(1);
// 	hp::QCollection<dim-1> face_quadrature;

// 	for (unsigned int i = 0 ; i < dof_handler.get_fe().size() ; i++)
// 		face_quadrature.push_back(face_quadrature_basic);

// 	PerCellError per_cell_error;
// 	PerCellErrorScratch per_cell_error_scratch(dof_handler.get_fe(),face_quadrature);

//     std::vector<std::vector<Vector<double>>> g(n_eqn_adjoint.size());

//     for (unsigned int i = 0 ; i < n_eqn_adjoint.size() ; i ++)
//     {
//           g[i].resize(4);
//           for (unsigned int id = 0 ; id < 4 ; id++)
//             ic_bc->bc_inhomo(system_matrices[i].B[id],id,g[i][id],t);
//      }

// 	WorkStream::run(cell,
// 		endc,
// 		std::bind(&run_problem<dim>::compute_error_per_cell_velocity,
// 			this,
// 			std::placeholders::_1,
// 			std::placeholders::_2,
// 			std::placeholders::_3,
// 			std::cref(g),
// 			std::cref(system_matrices),
// 			std::cref(n_eqn_primal),
// 			std::cref(n_eqn_adjoint),
// 			std::cref(adjoint_solution),
// 			std::cref(primal_solution),
//       std::cref(cell_index_primal),
//       std::cref(cell_index_adjoint)),
// 		std::bind(&run_problem<dim>::assemble_to_global,
// 			this,
// 			std::placeholders::_1),
// 		per_cell_error_scratch,
// 		per_cell_error);
// }


// template<int dim>
// void 
// run_problem<dim>::compute_error_per_cell_velocity(const typename hp::DoFHandler<dim>::active_cell_iterator &cell,
//                                       	PerCellErrorScratch &scratch,
//                                       	PerCellError &data,
//                                         const std::vector<std::vector<Vector<double>>> &g,
//                                         const std::vector<system_data> &system_matrices,
//                                         const std::vector<unsigned int> &n_eqn_primal,
// 										                    const std::vector<unsigned int> &n_eqn_adjoint,
//                                         const std::vector<Vector<double>> &adjoint_solution,
//                                         const std::vector<Vector<double>> &primal_solution,
//                                         const std::vector<Point<dim>> &cell_index_primal,
//                                         const std::vector<Point<dim>> &cell_index_adjoint)
//  {
//               // operations to avoid data races, computations are in steady state so we do not need to
//  			        // change the value of temp_g.
//               AssertThrow(cell_index_adjoint[cell->index()].distance(cell_index_primal[cell->index()])<1e-15,
//                         ExcOrderingChanged(cell_index_adjoint[cell->index()].distance(cell_index_primal[cell->index()])));
//               std::vector<std::vector<Vector<double>>> temp_g = g;
//               const unsigned int this_fe_index = cell->active_fe_index();

//               // no hanging nodes check
//               const unsigned int index = cell->index();
//               const double volume = cell->measure();

//               const Vector<double> adjoint_value = adjoint_solution[index];
//               const Vector<double> solution_value = primal_solution[index];
//               data.local_contri = 0;
//               data.index = index;

//               data.local_contri += xAy(adjoint_value,system_matrices[this_fe_index].P,solution_value) * volume;


//               // loop over the faces, we assume no hanging nodes 
//               for(unsigned int face  = 0; face < GeometryInfo<dim>::faces_per_cell; face++ )
//               {
//                 scratch.fe_v_face.reinit(cell,face);
//                 const FEFaceValues<dim> &fe_v_face_temp = scratch.fe_v_face.get_present_fe_values();
                
//                 // normal to the face assuming cartesian grid
//                 Tensor<1,dim> normal_vec = fe_v_face_temp.normal_vector(0);

//                 const double nx = normal_vec[0];
//                 const double ny = normal_vec[1];

//                 const typename Triangulation<dim>::face_iterator face_itr = cell->face(face);
//                 const double face_length = face_itr->measure();

//                   // construct An of the current cell
//                 Sparse_Matrix An_cell = system_matrices[this_fe_index].Ax * nx
//                                           + system_matrices[this_fe_index].Ay * ny;

//                 if (face_itr->at_boundary())
//                 {

//                   const unsigned int bc_id = face_itr->boundary_id();

//                   data.local_contri -= xAy(adjoint_value,An_cell,solution_value) * face_length;
//                   data.local_contri += xAy(adjoint_value,system_matrices[this_fe_index].penalty_B[bc_id],solution_value) * face_length;
//                   data.local_contri -= xAy(adjoint_value,system_matrices[this_fe_index].penalty[bc_id],temp_g[this_fe_index][bc_id]) * face_length;


//                 }
//                 else
//                 {
//                  typename hp::DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face);

//                  if (!neighbor->has_children())
//                  {
//                     const unsigned int index_neighbor = neighbor->index();
//                     const Vector<double> neighbor_value = primal_solution[index_neighbor];

//                     error_face(neighbor,
//                               this_fe_index,
//                               system_matrices,
//                               adjoint_value,
//                               solution_value,
//                               neighbor_value,
//                               face_length,
//                               nx,
//                               ny,
//                               data.local_contri);
                  
//                 } // end of if has childre
//                 else
//                 {
//                     for(unsigned int subface = 0 ; subface < face_itr->n_children() ; subface ++) // loop over the subfaces of the present cell
//                      {
//                         Assert(subface < 2,ExcInternalError());
//                         const typename hp::DoFHandler<dim>::active_cell_iterator neighbor_child 
//                                       = cell->neighbor_child_on_subface(face,subface);

//                         const unsigned int index_neighbor = neighbor_child->index();
//                         const Vector<double> neighbor_value = primal_solution[index_neighbor];

//                         error_face(neighbor_child,
//                               this_fe_index,
//                               system_matrices,
//                               adjoint_value,
//                               solution_value,
//                               neighbor_value,
//                               face_length/2,
//                               nx,
//                               ny,
//                               data.local_contri);
//                      }
//                 }// end of else for children


//                 } //end of else


//               }//end of loop over the faces


//  }

//  template<int dim>
//  void
//  run_problem<dim>::error_face(const typename hp::DoFHandler<dim>::cell_iterator &neighbor,
//                               const unsigned int &this_fe_index,
//                               const std::vector<system_data> &system_matrices,
//                               const Vector<double> &adjoint_value,
//                               const Vector<double> &solution_value,
//                               const Vector<double> &neighbor_value,
//                               const double &face_length,
//                               const double &nx,
//                               const double &ny,
//                               double &result)
//  {
//    const unsigned int neighbor_fe_index = neighbor->active_fe_index();
//    const unsigned int index_neighbor = neighbor->index();
//    const unsigned int effective_index = std::max(this_fe_index,neighbor_fe_index);

//    Sparse_Matrix An_effective = system_matrices[effective_index].Ax * nx
//                                 + system_matrices[effective_index].Ay * ny;

//    Sparse_Matrix Amod = system_matrices[effective_index].Ax_mod * fabs(nx)
//                         + system_matrices[effective_index].Ay_mod * fabs(ny);

//   result -= xAy(adjoint_value,An_effective,solution_value) * face_length/2;
//   result -= xAy(adjoint_value,Amod,solution_value) * face_length/2;
//   result -= xAy(adjoint_value,An_effective,neighbor_value) * face_length/2;
//   result += xAy(adjoint_value,Amod,neighbor_value) * face_length/2;  
// }


// template<int dim>
// void
// run_problem<dim>::assemble_to_global(const PerCellError &data)
//  {
//         error_per_cell(data.index) = fabs(data.local_contri);
//  }

//  template<int dim>
//  void 
//  run_problem<dim>::compute_error_h(const Triangulation<dim> &triangulation,
//                                    const hp::DoFHandler<dim> &dof_handler_adjoint,
//                                    const Vector<double> &value_adjoint,
//                                    const std::vector<Vector<double>> &primal_solution,
//                                    const std::vector<system_data> &system_matrices,
//                                    const std::vector<unsigned int> &n_eqn_adjoint,
//                                    ic_bc_base<dim> *ic_bc,
//                                    const std::vector<Point<dim>> &cell_index_center)
//  {
//     // a memory intensive way of interpolation
//     const unsigned int max_neqn = *std::max_element(std::begin(n_eqn_adjoint),std::end(n_eqn_adjoint));

//     // due to an unavailability of an inbuilt extrapolation function we give up the hp dof handler
//     DoFHandler<dim> temp_dof_handler(triangulation);
//     FE_DGQ<dim> temp_fe_basic(1);  
//     FESystem<dim> temp_fe(temp_fe_basic,max_neqn);

//     temp_dof_handler.distribute_dofs(temp_fe);
//     Vector<double> value_interpolated(temp_dof_handler.n_dofs());
//     Vector<int> cell_fe_index(triangulation.n_active_cells());


//     extrapolate(dof_handler_adjoint,
//                 value_adjoint,
//                 temp_dof_handler,
//                 value_interpolated,
//                 triangulation,
//                 cell_fe_index);

//     // we develop the boundary conditions
//     std::vector<std::vector<Vector<double>>> g(n_eqn_adjoint.size()); // the size of n_eqn_adjoint and n_eqn_primal is the same

//     for (unsigned int i = 0 ; i < n_eqn_adjoint.size() ; i ++)
//     {
//           g[i].resize(4);
//           for (unsigned int id = 0 ; id < 4 ; id++)
//             ic_bc->bc_inhomo(system_matrices[i].B[id],id,g[i][id],t);
//      }


//      // develop the scratch and the data variables
//     QGauss<dim> quadrature_basic(2);
//     QGauss<dim-1> face_quadrature_basic(2);
    

//     PerCellError per_cell_error;
//     PerCellErrorScratch_h per_cell_error_scratch(temp_fe,
//                                                  quadrature_basic,
//                                                  face_quadrature_basic);



//     // compute the system to component index for the adjoint dof_handler
                     
//     std::vector<Vector<double>> component_to_system_adjoint(temp_fe.n_components());
//     const unsigned int dofs_per_comp = 4;     // considering a first order DG scheme

//       // now we allocate the values for component_to_system
//         for (unsigned int k = 0 ; k < temp_fe.n_components() ;k ++)
//         {
//           component_to_system_adjoint[k].reinit(dofs_per_comp);
//             for (unsigned int j = 0 ; j < dofs_per_comp ; j++)
//                 component_to_system_adjoint[k](j) = temp_fe.component_to_system_index(k,j);
//         }


//     const typename DoFHandler<dim>::active_cell_iterator cell = temp_dof_handler.begin_active(),
//                                                          endc = temp_dof_handler.end();      

//     typename DoFHandler<dim>::active_cell_iterator cell_temp = temp_dof_handler.begin_active(),
//                                                          endc_temp = temp_dof_handler.end();      

//     FILE *fp;
//     std::string filename = "interpolated_adj" + std::to_string(triangulation.n_active_cells());
//     fp = fopen(filename.c_str(),"w+");

//     FEValues<dim> fe_values(temp_fe,quadrature_basic,update_quadrature_points);
//    for(; cell_temp != endc_temp ; cell_temp++)
//    {
//       fe_values.reinit(cell_temp);
//       Vector<double> value_q_point(temp_fe.n_components());
      

//       std::vector<Point<dim>> q_points;
//       q_points = fe_values.get_quadrature_points();

//       for(unsigned int i = 0 ; i < q_points.size() ; i++)
//       {
//          VectorTools::point_value(temp_dof_handler,value_interpolated,
//                                   q_points[i],value_q_point);  

//          for(unsigned int space = 0 ; space < dim; space ++)
//             fprintf(fp, "%f\t",q_points[i](space));

         
//          fprintf(fp, "%f\n",value_q_point(0));
//       }
//    }

//    fclose(fp);


//     WorkStream::run(cell,
//     endc,
//     std::bind(&run_problem<dim>::compute_error_per_cell_h,
//       this,
//       std::placeholders::_1,
//       std::placeholders::_2,
//       std::placeholders::_3,
//       std::cref(g),
//       std::cref(system_matrices),
//       std::cref(value_interpolated),
//       std::cref(primal_solution),
//       std::cref(component_to_system_adjoint),
//       std::cref(cell_fe_index),
//       std::cref(cell_index_center)),
//     std::bind(&run_problem<dim>::assemble_to_global,
//       this,
//       std::placeholders::_1),
//     per_cell_error_scratch,
//     per_cell_error);

//      temp_dof_handler.clear();
//   }

//  template<int dim>
//  void 
//  run_problem<dim>::compute_error_per_cell_h(const typename DoFHandler<dim>::active_cell_iterator &cell,
//                                         PerCellErrorScratch_h &scratch,
//                                         PerCellError &data,
//                                         const std::vector<std::vector<Vector<double>>> &g,
//                                         const std::vector<system_data> &system_matrices,
//                                         const Vector<double> &adjoint_solution,
//                                         const std::vector<Vector<double>> &primal_solution,
//                                         const std::vector<Vector<double>> &component_to_system_adjoint,
//                                         const Vector<int> &cell_fe_index,
//                                         const std::vector<Point<dim>> &cell_index_center)
//  {
//                 AssertThrow(cell->center().distance(cell_index_center[cell->index()])<1e-15,
//                         ExcOrderingChanged(cell->center().distance(cell_index_center[cell->index()])));

//               std::vector<std::vector<Vector<double>>> temp_g = g;
//               const unsigned int num_comp = cell->get_fe().n_components();
//               const unsigned int dofs_per_component = cell->get_fe().dofs_per_cell/num_comp;
//               const unsigned int this_fe_index = cell_fe_index(cell->index()); // fe index of the primal solution
//               // no hanging nodes check
//               const unsigned int index = cell->index();
//               const double volume = cell->measure();

//               std::vector<types::global_dof_index> local_dof_indices(cell->get_fe().dofs_per_cell);
//               cell->get_dof_indices(local_dof_indices);

//               data.local_contri = 0;
//               data.index = index;

//               scratch.fe_v.reinit(cell);

//               const std::vector<double> &Jacobians_interior = scratch.fe_v.get_JxW_values();
//               const unsigned int total_ngp = Jacobians_interior.size();

//              const Vector<double> primal_value = primal_solution[data.index]; // primal solution in the present cell
//              Vector<double> adjoint_value(num_comp);    // a vector to store the value of the adjoint solution 
//              adjoint_value = 0;


//               for(unsigned int q = 0 ; q < total_ngp ; q++)
//               {

//                   get_value_at_quad(local_dof_indices,
//                                     component_to_system_adjoint,
//                                     adjoint_value,
//                                     dofs_per_component,
//                                     adjoint_solution,
//                                     scratch.fe_v,
//                                     q);

//                   data.local_contri += xAy(adjoint_value,system_matrices[this_fe_index].P,primal_value)
//                                         * Jacobians_interior[q];
                  
//               }

//               for(unsigned int face = 0 ; face < GeometryInfo<dim>::faces_per_cell ; face++)
//               {
//                 scratch.fe_v_face.reinit(cell,face);  // reinitialise for the present face
//                 const std::vector<double> &Jacobians_face = scratch.fe_v_face.get_JxW_values();
//                 const unsigned int total_ngp_face = Jacobians_face.size();
                
//                 // normal to the face assuming cartesian grid
//                 Tensor<1,dim> normal_vec = scratch.fe_v_face.normal_vector(0);

//                 const double nx = normal_vec[0];
//                 const double ny = normal_vec[1];

//                 const typename Triangulation<dim>::face_iterator face_itr = cell->face(face);

//                 // construct An of the current cell
//                 Sparse_Matrix An_cell = system_matrices[this_fe_index].Ax * nx
//                                       + system_matrices[this_fe_index].Ay * ny;

//                 if(face_itr->at_boundary())
//                 {
//                   const unsigned int bc_id = face_itr->boundary_id();

//                   for(unsigned int q = 0 ; q < total_ngp_face ; q++)
//                   {

//                   get_value_at_quad(local_dof_indices,
//                                     component_to_system_adjoint,
//                                     adjoint_value,
//                                     dofs_per_component,
//                                     adjoint_solution,
//                                     scratch.fe_v_face,
//                                     q);

//                     // integral for SigmaB
//                     data.local_contri += xAy(adjoint_value,system_matrices[this_fe_index].penalty_B[bc_id],
//                                           primal_value)
//                                         * Jacobians_face[q];

//                     // integral for Sigmag
//                     data.local_contri -= xAy(adjoint_value,system_matrices[this_fe_index].penalty[bc_id],
//                                           temp_g[this_fe_index][bc_id])
//                                         * Jacobians_face[q];
                    
//                     }

//                     if (bc_id == 0 || bc_id == 2)
//                     std::cout << "boundary error: " << data.local_contri << "...bc id: " << bc_id <<
//                                  " norm adjoint: " << adjoint_value.l2_norm() <<
//                                  " norm primal: " << primal_value.l2_norm() << 
//                                  " location: " << cell->center() << std::endl;


//                   } // end of if
//                   else
//                   {
//                    typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face);
//                    if(!neighbor->has_children()) // having no children
//                    {

//                     error_face_h(neighbor,
//                                 primal_solution,
//                                 primal_value,
//                                 this_fe_index,
//                                 system_matrices,
//                                 local_dof_indices,
//                                 component_to_system_adjoint,
//                                 adjoint_solution,
//                                 scratch.fe_v_face,
//                                 cell_fe_index,
//                                 data.local_contri);
                     
//                   }//end over if neighbor children
//                   else
//                   {
//                     Assert(neighbor->has_children(),ExcInternalError());
                    
//                     for(unsigned int subface = 0 ; subface < face_itr->n_children() ; subface ++) // loop over the subfaces of the present cell
//                      {
//                         Assert(subface < 2,ExcInternalError());
//                         const typename DoFHandler<dim>::active_cell_iterator neighbor_child 
//                                           = cell->neighbor_child_on_subface(face,subface);

//                         scratch.fe_v_subface.reinit(cell,face,subface);

//                         error_face_h(neighbor_child,
//                                      primal_solution,
//                                      primal_value,
//                                      this_fe_index,
//                                      system_matrices,
//                                      local_dof_indices,
//                                      component_to_system_adjoint,
//                                      adjoint_solution,
//                                      scratch.fe_v_subface,
//                                      cell_fe_index,
//                                      data.local_contri);
//                      } // end of loop over the subfaces
//                   } // 
//                 }//end over else of interior edges
//               } // end over for
//  }

//  template<int dim>
//  void
//  run_problem<dim>::error_face_h(const typename DoFHandler<dim>::cell_iterator &neighbor,  // iterator of the neighbor
//                                 const std::vector<Vector<double>> &primal_solution,       // primal solution
//                                 const Vector<double> &primal_value,                       // value of the current cell
//                                 const unsigned int &this_fe_index,                        // fe index of the present cell
//                                 const std::vector<system_data> &system_matrices,          // matrices
//                                 const std::vector<types::global_dof_index> &local_dof_indices, // local dof indices for solution construction
//                                 const std::vector<Vector<double>> &component_to_system,         // component to system index of the adjoint
//                                 const Vector<double> &adjoint_solution,                         // the adjoint solution
//                                 const FEValuesBase<dim> &fe_v_face,                       
//                                 const Vector<int> &cell_fe_index,
//                                 double &result)
//  {

//                     const unsigned int dofs_per_component = neighbor->get_fe().dofs_per_cell/neighbor->get_fe().n_components();
//                     const unsigned int index_neighbor = neighbor->index();
//                     const unsigned int neighbor_fe_index = cell_fe_index(index_neighbor); 
//                     const Vector<double> neighbor_value = primal_solution[index_neighbor];

//                     Vector<double> adjoint_value_face(neighbor->get_fe().n_components());
                   
//                     Tensor<1,dim> normal_vec = fe_v_face.normal_vector(0);

//                     const double nx = normal_vec[0];
//                     const double ny = normal_vec[1];

//                     const std::vector<double> &Jacobians_face = fe_v_face.get_JxW_values();
//                     const unsigned int total_ngp_face = Jacobians_face.size();

//                     // take the bigger one during adaptivity
//                    const unsigned int effective_index = std::max(this_fe_index,neighbor_fe_index);

//                    Sparse_Matrix An_effective = system_matrices[effective_index].Ax * nx
//                                               + system_matrices[effective_index].Ay * ny;

//                    Sparse_Matrix Amod = system_matrices[effective_index].Ax_mod * fabs(nx)
//                                         + system_matrices[effective_index].Ay_mod * fabs(ny);


//                    // over the cell
//                    for(unsigned int q = 0 ; q < total_ngp_face ; q++) 
//                    {
//                       get_value_at_quad(local_dof_indices,
//                                     component_to_system,
//                                     adjoint_value_face,
//                                     dofs_per_component,
//                                     adjoint_solution,
//                                     fe_v_face,
//                                     q);

//                       result += xAy(adjoint_value_face,An_effective,primal_value)
//                                            * Jacobians_face[q]/2;

//                       result -= xAy(adjoint_value_face,Amod,primal_value)
//                                            * Jacobians_face[q]/2;

//                       result -= xAy(adjoint_value_face,An_effective,neighbor_value)
//                                            * Jacobians_face[q]/2;

//                       result += xAy(adjoint_value_face,Amod,neighbor_value)
//                                            * Jacobians_face[q]/2;

//                    }  // end of loop over quad points  
//  }


//  template<int dim>
// run_problem<dim>::PerCellErrorScratch::PerCellErrorScratch(const hp::FECollection<dim> &fe,
//                                                                  const hp::QCollection<dim-1> &quadrature)
// :
// fe_v_face(fe,quadrature,update_normal_vectors)
// {;}

// template<int dim>
// run_problem<dim>::PerCellErrorScratch::PerCellErrorScratch(const PerCellErrorScratch &scratch)
// :
// fe_v_face(scratch.fe_v_face.get_fe_collection(),
//           scratch.fe_v_face.get_quadrature_collection(),
//           update_normal_vectors)
// {;}

// template<int dim>
// run_problem<dim>::PerCellErrorScratch_h::PerCellErrorScratch_h(const FiniteElement<dim> &fe,
//                                                                 const QGauss<dim> &quadrature_int,
//                                                                  const QGauss<dim-1> &quadrature)
// :
// fe_v_face(fe,quadrature,update_normal_vectors | update_quadrature_points | update_JxW_values),
// fe_v_subface(fe,quadrature,update_normal_vectors | update_quadrature_points | update_JxW_values),
// fe_v(fe,quadrature_int,update_quadrature_points | update_JxW_values)
// {;}

// template<int dim>
// run_problem<dim>::PerCellErrorScratch_h::PerCellErrorScratch_h(const PerCellErrorScratch_h &scratch)
// :
// fe_v_face(scratch.fe_v_face.get_fe(),
//           scratch.fe_v_face.get_quadrature(),
//           update_normal_vectors | update_quadrature_points | update_JxW_values | update_values),
// fe_v_subface(scratch.fe_v_subface.get_fe(),
//           scratch.fe_v_subface.get_quadrature(),
//           update_normal_vectors | update_quadrature_points | update_JxW_values | update_values),
// fe_v(scratch.fe_v.get_fe(),
//           scratch.fe_v.get_quadrature(),
//           update_quadrature_points | update_JxW_values | update_values)
// {;}

// template<int dim>
// Sparse_Matrix
// run_problem<dim>::construct_An_effective(const Sparse_Matrix &An_cell,const Sparse_Matrix &An_neighbor)
// {
//   Assert(An_cell.rows() != 0 ,ExcNotInitialized());
//   Assert(An_cell.rows() == An_cell.cols(),ExcNotInitialized());

//   Assert(An_neighbor.rows() == An_neighbor.cols(),ExcNotInitialized());
//   Assert(An_neighbor.rows() != 0 ,ExcNotInitialized());

//   const unsigned int rows_cell = An_cell.rows();
//   const unsigned int rows_neighbor = An_neighbor.rows();

//   // first we copy the biggest one
//   Sparse_Matrix result = rows_cell >= rows_neighbor ? An_cell.block(0,0,rows_cell,rows_neighbor)
//                                                     : An_neighbor.block(0,0,rows_cell,rows_neighbor);

//   return(result);

// }

template<int dim>
void 
run_problem<dim>::write_error(const std::string &filename,const Triangulation<dim> &triangulation)
{
	FILE *fp;
	fp = fopen(filename.c_str(),"w+");

	typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();

	for(; cell != endc ; cell++)
	{
		 for (unsigned int space = 0 ; space < dim ; space ++)
        		fprintf(fp, "%0.16f\t",cell->center()(space));

         fprintf(fp, "%0.16f\n",error_per_cell(cell->index()));
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

// template<int dim>
// void
// run_problem<dim>::write_fe_index(const std::string &filename,
//                                 const hp::DoFHandler<dim> &dof_handler,
//                                 const std::vector<unsigned int> &n_eqn)
// {
//   FILE *fp;
//   fp = fopen(filename.c_str(),"w+");

//   typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
//                                                      endc = dof_handler.end();

//   for(; cell != endc ; cell++)
//          fprintf(fp, "%d %d\n",cell->active_fe_index(),n_eqn[cell->active_fe_index()]);
  

//   fclose(fp);
// }

// template<int dim>
// void
// run_problem<dim>::print_fe_index(const hp::DoFHandler<dim> &dof_handler,
//                                  const std::vector<unsigned int> &n_eqn)
// {
//   std::cout << "FE index data....." << std::endl;
//   typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), endc = dof_handler.end();
//   std::vector<unsigned int> fe_indices(n_eqn.size());

//   for(; cell!=endc ; cell++)
//     fe_indices[cell->active_fe_index()]++;
  

//   for(unsigned int i = 0 ; i < fe_indices.size(); i++)
//     std::cout << "FE index: " << i << 
//                 " appeared: " << fe_indices[i] << 
//                 " equations: " << n_eqn[i] << std::endl;


// }

// // computes x^T A y
// template<int dim>
// double
// run_problem<dim>::xAy(const Vector<double> &x,
//                                 const Sparse_Matrix &A,
//                                 const Vector<double> &y)
// {
//   const unsigned int size_x = x.size();
//   const unsigned int size_y = y.size();

//   const unsigned int rows_A = A.rows();
//   const unsigned int cols_A = A.cols();

//   Assert(size_x != 0,ExcNotInitialized());
//   Assert(size_y != 0,ExcNotInitialized());
//   Assert(A.rows() != 0,ExcNotInitialized());
//   Assert(A.cols() != 0,ExcNotInitialized());

//   const unsigned int effective_rows = std::min(rows_A,size_x);
//   const unsigned int effective_cols = std::min(cols_A,size_y);
  
//   double result = 0;
//   Sparse_Matrix temp = A.block(0,0,effective_rows,effective_cols);

//   // x^T A y
//   for (unsigned int m = 0 ; m < temp.outerSize() ; m++)
//           for (Sparse_Matrix::InnerIterator n(temp,m); n ; ++n)
//             result += x(n.row()) * n.value() * y(n.col());


//   return(result);

// }

// // return the value of a finite element solution at a given quad point
// template<int dim>
// void
// run_problem<dim>::get_value_at_quad(const std::vector<types::global_dof_index> &local_dof_indices,
//                                     const std::vector<Vector<double>> &component_to_system,
//                                     Vector<double> &adjoint_value,  // value to be filled
//                                     const unsigned int &dofs_per_component,
//                                     const Vector<double> &adjoint_solution,
//                                     const FEValuesBase<dim> &fe_v,
//                                     const unsigned int &q)
// {
//                   adjoint_value = 0;
//                   unsigned int dof_test;

//                   for(unsigned int i = 0 ; i < adjoint_value.size(); i++) // loop over the number of equations
//                     for(unsigned int j =  0 ; j < dofs_per_component; j++)  // loop over the dofs per component
//                     {

//                       dof_test = local_dof_indices[component_to_system[i](j)]; // dof value 
//                       adjoint_value(i) += adjoint_solution(dof_test) * fe_v.shape_value(j,q); // value of the adjoint
//                     }
// }

// template<int dim>
// void
// run_problem<dim>::extrapolate(const hp::DoFHandler<dim> &dofIn,
//                         const Vector<double> &InVec,
//                         const DoFHandler<dim>  &dofOut,
//                         Vector<double> &OutVec,
//                         const Triangulation<dim> &triangulation,
//                         Vector<int> &cell_fe_index)
// {
//       // interpolate to highest fe index so that we can stop using hp::DoFHandler
//       DoFHandler<dim> dof_handler(triangulation);
//       FE_DGQ<dim> fe_basic(0);              // first we copy the data onto the non-hp dof handler.
//       FESystem<dim> fe(fe_basic,dofOut.get_fe().n_components());

//       dof_handler.distribute_dofs(fe);

//       Vector<double> temp(dof_handler.n_dofs());

//       interpolate_to_highest_fe_index(dofIn,
//                              InVec,
//                              dof_handler,
//                              temp,
//                              cell_fe_index);
      
//       // interpolate to P1 DG space

//       FETools::extrapolate(dof_handler,
//                            temp,
//                            dofOut,
//                            OutVec);

//       dof_handler.clear();  // clear the old dof hander, not needed
// }

// template<int dim>
// void
// run_problem<dim>::interpolate_to_highest_fe_index(const hp::DoFHandler<dim> &dofIn,
//                                          const Vector<double> &InVec,
//                                          const DoFHandler<dim>  &dofOut,
//                                          Vector<double> &OutVec,
//                                          Vector<int> &cell_fe_index)
// { 
//   // dofOut already corresponds to the highest order one
//   const unsigned int max_neqn = dofOut.get_fe().n_components();
//   typename hp::DoFHandler<dim>::active_cell_iterator cellIn = dofIn.begin_active(), endc = dofIn.end();
//   typename DoFHandler<dim>::active_cell_iterator cellOut = dofOut.begin_active();

//   Assert(InVec.size() == dofIn.n_dofs(),ExcNotInitialized());
//   Assert(OutVec.size() == dofOut.n_dofs(),ExcNotInitialized());

//   OutVec = 0;

//   for(; cellIn != endc ; cellIn++)
//   {
//     unsigned int dofs_this_cell = cellIn->get_fe().dofs_per_cell;
//     unsigned int this_fe_index = cellIn->active_fe_index();
//     unsigned int neqn_this_cell = cellIn->get_fe().n_components();

//     std::vector<types::global_dof_index> local_dof_indicesIn(dofs_this_cell);
//     std::vector<types::global_dof_index> local_dof_indicesOut(max_neqn);    // the number of dofs are the same for every cell

//     cellIn->get_dof_indices(local_dof_indicesIn);
//     cellOut->get_dof_indices(local_dof_indicesOut);

//     for(unsigned int i = 0 ; i < neqn_this_cell; i++)    
//     {
//       // we assume that the dofs are order in the increasing order
//       const unsigned int dofIn = local_dof_indicesIn[i];
//       const unsigned int dofOut = local_dof_indicesOut[i];

//       // all the higher order components are just set to zero
//       OutVec(dofOut) = InVec(dofIn);
//     }  
    
//     // we store the current fe index  
//     cell_fe_index(cellIn->index()) = this_fe_index;

//     cellOut++;      
//   }

// }

template<int dim>
void
run_problem<dim>::develop_convergence_table(const double &error_primal,
                                            const double &error_adjoint,
                                            const double &min_h)
{
  std::string col1 = "primal error";
  std::string col2 = "adjoint error";
  std::string col3 = "min h";

  convergence_table.add_value(col1,error_primal);
  convergence_table.add_value(col2,error_adjoint);
  convergence_table.add_value(col3,min_h);

  convergence_table.set_scientific(col1,true);
  convergence_table.set_scientific(col2,true);
  convergence_table.set_scientific(col3,true);
}

template<int dim>
void
run_problem<dim>::print_convergence_table()
{
      std::ofstream output_convergence("convergence_table.txt");
      convergence_table.evaluate_convergence_rates("primal error",
                                                  "min h",
                                                  ConvergenceTable::reduction_rate_log2,1);

      convergence_table.evaluate_convergence_rates("adjoint error",
                                                  "min h",
                                                  ConvergenceTable::reduction_rate_log2,1);

      convergence_table.write_text(output_convergence);

}

// template<int dim>
// void
// run_problem<dim>::construct_fe_collection(const FiniteElement<dim> &fe_basic,
//                                                        const std::vector<unsigned int> &n_eqn,
//                                                        hp::FECollection<dim> &fe)
// {
//     std::vector<int> block_structure;

//     Assert(n_eqn.size() != 0,ExcNotInitialized());

//     // first create the block structure for the finite element object 
//     // which will then be used to construct the fe system

//     construct_block_structure(block_structure,n_eqn);

//     std::vector< std::vector< const FiniteElement<dim>* > > fe_list;
//     std::vector<unsigned int> multiplicities;

//     fe_list.resize(block_structure.size());
    
//     for(unsigned int i = 0 ; i < block_structure.size() ; i++)
//       multiplicities.push_back(block_structure[i]);

//     switch (block_structure.size())
//     {
//       case 1:
//       {
//         fe_list[0].push_back(&fe_basic);
//         break;
//       }

//       case 2:
//       {
//         fe_list[0].push_back(&fe_basic);
//         fe_list[0].push_back(new FE_Nothing<dim>());

//         fe_list[1].push_back(&fe_basic);
//         fe_list[1].push_back(&fe_basic);

//         break;
//       }

//       case 3:
//       {
//         fe_list[0].push_back(&fe_basic);
//         fe_list[0].push_back(new FE_Nothing<dim>());
//         fe_list[0].push_back(new FE_Nothing<dim>());

//         fe_list[1].push_back(&fe_basic);
//         fe_list[1].push_back(&fe_basic);
//         fe_list[1].push_back(new FE_Nothing<dim>());

//         fe_list[2].push_back(&fe_basic);
//         fe_list[2].push_back(&fe_basic);
//         fe_list[2].push_back(&fe_basic);
//         break;
//       }

//       case 4:
//       {
//         fe_list[0].push_back(&fe_basic);
//         fe_list[0].push_back(new FE_Nothing<dim>());
//         fe_list[0].push_back(new FE_Nothing<dim>());
//         fe_list[0].push_back(new FE_Nothing<dim>());

//         fe_list[1].push_back(&fe_basic);
//         fe_list[1].push_back(&fe_basic);
//         fe_list[1].push_back(new FE_Nothing<dim>());
//         fe_list[1].push_back(new FE_Nothing<dim>());

//         fe_list[2].push_back(&fe_basic);
//         fe_list[2].push_back(&fe_basic);
//         fe_list[2].push_back(&fe_basic);
//         fe_list[2].push_back(new FE_Nothing<dim>());

//         fe_list[3].push_back(&fe_basic);
//         fe_list[3].push_back(&fe_basic);
//         fe_list[3].push_back(&fe_basic);
//         fe_list[3].push_back(&fe_basic);
//         break;
//       }

//       case 5:
//       {
//         fe_list[0].push_back(&fe_basic);
//         fe_list[0].push_back(new FE_Nothing<dim>());
//         fe_list[0].push_back(new FE_Nothing<dim>());
//         fe_list[0].push_back(new FE_Nothing<dim>());
//         fe_list[0].push_back(new FE_Nothing<dim>());

//         fe_list[1].push_back(&fe_basic);
//         fe_list[1].push_back(&fe_basic);
//         fe_list[1].push_back(new FE_Nothing<dim>());
//         fe_list[1].push_back(new FE_Nothing<dim>());
//         fe_list[1].push_back(new FE_Nothing<dim>());

//         fe_list[2].push_back(&fe_basic);
//         fe_list[2].push_back(&fe_basic);
//         fe_list[2].push_back(&fe_basic);
//         fe_list[2].push_back(new FE_Nothing<dim>());
//         fe_list[2].push_back(new FE_Nothing<dim>());

//         fe_list[3].push_back(&fe_basic);
//         fe_list[3].push_back(&fe_basic);
//         fe_list[3].push_back(&fe_basic);
//         fe_list[3].push_back(&fe_basic);
//         fe_list[3].push_back(new FE_Nothing<dim>());


//         fe_list[4].push_back(&fe_basic);
//         fe_list[4].push_back(&fe_basic);
//         fe_list[4].push_back(&fe_basic);
//         fe_list[4].push_back(&fe_basic);
//         fe_list[4].push_back(&fe_basic);
//         break;
//       }

//       case 6:
//       {
//         fe_list[0].push_back(&fe_basic);
//         fe_list[0].push_back(new FE_Nothing<dim>());
//         fe_list[0].push_back(new FE_Nothing<dim>());
//         fe_list[0].push_back(new FE_Nothing<dim>());
//         fe_list[0].push_back(new FE_Nothing<dim>());
//         fe_list[0].push_back(new FE_Nothing<dim>());

//         fe_list[1].push_back(&fe_basic);
//         fe_list[1].push_back(&fe_basic);
//         fe_list[1].push_back(new FE_Nothing<dim>());
//         fe_list[1].push_back(new FE_Nothing<dim>());
//         fe_list[1].push_back(new FE_Nothing<dim>());
//         fe_list[1].push_back(new FE_Nothing<dim>());

//         fe_list[2].push_back(&fe_basic);
//         fe_list[2].push_back(&fe_basic);
//         fe_list[2].push_back(&fe_basic);
//         fe_list[2].push_back(new FE_Nothing<dim>());
//         fe_list[2].push_back(new FE_Nothing<dim>());
//         fe_list[2].push_back(new FE_Nothing<dim>());

//         fe_list[3].push_back(&fe_basic);
//         fe_list[3].push_back(&fe_basic);
//         fe_list[3].push_back(&fe_basic);
//         fe_list[3].push_back(&fe_basic);
//         fe_list[3].push_back(new FE_Nothing<dim>());
//         fe_list[3].push_back(new FE_Nothing<dim>());


//         fe_list[4].push_back(&fe_basic);
//         fe_list[4].push_back(&fe_basic);
//         fe_list[4].push_back(&fe_basic);
//         fe_list[4].push_back(&fe_basic);
//         fe_list[4].push_back(&fe_basic);
//   fe_list[4].push_back(new FE_Nothing<dim>());
 

//         fe_list[5].push_back(&fe_basic);
//         fe_list[5].push_back(&fe_basic);
//         fe_list[5].push_back(&fe_basic);
//         fe_list[5].push_back(&fe_basic);
//         fe_list[5].push_back(&fe_basic);
//   fe_list[5].push_back(&fe_basic);

//        break;
//       }

//       case 7:
//       {
//         fe_list[0].push_back(&fe_basic);
//         fe_list[0].push_back(new FE_Nothing<dim>());
//         fe_list[0].push_back(new FE_Nothing<dim>());
//         fe_list[0].push_back(new FE_Nothing<dim>());
//         fe_list[0].push_back(new FE_Nothing<dim>());
//         fe_list[0].push_back(new FE_Nothing<dim>());
//         fe_list[0].push_back(new FE_Nothing<dim>());

//         fe_list[1].push_back(&fe_basic);
//         fe_list[1].push_back(&fe_basic);
//         fe_list[1].push_back(new FE_Nothing<dim>());
//         fe_list[1].push_back(new FE_Nothing<dim>());
//         fe_list[1].push_back(new FE_Nothing<dim>());
//         fe_list[1].push_back(new FE_Nothing<dim>());
//         fe_list[1].push_back(new FE_Nothing<dim>());

//         fe_list[2].push_back(&fe_basic);
//         fe_list[2].push_back(&fe_basic);
//         fe_list[2].push_back(&fe_basic);
//         fe_list[2].push_back(new FE_Nothing<dim>());
//         fe_list[2].push_back(new FE_Nothing<dim>());
//         fe_list[2].push_back(new FE_Nothing<dim>());
//         fe_list[2].push_back(new FE_Nothing<dim>());

//         fe_list[3].push_back(&fe_basic);
//         fe_list[3].push_back(&fe_basic);
//         fe_list[3].push_back(&fe_basic);
//         fe_list[3].push_back(&fe_basic);
//         fe_list[3].push_back(new FE_Nothing<dim>());
//         fe_list[3].push_back(new FE_Nothing<dim>());
//         fe_list[3].push_back(new FE_Nothing<dim>());


//         fe_list[4].push_back(&fe_basic);
//         fe_list[4].push_back(&fe_basic);
//         fe_list[4].push_back(&fe_basic);
//         fe_list[4].push_back(&fe_basic);
//         fe_list[4].push_back(&fe_basic);
//         fe_list[4].push_back(new FE_Nothing<dim>());
//         fe_list[4].push_back(new FE_Nothing<dim>());
 

//         fe_list[5].push_back(&fe_basic);
//         fe_list[5].push_back(&fe_basic);
//         fe_list[5].push_back(&fe_basic);
//         fe_list[5].push_back(&fe_basic);
//         fe_list[5].push_back(&fe_basic);
//         fe_list[5].push_back(&fe_basic);
//         fe_list[5].push_back(new FE_Nothing<dim>());

//         fe_list[6].push_back(&fe_basic);
//         fe_list[6].push_back(&fe_basic);
//         fe_list[6].push_back(&fe_basic);
//         fe_list[6].push_back(&fe_basic);
//         fe_list[6].push_back(&fe_basic);
//         fe_list[6].push_back(&fe_basic);
//         fe_list[6].push_back(&fe_basic);

//        break;
//       }

//       case 8:
//       {
//         fe_list[0].push_back(&fe_basic);
//         fe_list[0].push_back(new FE_Nothing<dim>());
//         fe_list[0].push_back(new FE_Nothing<dim>());
//         fe_list[0].push_back(new FE_Nothing<dim>());
//         fe_list[0].push_back(new FE_Nothing<dim>());
//         fe_list[0].push_back(new FE_Nothing<dim>());
//         fe_list[0].push_back(new FE_Nothing<dim>());
//         fe_list[0].push_back(new FE_Nothing<dim>());

//         fe_list[1].push_back(&fe_basic);
//         fe_list[1].push_back(&fe_basic);
//         fe_list[1].push_back(new FE_Nothing<dim>());
//         fe_list[1].push_back(new FE_Nothing<dim>());
//         fe_list[1].push_back(new FE_Nothing<dim>());
//         fe_list[1].push_back(new FE_Nothing<dim>());
//         fe_list[1].push_back(new FE_Nothing<dim>());
//         fe_list[1].push_back(new FE_Nothing<dim>());

//         fe_list[2].push_back(&fe_basic);
//         fe_list[2].push_back(&fe_basic);
//         fe_list[2].push_back(&fe_basic);
//         fe_list[2].push_back(new FE_Nothing<dim>());
//         fe_list[2].push_back(new FE_Nothing<dim>());
//         fe_list[2].push_back(new FE_Nothing<dim>());
//         fe_list[2].push_back(new FE_Nothing<dim>());
//         fe_list[2].push_back(new FE_Nothing<dim>());

//         fe_list[3].push_back(&fe_basic);
//         fe_list[3].push_back(&fe_basic);
//         fe_list[3].push_back(&fe_basic);
//         fe_list[3].push_back(&fe_basic);
//         fe_list[3].push_back(new FE_Nothing<dim>());
//         fe_list[3].push_back(new FE_Nothing<dim>());
//         fe_list[3].push_back(new FE_Nothing<dim>());
//         fe_list[3].push_back(new FE_Nothing<dim>());


//         fe_list[4].push_back(&fe_basic);//1
//         fe_list[4].push_back(&fe_basic);//2
//         fe_list[4].push_back(&fe_basic);//3
//         fe_list[4].push_back(&fe_basic);//4
//         fe_list[4].push_back(&fe_basic);//5
//         fe_list[4].push_back(new FE_Nothing<dim>());//6
//         fe_list[4].push_back(new FE_Nothing<dim>());//7
//         fe_list[4].push_back(new FE_Nothing<dim>());//8
 

//         fe_list[5].push_back(&fe_basic);
//         fe_list[5].push_back(&fe_basic);
//         fe_list[5].push_back(&fe_basic);
//         fe_list[5].push_back(&fe_basic);
//         fe_list[5].push_back(&fe_basic);
//         fe_list[5].push_back(&fe_basic);
//         fe_list[5].push_back(new FE_Nothing<dim>());
//         fe_list[5].push_back(new FE_Nothing<dim>());

//         fe_list[6].push_back(&fe_basic);
//         fe_list[6].push_back(&fe_basic);
//         fe_list[6].push_back(&fe_basic);
//         fe_list[6].push_back(&fe_basic);
//         fe_list[6].push_back(&fe_basic);
//         fe_list[6].push_back(&fe_basic);
//         fe_list[6].push_back(&fe_basic);
//         fe_list[6].push_back(new FE_Nothing<dim>());


//         fe_list[7].push_back(&fe_basic);
//         fe_list[7].push_back(&fe_basic);
//         fe_list[7].push_back(&fe_basic);
//         fe_list[7].push_back(&fe_basic);
//         fe_list[7].push_back(&fe_basic);
//         fe_list[7].push_back(&fe_basic);
//         fe_list[7].push_back(&fe_basic);
//         fe_list[7].push_back(&fe_basic);

//        break;
//       }

//       default:
//       {
//         AssertThrow(1 == 0, ExcMessage("Should not have reached here"));
//         break;
//       }
//     }

//     for(unsigned int i = 0 ; i < fe_list.size() ;i++)
//       fe.push_back(FESystem<dim>(fe_list[i],multiplicities));
// }

// template<int dim>
// void 
// run_problem<dim>::construct_block_structure(std::vector<int> &block_structure,
//                                             const std::vector<unsigned int> &n_eqn)
// {
//     AssertThrow(n_eqn.size() != 0, ExcNotInitialized());
//     AssertThrow(std::is_sorted(std::begin(n_eqn),std::end(n_eqn)),
//                 ExcMessage("number of equations not sorted"));

//     // the very first entry should be the number of equations in the first system
//     block_structure.push_back(n_eqn[0]);

//     for (unsigned long int i = 1 ; i < n_eqn.size() ; i++)
//       block_structure.push_back(n_eqn[i]-n_eqn[i-1]);

//     AssertDimension(block_structure.size(),n_eqn.size());
// }

// template<int dim>
// void 
// run_problem<dim>::count_bc_id(const Triangulation<2> &triangulation)
// {
//     typename Triangulation<2>::active_cell_iterator cell = triangulation.begin_active(),
//                                endc = triangulation.end();
    
//     int id0 = 0,id1 = 0, id2 = 0, id3 = 0; 

//     for(; cell != endc ; cell++)   
//       if (cell->is_locally_owned())
//     {

//       for(unsigned int face = 0 ; face < GeometryInfo<2>::faces_per_cell ; face++)
//         if (cell->face(face)->at_boundary())
//           switch (cell->face(face)->boundary_id())
//           {
//             case 0:
//             {
//               id0++;
//               break;
//             }

//             case 1:
//             {
//               id1++;
//               break;
//             }

//             case 2:
//             {
//               id2++;
//               break;
//             }

//             case 3:
//             {
//               id3++;
//               break;
//             }
//           }
//         }
        
    

//           std::cout << "id0 " << id0 
//           << " id1 " << id1 
//           << " id2 " << id2 
//           << " id3 " << id3 << std::endl;

//         fflush(stdout);


// }

template class run_problem<2>;
