#include "solve_system_SS_adaptive.h"


template<int dim>
Solve_System_SS_adaptive<dim>::Solve_System_SS_adaptive(std::vector<system_data> &system_mat,
                                                        Triangulation<dim> &triangulation,
							                                         	const int poly_degree,
								                                        ic_bc_base<dim> *ic_bc,
                                                        const unsigned int &maximum_eqn)
:
mpi_comm(MPI_COMM_WORLD),
dof_handler(triangulation),
fe_basic(poly_degree),
fe(fe_basic,maximum_eqn),
initial_boundary(ic_bc),
max_neqn(maximum_eqn),
system_matrices(system_mat),
this_mpi_process(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)),
pout(std::cout,this_mpi_process==0)
{

      develop_neqn();
          
      max_speed = compute_max_speed();
}


template<int dim>
void 
Solve_System_SS_adaptive<dim>::distribute_dofs()
{
      dof_handler.distribute_dofs(fe);
}

template<int dim>
void
Solve_System_SS_adaptive<dim>::allocate_memory()
{
      locally_relevant_solution.reinit(dof_handler.n_dofs());
      locally_owned_solution.reinit(dof_handler.n_dofs());
      locally_owned_residual.reinit(dof_handler.n_dofs());
}

template<int dim>
int
Solve_System_SS_adaptive<dim>::solve_steady_state(Triangulation<dim> &triangulation,
                                                  double &t,
                                                  const std::vector<Vector<double>> &force_vector) 
{
        // an approximation to delta_t
        double min_length = min_h(triangulation);
        Utilities::MPI::min(min_length,mpi_comm);

        Assert(max_speed.size() != 0,ExcNotInitialized());
        current_max_index = current_max_fe_index();
        current_max_speed = max_speed[current_max_fe_index()]; // speed corresponding to the current maximum system we have

        dt = CFL * min_length/(dim * current_max_speed);

        pout << "total cells: " << triangulation.n_global_active_cells() << std::endl;
        pout << "total dofs: " << dof_handler.n_dofs() << std::endl;

        Vector<double> component = return_component();
        std::vector<Vector<double>> component_to_system = return_component_to_system();

        Assert(component.size() !=0 , ExcNotInitialized());
        Assert(component_to_system.size() !=0 , ExcNotInitialized());

        QGauss<dim-1> face_quadrature(1);

        PerCellAssemble per_cell_assemble;
        PerCellAssembleScratch per_cell_assemble_scratch(fe,face_quadrature);

        int step_count = 0;
        double residual_ss = 100;

        std::vector<Vector<double>> g(4);

        for (unsigned int id = 0 ; id < 4 ; id++)
            initial_boundary->bc_inhomo(system_matrices[current_max_index].B[id],id,g[id],t);
        

        while (step_count < 100 || residual_ss > 1e-7 ) // we atleast run till t_end || residual_ss > 1e-8
        {
                WorkStream::run (dof_handler.begin_active(),
                dof_handler.end(),
                std::bind(&Solve_System_SS_adaptive<dim>::assemble_per_cell,
                  this,
                  std::placeholders::_1,
                  std::placeholders::_2,
                  std::placeholders::_3,
                  std::cref(component_to_system),
                  std::cref(t),
                  std::cref(g),
                  std::cref(force_vector)),
                std::bind(&Solve_System_SS_adaptive<dim>::assemble_to_global,
                  this,
                  std::placeholders::_1,
                  std::cref(component_to_system)),
                per_cell_assemble_scratch,
                per_cell_assemble);

                residual_ss = locally_owned_residual.l2_norm();
                locally_relevant_solution = locally_owned_solution;

                t += dt;
                step_count++;

      if (step_count%500 == 0)
        pout<< "residual: " << residual_ss << "time: " << t << std::endl;

    }

    return(step_count);

}

template<int dim>
void
Solve_System_SS_adaptive<dim>::run_time_loop(Triangulation<dim> &triangulation,
                                             const unsigned int &cycle,
                                             const unsigned int &refine_cycles,
                                             double &t,
                                             const std::vector<Vector<double>> &force_vector)
{      
      std::vector<Vector<double>> component_to_system = return_component_to_system();   

      discretization_error = 0;
      locally_relevant_solution.reinit(locally_owned_solution.size());
      locally_owned_residual.reinit(locally_owned_solution.size());
      locally_relevant_solution = locally_owned_solution; // the value set to locally owned solution is set here
    

        int total_steps = 0;
        std::cout << "cycle: " << cycle << std::endl;
        fflush(stdout);
      
        total_steps = solve_steady_state(triangulation,t,force_vector);

        pout << "Steps taken: " << total_steps << std::endl;

}


// prescribe the initial conditions
template<int dim>
void
Solve_System_SS_adaptive<dim>::prescribe_initial_conditions()
{

  PerCellIC percellIC;
  PerCellICScratch percellICScratch;

  Vector<double> component = return_component();
  

  WorkStream::run (dof_handler.begin_active(),
                   dof_handler.end(),
                   std::bind(&Solve_System_SS_adaptive<dim>::compute_ic_per_cell,
                            this,
                            std::placeholders::_1,
                            std::placeholders::_2,
                            std::placeholders::_3,
                            std::cref(component)),
                   std::bind(&Solve_System_SS_adaptive<dim>::copy_ic_to_global,
                            this,
                            std::placeholders::_1),
                   percellICScratch,
                   percellIC);


}


// initial conditions for one given cell
template<int dim>
void
Solve_System_SS_adaptive<dim>::compute_ic_per_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                        PerCellICScratch &scratch,
                                        PerCellIC &data,
                                        const Vector<double> &component)
{
  
          data.dofs_per_cell = fe.dofs_per_cell;
          data.local_dof_indices.resize(data.dofs_per_cell);
          data.local_contri.reinit(data.dofs_per_cell);

          const Point<dim> location = cell->center();

          for (unsigned int dof = 0 ; dof < data.dofs_per_cell ; dof++)
            data.local_contri[dof] = initial_boundary->ic(location,component(dof));
          
          cell->get_dof_indices(data.local_dof_indices);

  
}

template<int dim>
void
Solve_System_SS_adaptive<dim>::copy_ic_to_global(const PerCellIC &data)
{
    for(unsigned int dof = 0 ; dof< data.dofs_per_cell ; dof++)
      locally_owned_solution(data.local_dof_indices[dof]) = data.local_contri(dof);
}


template<int dim>
void 
Solve_System_SS_adaptive<dim>::assemble_per_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                      PerCellAssembleScratch &scratch,
                                      PerCellAssemble &data,
                                      const std::vector<Vector<double>> &component_to_system,
                                      const double &t,
                                      const std::vector<Vector<double>> &g,
                                      const std::vector<Vector<double>> &force_vector)
 {

              // operations to avoid data races
              std::vector<Vector<double>> temp_g = g;

              // no hanging nodes check
              Assert(!cell->has_children(), ExcInternalError());
              
              const unsigned int this_fe_index = cell->user_index();

              data.this_neqn = n_eqn[this_fe_index];
              data.local_dof_indices.resize(fe.dofs_per_cell);
              cell->get_dof_indices(data.local_dof_indices);
              const double volume = cell->measure();

              data.local_contri.reinit(fe.dofs_per_cell);
              data.local_contri = 0;


              // contribution from collisions P
              for (unsigned int m = 0 ; m < system_matrices[this_fe_index].P.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(system_matrices[this_fe_index].P,m); n ; ++n)
                  {
              // 0 because of finite volume
                    unsigned int dof_sol = component_to_system[n.row()](0);

              // the solution id which meets An
                    unsigned int dof_sol_col = data.local_dof_indices[component_to_system[n.col()](0)];

              // explicit euler update, 
                    data.local_contri(dof_sol) +=  n.value() * dt
                                             * locally_relevant_solution(dof_sol_col);

                  }

              if(system_matrices[this_fe_index].have_force)
              {
                Vector<double> force_value(n_eqn[this_fe_index]);
                initial_boundary->force(force_value,force_vector[cell->active_cell_index()],
                                        cell->center(),t);

                for(unsigned int i = 0 ; i < n_eqn[this_fe_index] ; i++)
                    data.local_contri(component_to_system[i](0)) += dt * force_value(i);                
              }



              // loop over the faces, we assume no hanging nodes 
              for(unsigned int face  = 0; face < GeometryInfo<dim>::faces_per_cell; face++ )
              {
                scratch.fe_v_face.reinit(cell,face);                
                // normal to the face assuming cartesian grid
                Tensor<1,dim> normal_vec = scratch.fe_v_face.normal_vector(0);
                Vector<double> temp_normal_vec(2);      // 2 is the current maximum dimension

                for(unsigned int space = 0 ; space < dim; space ++)
                  temp_normal_vec(space) = normal_vec[space];

                const double nx = temp_normal_vec[0];
                const double ny = temp_normal_vec[1];

                const typename DoFHandler<dim>::face_iterator face_itr = cell->face(face);
                const double face_length = return_face_length(face_itr);

                  // construct An of the current cell
                  Sparse_Matrix An_cell = system_matrices[this_fe_index].Ax * nx
                                          + system_matrices[this_fe_index].Ay * ny;

                if (face_itr->at_boundary())
                {

                  const unsigned int bc_id = face_itr->boundary_id();
                  // only compute for times greater than zero, already computed for t= 0 before
                    if (system_matrices[this_fe_index].bc_inhomo_time && t > 1e-16)
                                initial_boundary->bc_inhomo(system_matrices[this_fe_index].B[bc_id],bc_id,
                                                          temp_g[bc_id],t);


                  for (unsigned int m = 0 ; m < An_cell.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(An_cell,m); n ; ++n)
                  {
                      // 0 because of finite volume
                      unsigned int dof_sol = component_to_system[n.row()](0);

                     // the solution id which meets An
                      unsigned int dof_sol_col = data.local_dof_indices[component_to_system[n.col()](0)];

                    // explicit euler update
                      data.local_contri(dof_sol) -=  n.value() * dt
                                             * locally_relevant_solution(dof_sol_col) * face_length/volume;

                  }

                  // contribution from penalty_B
                  for (unsigned int m = 0 ; m < system_matrices[this_fe_index].penalty_B[bc_id].outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(system_matrices[this_fe_index].penalty_B[bc_id],m); n ; ++n)
                  {
                    unsigned int dof_sol = component_to_system[n.row()](0);
                    unsigned int dof_sol_col = data.local_dof_indices[component_to_system[n.col()](0)];
                    data.local_contri(dof_sol) +=  n.value() * dt
                                             * locally_relevant_solution(dof_sol_col) * face_length/volume;

                  }
                  // contribution from penalty * g
                  for (unsigned int m = 0 ; m < system_matrices[this_fe_index].penalty[bc_id].outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(system_matrices[this_fe_index].penalty[bc_id],m); n ; ++n)
                  {
              
                    unsigned int dof_sol = component_to_system[n.row()](0);
                    unsigned int dof_sol_col = data.local_dof_indices[component_to_system[n.col()](0)];
                    data.local_contri(dof_sol) -=  n.value() * dt
                                                  * temp_g[bc_id](n.col()) * face_length/volume;

                  }


                }
                else
                {

                    std::cout << "neighbor in " << std::endl;
                    fflush(stdout);

                   typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face);

                   if(!neighbor->has_children())  // if neighbor is active then either its on the same level or its 
                                              // coarser. 
                   {

                    Assert(!face_itr->has_children(),ExcInternalError());
                    integrate_face(data.local_contri,
                                 neighbor,
                                 system_matrices,
                                 component_to_system,
                                 this_fe_index,
                                 face_length,
                                 volume,
                                 nx,
                                 ny,
                                 An_cell,
                                 data.local_dof_indices);

                   }
                   else
                   {
                     Assert(neighbor->has_children(),ExcInternalError()); // the neighbor should have children in this case
                     Assert(neighbor->n_children() > 0, ExcInternalError());
                     if (dim != 1)
                     {
                     for(unsigned int subface = 0 ; subface < face_itr->n_children() ; subface ++) // loop over the subfaces of the present cell
                     {
                      Assert(subface < 2,ExcInternalError());
                      Assert(dim != 1,ExcMessage("should not have reached here"));

                      const typename DoFHandler<dim>::active_cell_iterator neighbor_child 
                      = cell->neighbor_child_on_subface(face,subface);

                      integrate_face(data.local_contri,
                        neighbor_child,
                        system_matrices,
                        component_to_system,
                        this_fe_index,
                                      face_length/2,      // only for isotropic refinement, the length of the subface
                                      volume,
                                      nx,
                                      ny,
                                      An_cell,
                                      data.local_dof_indices);
                     }
                   }

                     if(dim == 1)
                     {
                      std::cout << "neighbor children in " << std::endl;
                      fflush(stdout);

                      std::cout << "find neighbor in" << std::endl;
                      fflush(stdout);
                      const typename DoFHandler<dim>::cell_iterator 
                                                      neighbor_child = return_child_refined_neighbor(neighbor,cell); 

                      Assert(!neighbor_child->has_children(),ExcOrderingChanged(neighbor_child->n_children()));
                      std::cout << "find neighbor out" << std::endl;
                      fflush(stdout);

                      integrate_face(data.local_contri,
                                      neighbor_child,
                                      system_matrices,
                                      component_to_system,
                                      this_fe_index,
                                      face_length,      // only for isotropic refinement, the length of the subface
                                      volume,
                                      nx,
                                      ny,
                                      An_cell,
                                      data.local_dof_indices);

                      std::cout << "neighbor children out " << std::endl;
                      fflush(stdout);

                    }
                     
                   }

                  std::cout << "neighbor out " << std::endl;
                  fflush(stdout);

                } //end of else


              } 
              //end of loop over the faces


 }


template<int dim>
 void
 Solve_System_SS_adaptive<dim>::integrate_face(Vector<double> &result,
                                               const typename DoFHandler<dim>::cell_iterator &neighbor,
                                              const std::vector<system_data> &system_matrices,
                                              const std::vector<Vector<double>> &component_to_system,
                                              const unsigned int &this_fe_index,
                                              const double &face_length,
                                              const double &volume,
                                              const double &nx,
                                              const double &ny,
                                              const Sparse_Matrix &An_cell,
                                              const std::vector<types::global_dof_index> &local_dof_indices)
 {

                   const unsigned int neighbor_fe_index = neighbor->user_index();
                   const unsigned int dofs_per_cell_neighbor = fe.dofs_per_cell;
                   std::vector<types::global_dof_index> local_dof_indices_neighbor(dofs_per_cell_neighbor);

                   neighbor->get_dof_indices(local_dof_indices_neighbor);


                   Sparse_Matrix An_neighbor = system_matrices[neighbor_fe_index].Ax * nx
                                              + system_matrices[neighbor_fe_index].Ay * ny;

                   Sparse_Matrix An_effective = construct_An_effective(An_cell,An_neighbor);

                   Sparse_Matrix Amod = system_matrices[std::max(this_fe_index,neighbor_fe_index)].Ax_mod * fabs(nx) 
                                        +system_matrices[std::max(this_fe_index,neighbor_fe_index)].Ay_mod * fabs(ny);


                   // contribution from the present cell
                  for (unsigned int m = 0 ; m < An_cell.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(An_cell,m); n ; ++n)
                  {
                    unsigned int dof_sol = component_to_system[n.row()](0);
                    unsigned int dof_sol_col = local_dof_indices[component_to_system[n.col()](0)];
                    result(dof_sol) -=  dt * n.value() * (locally_relevant_solution(dof_sol_col))
                                                   * face_length/(2 * volume);

                  }

                  // we add the diffusion now
                  for (unsigned int m = 0 ; m < Amod.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(Amod,m); n ; ++n)
                      if(n.row() < n_eqn[this_fe_index] && n.col() < n_eqn[this_fe_index])
                  {
                    unsigned int dof_sol = component_to_system[n.row()](0);
                    unsigned int dof_sol_col = local_dof_indices[component_to_system[n.col()](0)];
                    result(dof_sol) -=  dt * n.value() * (locally_relevant_solution(dof_sol_col))
                                                   * face_length/(2 * volume);

                  }

                  // contribution from the neighboring cell
                  for (unsigned int m = 0 ; m < An_effective.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(An_effective,m); n ; ++n)
                  {
                    unsigned int dof_sol = component_to_system[n.row()](0);
                    unsigned int dof_neighbor_col = local_dof_indices_neighbor[component_to_system[n.col()](0)];
                    result(dof_sol) -=  dt * n.value() * (locally_relevant_solution(dof_neighbor_col))
                                                   * face_length/(2 * volume);

                  }

                  // diffusion with the neighbouring cell
                  for (unsigned int m = 0 ; m < Amod.outerSize() ; m++)
                    for (Sparse_Matrix::InnerIterator n(Amod,m); n ; ++n)
                      if(n.row() < n_eqn[this_fe_index] && n.col() < n_eqn[neighbor_fe_index])
                  {
                    unsigned int dof_sol = component_to_system[n.row()](0);
                    unsigned int dof_neighbor_col = local_dof_indices_neighbor[component_to_system[n.col()](0)];
                    result(dof_sol) -=  dt * n.value() * (-locally_relevant_solution(dof_neighbor_col))
                                                   * face_length/(2 * volume);

                  }
 }


template<int dim>
void
Solve_System_SS_adaptive<dim>::assemble_to_global(const PerCellAssemble &data,
                                                  const std::vector<Vector<double>> &component_to_system)
 {
        for (unsigned int i = 0 ; i < data.this_neqn ; i ++)
        {
              locally_owned_solution(data.local_dof_indices[component_to_system[i](0)]) = locally_relevant_solution(data.local_dof_indices[component_to_system[i](0)])
                                                               + data.local_contri(component_to_system[i](0));

              if (i <= 7) // only count till tensor degree 2
                locally_owned_residual(data.local_dof_indices[component_to_system[i](0)]) = fabs(data.local_contri(component_to_system[i](0)))/dt;          
        }

 }

template<int dim>
void 
Solve_System_SS_adaptive<dim>::create_output(const std::string &filename)
{

  typename DoFHandler<dim>::active_cell_iterator endc = dof_handler.end(), cell = dof_handler.begin_active();
  std::vector<Vector<double>> component_to_system = return_component_to_system();
  const unsigned int to_print = 13;   // print till the third order tensor

  FILE *fp;
  fp = fopen(filename.c_str(),"w+");

  AssertThrow(fp!=NULL,ExcFileNotOpen(filename));
  for(; cell != endc ; cell++)
    {
      std::vector<types::global_dof_index> local_dof_indices(fe.dofs_per_cell);
      cell->get_dof_indices(local_dof_indices);

      for (unsigned int space = 0 ; space < dim ; space ++)
        fprintf(fp, "%0.16f\t",cell->center()(space));

      for (unsigned int i = 0 ; i < max_neqn ; i++) 
        if(i < to_print)
      {
        const double sol_value =  locally_owned_solution(local_dof_indices[component_to_system[i](0)]);
        fprintf(fp, "%0.16f\t",sol_value);
      }

      fprintf(fp, "\n");
    }

  fclose(fp);

}

template<int dim>
void 
Solve_System_SS_adaptive<dim>::compute_error_per_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                                      PerCellErrorScratch &scratch,
                                                      PerCellError &data)
{

    data.error_value = 0;
    data.solution_value.reinit(fe.n_components());
    data.exact_solution.reinit(fe.n_components());

    VectorTools::point_value(cell->get_dof_handler(),locally_owned_solution, 
                             cell->center(),data.solution_value); 

    initial_boundary->exact_solution(cell->center(),data.exact_solution,t_end);

    // overwrite the solution value with the error in this cell
    data.solution_value.sadd(1,-1,data.exact_solution);
    
    for (unsigned int i = 0 ; i < data.solution_value.size(); i++) 
      data.error_value = pow(data.solution_value(i),2) +data.error_value;

    data.cell_index = cell->active_cell_index();

    data.volume = cell->measure();
    
  
}

template<int dim>
void
Solve_System_SS_adaptive<dim>::copy_error_to_global(const PerCellError &data)
{
   discretization_error += data.error_value * data.volume;
}


template<int dim>
void
Solve_System_SS_adaptive<dim>::compute_error()
{
  PerCellError per_cell_error;
  PerCellErrorScratch per_cell_error_scratch;

  WorkStream::run ( dof_handler.begin_active(),
                    dof_handler.end(),
                    std::bind(&Solve_System_SS_adaptive<dim>::compute_error_per_cell,
                            this,
                            std::placeholders::_1,
                            std::placeholders::_2,
                            std::placeholders::_3),
                    std::bind(&Solve_System_SS_adaptive<dim>::copy_error_to_global,
                            this,
                            std::placeholders::_1),
                   per_cell_error_scratch,
                   per_cell_error);
  
  discretization_error = sqrt(discretization_error);
}



template<int dim>
std::vector<double>
Solve_System_SS_adaptive<dim>::compute_max_speed()
{
    std::vector<double> max_speed(n_eqn.size());

    for(unsigned int i = 0 ; i < n_eqn.size(); i++)
    {
      EigenSolver<MatrixXd> ES(system_matrices[i].Ax);
      VectorXd vals = ES.pseudoEigenvalueMatrix().diagonal().cwiseAbs();      
      max_speed[i] = vals.maxCoeff();
    }


    return(max_speed);
}


template<>
double
Solve_System_SS_adaptive<2>::min_h(Triangulation<2> &triangulation)
{

  typename Triangulation<2>::active_cell_iterator cell = triangulation.begin_active(),
                                                  endc = triangulation.end();


  double min_length = 10000;

  for(; cell != endc ; cell++)
      for(unsigned int face = 0 ; face < GeometryInfo<2>::faces_per_cell ; face++)
        if (cell->face(face)->measure() < min_length)
          min_length = cell->face(face)->measure();


  return(min_length);
}

template<>
double
Solve_System_SS_adaptive<1>::min_h(Triangulation<1> &triangulation)
{
  return(GridTools::minimal_cell_diameter(triangulation));
}

template<>
double
Solve_System_SS_adaptive<2>::return_face_length(const typename DoFHandler<2>::face_iterator &face_itr)
{
  return(face_itr->measure());
}

template<>
double
Solve_System_SS_adaptive<1>::return_face_length(const typename DoFHandler<1>::face_iterator &face_itr)
{
  return(1);
}

template<>
typename DoFHandler<1>::cell_iterator
Solve_System_SS_adaptive<1>::return_child_refined_neighbor(const typename DoFHandler<1>::cell_iterator &neighbor,
                                                             const typename DoFHandler<1>::active_cell_iterator &cell)
{
  std::vector<double> distances(2);
  Assert(neighbor->n_children() == 2,ExcInternalError());
  std::vector<typename DoFHandler<1>::cell_iterator> neighbor_child(2);


  neighbor_child[0] = neighbor->child(0);
  
  distances[0] = cell->center().distance(neighbor_child[0]->center());

  neighbor_child[1] = neighbor->child(1);
  
  distances[1] = cell->center().distance(neighbor_child[1]->center());

  std::cout << "center neighbor: " << neighbor->center() 
            << " center child 1: " << neighbor_child[0]->center()
            << " cener child 2: " << neighbor_child[1]->center() << std::endl;

  std::cout << "level neighbor: " << neighbor->level() 
            << " level child 1: " << neighbor_child[0]->level()
            << " level child 2: " << neighbor_child[1]->level() << std::endl;
            
  Assert(fabs(distances[0]-distances[1]) > 1e-15,ExcInternalError());
  return(distances[0] >= distances[1] ? neighbor_child[1] : neighbor_child[0]);
  
}

template<int dim>
Solve_System_SS_adaptive<dim>::PerCellAssembleScratch::PerCellAssembleScratch(const FESystem<dim> &fe,
                                                                              const QGauss<dim-1> &quadrature)
:
fe_v_face(fe,quadrature,update_normal_vectors)
{;}

template<int dim>
Solve_System_SS_adaptive<dim>::PerCellAssembleScratch::PerCellAssembleScratch(const PerCellAssembleScratch &scratch)
:
fe_v_face(scratch.fe_v_face.get_fe(),
          scratch.fe_v_face.get_quadrature(),
          update_normal_vectors)
{;}

template<int dim>
void 
Solve_System_SS_adaptive<dim>::develop_neqn()
{
  n_eqn.resize(system_matrices.size());

  for(unsigned int i = 0 ; i < system_matrices.size() ; i++)
    n_eqn[i] = system_matrices[i].Ax.rows();

}

template<int dim>
unsigned int 
Solve_System_SS_adaptive<dim>::current_max_fe_index()
{
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();

    unsigned int fe_index_max = 0;

    for(; cell != endc ; cell++)
        if(cell->user_index() > fe_index_max)
          fe_index_max = cell->user_index();

    return(fe_index_max);

}


template<int dim>
Vector<double>
Solve_System_SS_adaptive<dim>::return_component()
{
  Vector<double> component(fe.dofs_per_cell);
  
  for (unsigned int j = 0 ; j < fe.dofs_per_cell ; j++)
        component(j) = fe.system_to_component_index(j).first;
  

  return(component);
}

template<int dim>
std::vector<Vector<double>>
Solve_System_SS_adaptive<dim>::return_component_to_system()
{
        std::vector<Vector<double>> component_to_system(max_neqn);
        const unsigned int dofs_per_comp = 1;     // considering a finite volume scheme

        for (unsigned int k = 0 ; k < max_neqn ;k ++)
        {
          component_to_system[k].reinit(dofs_per_comp);
              for (unsigned int j = 0 ; j < dofs_per_comp ; j++)
                component_to_system[k](j) = fe.component_to_system_index(k,j);          
        }

      return(component_to_system);
}

// returns the effective An of the neighbor
template<int dim>
Sparse_Matrix
Solve_System_SS_adaptive<dim>::construct_An_effective(const Sparse_Matrix &An_cell,const Sparse_Matrix &An_neighbor)
{
  Assert(An_cell.rows() != 0 ,ExcNotInitialized());
  Assert(An_cell.rows() == An_cell.cols(),ExcNotInitialized());

  Assert(An_neighbor.rows() == An_neighbor.cols(),ExcNotInitialized());
  Assert(An_neighbor.rows() != 0 ,ExcNotInitialized());

  const unsigned int rows_cell = An_cell.rows();
  const unsigned int rows_neighbor = An_neighbor.rows();

  // first we copy the biggest one
  Sparse_Matrix result = rows_cell >= rows_neighbor ? An_cell.block(0,0,rows_cell,rows_neighbor)
                                                    : An_neighbor.block(0,0,rows_cell,rows_neighbor);

  return(result);

}


template<int dim>
Solve_System_SS_adaptive<dim>::~Solve_System_SS_adaptive()
{
  dof_handler.clear();
}

// explicit initiation to avoid linker error
template class Solve_System_SS_adaptive<2>;
template class Solve_System_SS_adaptive<1>;


