#include "include_dealii.h"
using namespace dealii;


template<int dim>
double
function_value(const Point<dim> &p);

template<int dim>
void 
write_solution(const std::string &filename,
               const FiniteElement<dim> &fe,
               const DoFHandler<dim> &dof_handler,
               const Vector<double> &value);

template<int dim>
Vector<double>
get_interpolated_values(const DoFHandler<dim> &dof1,
                        const Vector<double> &value1,
                        const Triangulation<dim> &triangulation);

int main(int argc, char *argv[])
{
	const int dim = 1;

	Triangulation<dim> triangulation;

      //The diagonal of the rectangle is the line joining p1 and p2
    GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(atoi(argv[1]));


   DoFHandler<dim> dof_handler1(triangulation);
   DoFHandler<dim> dof_handler2(triangulation);

   FE_DGQ<dim> fe1(0);
   FE_DGQ<dim> fe2(1);

   dof_handler1.distribute_dofs(fe1);
   dof_handler2.distribute_dofs(fe2);

   Vector<double> value(dof_handler1.n_dofs());

   
   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler1.begin_active(),
   											      endc = dof_handler1.end();

   
   for(; cell != endc ; cell++)
   {	
         std::vector<types::global_dof_index> local_dof_indices(cell->get_fe().dofs_per_cell);
         cell->get_dof_indices(local_dof_indices);
   		value(local_dof_indices[0]) = function_value(cell->center());
   }

    Vector<double> value_interpolate = get_interpolated_values(dof_handler1,value,triangulation);

   write_solution<dim>("lower_order",
               fe1,
               dof_handler1,
               value);

   write_solution<dim>("higher_order",
               fe2,
               dof_handler2,
               value_interpolate);

   dof_handler1.clear();
   dof_handler2.clear();

}

template<int dim>
double
function_value(const Point<dim> &p)
{
	const double x = p[0];
	return(pow(x,2));
}

template<int dim>
void 
write_solution(const std::string &filename,
               const FiniteElement<dim> &fe,
               const DoFHandler<dim> &dof_handler,
               const Vector<double> &value)
{
   QGauss<dim> quadrature(3);
   UpdateFlags update_flag = update_quadrature_points;
   FEValues<dim> fe_values(fe,quadrature,update_flag);


   FILE *fp;

   fp = fopen(filename.c_str(),"w+");

   typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(), 
                                                  endc = dof_handler.end();

   for(; cell != endc ; cell++)
   {
      fe_values.reinit(cell);
      Vector<double> value_q_point(fe.n_components());
      

      std::vector<Point<dim>> q_points;
      q_points = fe_values.get_quadrature_points();

      for(unsigned int i = 0 ; i < q_points.size() ; i++)
      {
         VectorTools::point_value(dof_handler,value,
                                  q_points[i],value_q_point);  

         fprintf(fp, "%f\t",q_points[i](0));

         for(unsigned int j = 0 ; j < value_q_point.size(); i++)
            fprintf(fp, "%f\n",value_q_point(0));
      }
   }

   fclose(fp);

}

template<int dim>
Vector<double>
get_interpolated_values(const DoFHandler<dim> &dof1,
                        const Vector<double> &value1,
                        const Triangulation<dim> &triangulation)
{
   DoFHandler<dim> dof2(triangulation);
   FE_DGQ<dim> fe(1);
   dof2.distribute_dofs(fe);

   Vector<double> value_interpolate(dof2.n_dofs());

   FETools::extrapolate(dof1,
                        value1,
                        dof2,
                        value_interpolate);   

   dof2.clear();
   return(value_interpolate);
}


