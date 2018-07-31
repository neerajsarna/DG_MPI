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


    hp::DoFHandler<dim> dof_handler(triangulation);
    hp::FECollection<dim> fe;

    fe.push_back(FE_DGQ<dim>(0));
    
    typename hp::DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                       endc = dof_handler.end();

    dof_handler.distribute_dofs(fe);

    Vector<double> value(dof_handler.n_dofs());

    for(; cell != endc ; cell++)
    {
      
      std::vector<types::global_dof_index> local_dof_indices(1);
      cell->get_dof_indices(local_dof_indices);

      value(local_dof_indices[0]) = function_value(cell->center());
    }


    hp::DoFHandler<dim> dof2(triangulation);
   hp::FECollection<dim> fe2;

   FE_DGQ<dim> fe_basic2(1);

   fe2.push_back(fe_basic2);
   dof2.distribute_dofs(fe2);

   Vector<double> value_interpolate(dof2.n_dofs());

   FETools::interpolate(dof_handler,
                        value,
                        dof2,
                        value_interpolate);   

   std::cout << "lower_order" << std::endl;
   std::cout << value << std::endl;

   std::cout << "higher_order" << std::endl;
   std::cout << value_interpolate << std::endl;

    dof_handler.clear();
    dof2.clear();


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
   hp::DoFHandler<dim> dof2(triangulation);
   hp::FECollection<dim> fe;

   FE_DGQ<dim> fe_basic(1);

   fe.push_back(fe_basic);
   dof2.distribute_dofs(fe);

   Vector<double> value_interpolate(dof2.n_dofs());

   FETools::extrapolate(dof1,
                        value1,
                        dof2,
                        value_interpolate);   

   dof2.clear();
   return(value_interpolate);
}


