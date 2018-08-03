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
	const int dim = 2;

	Triangulation<dim> triangulation;

      //The diagonal of the rectangle is the line joining p1 and p2
    GridGenerator::subdivided_hyper_cube(triangulation,atoi(argv[1]));
  

    DoFHandler<dim> dof_handler(triangulation);
    FE_DGQ<dim> fe(30);

    dof_handler.distribute_dofs(fe);

    std::cout << "memory consumed(GB): " << dof_handler.memory_consumption()/1e9 << std::endl;
    Vector<double> vec1(dof_handler.n_dofs());
    Vector<double> vec2(dof_handler.n_dofs());
    Vector<double> vec3(dof_handler.n_dofs());
    Vector<double> vec4(dof_handler.n_dofs());
    Vector<double> vec5(dof_handler.n_dofs());
    Vector<double> vec6(dof_handler.n_dofs());

    std::cout << "memory consumed(GB): " << vec1.memory_consumption()/1e9 << std::endl;

    dof_handler.clear();
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


