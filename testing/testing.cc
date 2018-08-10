#include "include_dealii.h"
using namespace dealii;


template<int dim>
double
function_value(const Point<dim> &p);

template<int dim>
void
write_grid(const std::string &filename,const Triangulation<dim> &triangulation)
{
      std::ofstream out (filename.c_str());
      GridOut grid_out;
      grid_out.write_eps (triangulation, out);
}

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

bool
compare(int i,int j,const std::vector<int> &x)
{
  return(x[i]<=x[j]);
}

double my_divide (double x, double y) {return x/y;}

int main(int argc, char *argv[])
{
	// const int dim = 2;

	// Triangulation<dim> triangulation;

 //      //The diagonal of the rectangle is the line joining p1 and p2
 //  GridGenerator::subdivided_hyper_cube(triangulation,atoi(argv[1]));

 //  typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(),endc = triangulation.end();

 //  for(; cell != endc ; cell++)
 //  {
 //    if(cell->center()(0) > 0.5)
 //    cell->set_refine_flag(RefinementCase<dim>::cut_axis(0));
 //  }

 //  triangulation.execute_coarsening_and_refinement();

 //  cell = triangulation.begin_active();
 //  endc = triangulation.end();

 //  for(; cell != endc ; cell++)
 //   {
 //    std::cout << cell->center() << std::endl;
 //    for(unsigned int face = 0 ; face < GeometryInfo<dim>::faces_per_cell ; face++)
 //      if(!cell->face(face)->at_boundary())
 //    {
 //        typename Triangulation<dim>::cell_iterator neighbor = cell->neighbor(face);
 //        std::cout << "children" << neighbor->n_children() << std::endl; 
 //    }    
 //  }

 //  write_grid<dim>("grid",triangulation);

  std::vector<int> error;
  std::vector<int> index;


  error.push_back(0);
  error.push_back(-1);
  error.push_back(2);
  error.push_back(1);

  index.push_back(0);
  index.push_back(1);
  index.push_back(2);
  index.push_back(3);

  std::vector<int> error_sorted = error;

  std::sort(std::begin(error_sorted),std::end(error_sorted),std::less_equal<int>());

  auto fn_half = std::bind (compare,std::placeholders::_1,std::placeholders::_2,std::cref(error));  

  std::sort(std::begin(index),std::end(index),fn_half);

  for(unsigned int i = 0 ; i < 4 ; i++)
    std::cout << "error: " << error[i] << " sorted: " << error_sorted[i] << " index " << index[i] << std::endl;
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


