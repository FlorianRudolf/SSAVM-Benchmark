#include "viennagrid/algorithm/quantity_interpolate.hpp"
#include "viennagrid/algorithm/intersect.hpp"
#include "viennagrid/algorithm/inclusion.hpp"
#include "viennagrid/algorithm/distance.hpp"
#include "viennagrid/algorithm/closest_points.hpp"
#include "viennagrid/io/vtk_reader.hpp"

#include <tclap/CmdLine.h>
#include "common.hpp"


#include "vtk_reader.hpp"






class triangle_wrapper
{
public:

  typedef ElementType element_type;
  typedef viennagrid::result_of::point<ElementType>::type point_type;

  triangle_wrapper(ElementType const & triangle_) : triangle(triangle_)
  {
    point_type p[3];
    p[0] = viennagrid::get_point( viennagrid::vertices(triangle)[0] );
    p[1] = viennagrid::get_point( viennagrid::vertices(triangle)[1] );
    p[2] = viennagrid::get_point( viennagrid::vertices(triangle)[2] );
    std::pair<point_type, point_type> bounding_box = viennagrid::bounding_box(p, p+3);

    center = (bounding_box.first + bounding_box.second)/2.0;
    offset = (bounding_box.second - center) * 1.01;
  }

  ElementType const & operator()() const { return triangle; }

  template<typename BoxT>
  bool intersect(BoxT const & box) const
  {
    point_type other_center = (box.max() + box.min())/2.0;
    point_type other_offset = (other_center - box.min())*1.01;

    if (triangle.is_triangle())
    {
      if ( std::abs(center[0] - other_center[0]) > (offset[0] + other_offset[0]) ) return false;
      if ( std::abs(center[1] - other_center[1]) > (offset[1] + other_offset[1]) ) return false;

      return true;
    }

    assert(false);
    return false;
  }

  template<typename StreamT>
  void print(StreamT & s) const
  {
    s << triangle;
  }

private:
  ElementType triangle;
  point_type center;
  point_type offset;
  viennagrid_numeric tolerance;
};


typedef triangle_wrapper WrapperType;
typedef viennagrid::ntree_node<WrapperType> NTreeNode;


viennagrid_numeric get_value(viennagrid::mesh const & mesh, NTreeNode * ntree, viennagrid::quantity_field const & solution,
                             PointType const & point, viennagrid_numeric epsilon)
{

  NTreeNode const * cur_node = ntree->get(point);
  std::vector<triangle_wrapper> const & wrappers = cur_node->elements();

  std::pair<PointType, PointType> best_cp;
  ElementType best_cell;
  viennagrid_numeric best_distance = -1;

  for (std::size_t i = 0; i != wrappers.size(); ++i)
  {
    std::pair<PointType, PointType> cp = viennagrid::closest_points(wrappers[i](), point);
    viennagrid_numeric distance = viennagrid::distance(cp.first, cp.second);
    if ((best_distance < 0) || (distance < best_distance))
    {
      best_cp = cp;
      best_distance = distance;
      best_cell = wrappers[i]();
    }
  }


//   typedef viennagrid::result_of::const_element_range<MeshType>::type          ConstElementRangeType;
//   typedef viennagrid::result_of::iterator<ConstElementRangeType>::type        ConstElementRangeIterator;
//
//   std::pair<PointType, PointType> best_cp;
//   ElementType best_cell;
//   viennagrid_numeric best_distance = -1;
//
// //   for (viennagrid_dimension dim = 0; dim != viennagrid::cell_dimension(mesh)+1; ++dim)
//   {
//     ConstElementRangeType elements(mesh, viennagrid::cell_dimension(mesh));
//     for (ConstElementRangeIterator cit = elements.begin(); cit != elements.end(); ++cit)
//     {
//       std::pair<PointType, PointType> cp = viennagrid::closest_points(*cit, point);
//       viennagrid_numeric distance = viennagrid::distance(cp.first, cp.second);
//       if ((best_distance < 0) || (distance < best_distance))
//       {
//         best_cp = cp;
//         best_distance = distance;
//         best_cell = *cit;
//       }
//     }
//   }

  return viennagrid::interpolate(best_cell, (best_cp.first+best_cp.second)/2, solution);
}


viennagrid_numeric compare_on_vertices(viennagrid::mesh const & mesh1, NTreeNode * ntree1, viennagrid::quantity_field const & solution1,
                                       viennagrid::mesh const & mesh2, NTreeNode * ntree2, viennagrid::quantity_field const & solution2,
                                       viennagrid_numeric inside_epsilon)
{
  typedef viennagrid::result_of::const_element_range<MeshType>::type          ConstElementRangeType;
  typedef viennagrid::result_of::iterator<ConstElementRangeType>::type        ConstElementRangeIterator;

  viennagrid_numeric worst_distance = -1;

  ConstElementRangeType vertices1(mesh1, 0);
  for (ConstElementRangeIterator vit = vertices1.begin(); vit != vertices1.end(); ++vit)
  {
    PointType point = viennagrid::get_point(*vit);
    viennagrid_numeric value1 = solution1.get(*vit);
    viennagrid_numeric value2 = get_value(mesh2, ntree2, solution2, point, inside_epsilon);

    if ((value1 < 0) || (value2 < 0))
      continue;

    viennagrid_numeric distance = std::abs(value1 - value2);

    if (worst_distance < distance)
      worst_distance = distance;
  }

  return worst_distance;
}


viennagrid_numeric compare(viennagrid::mesh const & mesh1, NTreeNode * ntree1, viennagrid::quantity_field const & solution1,
                           viennagrid::mesh const & mesh2, NTreeNode * ntree2, viennagrid::quantity_field const & solution2,
                           viennagrid_numeric inside_epsilon, viennagrid_numeric epsilon)
{

  typedef viennagrid::result_of::const_element_range<MeshType>::type          ConstElementRangeType;
  typedef viennagrid::result_of::iterator<ConstElementRangeType>::type        ConstElementRangeIterator;


  ConstElementRangeType vertices1(mesh1, 0);
  ConstElementRangeType vertices2(mesh2, 0);

  ConstElementRangeType lines1(mesh1, 1);
  ConstElementRangeType lines2(mesh2, 1);

  ConstElementRangeType triangles1(mesh1, 2);
  ConstElementRangeType triangles2(mesh2, 2);


  viennagrid_numeric distance_vertices1 = compare_on_vertices(mesh1, ntree1, solution1, mesh2, ntree2, solution2, inside_epsilon);
  viennagrid_numeric distance_vertices2 = compare_on_vertices(mesh2, ntree2, solution2, mesh1, ntree1, solution1, inside_epsilon);


  viennagrid_numeric worst_distance = -1;

  for (ConstElementRangeIterator lit1 = lines1.begin(); lit1 != lines1.end(); ++lit1)
  {
    PointType l1p1 = viennagrid::get_point( viennagrid::vertices(*lit1)[0] );
    PointType l1p2 = viennagrid::get_point( viennagrid::vertices(*lit1)[1] );

    for (ConstElementRangeIterator lit2 = lines2.begin(); lit2 != lines2.end(); ++lit2)
    {
      PointType l2p1 = viennagrid::get_point( viennagrid::vertices(*lit2)[0] );
      PointType l2p2 = viennagrid::get_point( viennagrid::vertices(*lit2)[1] );

      bool intersect = viennagrid::element_line_intersect(*lit2, l1p1, l1p2, epsilon);
      if (intersect)
      {
        std::pair<PointType, PointType> cp = viennagrid::closest_points(*lit1, *lit2);
        if (!(viennagrid::detail::is_equal(epsilon, cp.first, l1p1) ||
              viennagrid::detail::is_equal(epsilon, cp.first, l1p2) ||
              viennagrid::detail::is_equal(epsilon, cp.first, l2p1) ||
              viennagrid::detail::is_equal(epsilon, cp.first, l2p2)))
        {
          viennagrid_numeric value1 = get_value(mesh1, ntree1, solution1, cp.first, inside_epsilon);
          viennagrid_numeric value2 = get_value(mesh2, ntree2, solution2, cp.first, inside_epsilon);

          if ((value1 < 0) || (value2 < 0))
            continue;

          viennagrid_numeric distance = std::abs(value1 - value2);
          if (worst_distance < distance)
            worst_distance = distance;
        }
      }
    }
  }

  return std::max( std::max(distance_vertices1, distance_vertices2), worst_distance );
}









int main(int argc, char **argv)
{
  TCLAP::ValueArg<int> rotation1("","rotation1", "The rotation of mesh 1", false, 0, "int");
  TCLAP::ValueArg<int> rotation2("","rotation2", "The rotation of mesh 2", false, 0, "int");
  TCLAP::SwitchArg dual1("","dual1", "The rotation of mesh 1", false);
  TCLAP::SwitchArg dual2("","dual2", "The rotation of mesh 2", false);

  TCLAP::UnlabeledValueArg<std::string> solution1_filename("solution1", "File name of first solution", true, "", "string");
  TCLAP::UnlabeledValueArg<std::string> solution2_filename("solution2", "File name of first solution", true, "", "string");


  try
  {
    TCLAP::CmdLine cmd("Compare solutions", ' ', "1.0");

    cmd.add( rotation1 );
    cmd.add( rotation2 );
    cmd.add( dual1 );
    cmd.add( dual2 );
    cmd.add( solution1_filename );
    cmd.add( solution2_filename );

    cmd.parse( argc, argv );
  }
  catch (TCLAP::ArgException &e)  // catch any exceptions
  {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return false;
  }


  viennagrid_numeric rotation_angle1 = 2*M_PI / 5 * rotation1.getValue();
  viennagrid_numeric rotation_angle2 = 2*M_PI / 5 * rotation2.getValue();


  viennagrid_mesh cmesh1;
  viennagrid_quantity_field solution1;
  viennagrid_mesh_io reader1;

  viennagrid_mesh_io_create(&reader1);
  viennagrid_error err1 = viennagrid_mesh_io_read_pvd(reader1, solution1_filename.getValue().c_str());
  viennagrid_mesh_io_mesh_get(reader1, &cmesh1);
  viennagrid_mesh_io_quantity_field_get_by_index(reader1, 0, &solution1);


  viennagrid_mesh cmesh2;
  viennagrid_quantity_field solution2;
  viennagrid_mesh_io reader2;

  viennagrid_mesh_io_create(&reader2);
  viennagrid_error err2 = viennagrid_mesh_io_read_pvd(reader2, solution2_filename.getValue().c_str());
  viennagrid_mesh_io_mesh_get(reader2, &cmesh2);
  viennagrid_mesh_io_quantity_field_get_by_index(reader2, 0, &solution2);

  viennagrid::mesh mesh1(cmesh1);
  viennagrid::mesh mesh2(cmesh2);


//   viennagrid_numeric matrix[4];

  Transformation transformation = Transformation::make_rotation_2(- (rotation_angle1-rotation_angle2));
//   Transformation rotation =

  if (dual1.getValue() != dual2.getValue())
  {
    viennagrid_numeric angle = rotation_angle2 - 2*2*M_PI/5;
    viennagrid_numeric back_angle = rotation_angle2 - 3*2*M_PI/5;

    Transformation de_mirror =
      composition(
        Transformation::make_rotation_2(2*M_PI/5),
         Transformation::make_reflection_x(2)
      );

//       composition(
//         Transformation::make_rotation_2(-back_angle),
//         composition(
//           Transformation::make_reflection_x(2),
//           Transformation::make_rotation_2(angle)
//         )
//       );

    transformation = composition(
      transformation,
      composition(
        Transformation::make_rotation_2(M_PI),
        Transformation::make_reflection_2(rotation_angle2 - M_PI/5)
      )
    );
  }

//
//   if (dual.isSet())
//   {
//
//   }
//   else
//   {
//     matrix[0] =   std::cos(rotation_angle);
//     matrix[1] = - std::sin(rotation_angle);
//     matrix[2] =   std::sin(rotation_angle);
//     matrix[3] =   std::cos(rotation_angle);
//   }


  viennagrid_mesh_affine_transform(mesh2.internal(), 2, &transformation.matrix.values[0], 0);

  viennagrid_mesh_io_write(reader1, "mesh1.vtu");
  viennagrid_mesh_io_write(reader2, "mesh2.vtu");


  int max_element_per_node = 10;
  int max_depth = 1000;

  NTreeNode * ntree1 = new NTreeNode(viennagrid::make_point(-1.1, -1.1),
                                     viennagrid::make_point( 1.1,  1.1));
  {
    ElementRangeType cells(mesh1, 2);
    for (ElementRangeIterator cit = cells.begin(); cit != cells.end(); ++cit)
      ntree1->add(*cit, max_element_per_node, max_depth);
  }

  NTreeNode * ntree2 = new NTreeNode(viennagrid::make_point(-1.1, -1.1),
                                     viennagrid::make_point( 1.1,  1.1));
  {
    ElementRangeType cells(mesh2, 2);
    for (ElementRangeIterator cit = cells.begin(); cit != cells.end(); ++cit)
      ntree2->add(*cit, max_element_per_node, max_depth);
  }


  std::cout << compare(mesh1, ntree1, viennagrid::quantity_field(solution1),
                       mesh2, ntree2, viennagrid::quantity_field(solution2), 1e-6, 1e-6) << std::endl;;

  viennagrid_mesh_io_release(reader1);
  viennagrid_mesh_io_release(reader2);

  return 0;
}
