#ifndef DISS_COMMON_HPP_
#define DISS_COMMON_HPP_

#include <string>

#include "viennagrid/algorithm/inclusion.hpp"
#include "viennagrid/algorithm/geometry.hpp"
#include "viennagrid/algorithm/geometric_transform.hpp"
#include "viennagrid/io/vtk_writer.hpp"

#include "viennameshpp/core.hpp"
#include "viennameshpp/algorithm_pipeline.hpp"

#include "statistics/element_metrics.hpp"
#include "statistics/statistic.hpp"

#include "multimesh_element_copy_map.hpp"


std::vector<std::string> split_string( std::string const & str, std::string const & delimiter );
std::vector<std::string> split_string_brackets( std::string const & str, std::string const & delimiter );

template<typename T>
std::vector<T> split(std::string const & str, std::string const & delimiter = ";")
{
  std::vector<std::string> tokens = split_string(str, delimiter);
  std::vector<T> result(tokens.size());

  for (std::size_t i = 0; i != tokens.size(); ++i)
  {
    result[i] = boost::lexical_cast<T>(tokens[i]);
  }

  return result;
}



typedef viennagrid::mesh                                              MeshType;
typedef viennagrid::result_of::region<MeshType>::type                 RegionType;

typedef viennagrid::result_of::element<MeshType>::type                ElementType;
typedef viennagrid::result_of::point<MeshType>::type                  PointType;

typedef viennagrid::result_of::element_range<MeshType>::type          ElementRangeType;
typedef viennagrid::result_of::iterator<ElementRangeType>::type       ElementRangeIterator;

typedef viennagrid::result_of::element_range<ElementType>::type           BoundaryElementRangeType;
typedef viennagrid::result_of::iterator<BoundaryElementRangeType>::type   BoundaryElementRangeIterator;

typedef viennagrid::result_of::vertex_range<MeshType>::type           VertexRangeType;
typedef viennagrid::result_of::iterator<VertexRangeType>::type        VertexRangeIterator;

typedef viennagrid::result_of::cell_range<MeshType>::type             CellRangeType;
typedef viennagrid::result_of::iterator<CellRangeType>::type          CellRangeIterator;

typedef viennagrid::result_of::coboundary_range<MeshType>::type       CoboundaryRangeType;
typedef viennagrid::result_of::iterator<CoboundaryRangeType>::type    CoboundaryRangeIterator;

typedef viennagrid::result_of::neighbor_range<MeshType>::type         NeighborRangeType;
typedef viennagrid::result_of::iterator<NeighborRangeType>::type      NeighborRangeIterator;

typedef viennagrid::result_of::region_range<MeshType>::type           RegionRangeType;
typedef viennagrid::result_of::iterator<RegionRangeType>::type        RegionRangeIterator;


typedef viennagrid::quantity_field QuantityField;


typedef viennagrid::multimesh_element_copy_map    CopyMapType;
// typedef viennagrid::result_of::element_copy_map<>::type    CopyMapType;


using boost::shared_ptr;





typedef viennamesh::statistic<viennagrid_numeric>         StatisticType;
typedef StatisticType::histogram_type                     HistogramType;



viennamesh::context_handle context();









template<typename Iterator>
PointType centroid(Iterator begin, Iterator const & end)
{
  if (begin == end)
    return PointType();

  PointType centroid = *begin; ++begin;
  int count = 1;
  for (; begin != end; ++begin, ++count)
    centroid += *begin;

  return centroid / count;
}

inline viennagrid_numeric angle(PointType vec)
{
  vec.normalize();
  assert(vec.size() == 2);
  return atan2(vec[1], vec[0]);
}

inline viennagrid_numeric determinant(PointType const & x, PointType const & y, PointType const & z)
{
  return x[0]*y[1]*z[2] + y[0]*z[1]*x[2] + z[0]*x[1]*y[2] - x[2]*y[1]*z[0] - y[2]*z[1]*x[0] - z[2]*x[1]*y[0];
}

struct MatrixType
{
  MatrixType() : row_count(0), column_count(0) {}
  MatrixType(int row_count_, int column_count_) : row_count(row_count_), column_count(column_count_)
  {
    values.resize(row_count * column_count);
  }

  std::vector<viennagrid_numeric> values;

  viennagrid_numeric operator()(int row, int col) const { return values[column_count*row+col]; }
  viennagrid_numeric & operator()(int row, int col) { return values[column_count*row+col]; }

  viennagrid_numeric * ptr() { return &values[0]; }
  viennagrid_numeric const * ptr() const { return &values[0]; }

  bool valid() const { return !values.empty(); }

  int row_count;
  int column_count;
};

void prod(MatrixType const & matrix, viennagrid_numeric const * pt, viennagrid_numeric * result);
PointType prod(MatrixType const & matrix, PointType const & pt);
MatrixType prod(MatrixType const & A1, MatrixType const & A2);



struct Transformation;

Transformation composition(Transformation const & T1, Transformation const & T2);

struct Transformation
{
  static Transformation make_identity(int geo_dim)
  {
    Transformation result;
    result.matrix.row_count = geo_dim;
    result.matrix.column_count = geo_dim;
    result.matrix.values.resize(geo_dim*geo_dim, 0);
    for (int i = 0; i != geo_dim; ++i)
      result.matrix(i,i) = 1;
    result.translate.resize(geo_dim, 0);
    return result;
  }

  static Transformation make_translation(PointType const & translate_)
  {
    Transformation result = make_identity(translate_.size());
    result.translate = translate_;
    return result;
  }

  static Transformation make_reflection(int geo_dim, int dim)
  {
    Transformation result = make_identity(geo_dim);
    result.matrix(dim,dim) = -1;
    return result;
  }

  static Transformation make_reflection_x(int geo_dim)
  {
    return make_reflection(geo_dim, 0);
  }

  static Transformation make_reflection_y(int geo_dim)
  {
    return make_reflection(geo_dim, 1);
  }

  static Transformation make_reflection_xy(int geo_dim)
  {
    Transformation result = make_identity(geo_dim);
    result.matrix(0,0) = -1;
    result.matrix(1,1) = -1;
    return result;
  }

  static Transformation make_reflection_2(PointType const & p)
  {
    // https://en.wikipedia.org/wiki/Transformation_matrix#Reflection
    Transformation result = make_identity(2);
    result.matrix(0,0) = p[0]*p[0] - p[1]*p[1];
    result.matrix(1,0) = 2*p[0]*p[1];
    result.matrix(0,1) = 2*p[0]*p[1];
    result.matrix(1,1) = p[1]*p[1] - p[0]*p[0];
    return result;
  }

  static Transformation make_reflection_2(double angle)
  {
    // http://planetmath.org/derivationof2dreflectionmatrix
    Transformation result = make_identity(2);
    result.matrix(0,0) = std::cos(2*angle);
    result.matrix(1,0) = std::sin(2*angle);
    result.matrix(0,1) = std::sin(2*angle);
    result.matrix(1,1) = -std::cos(2*angle);
    return result;
  }

  static Transformation make_rotation_2(double angle)
  {
    Transformation result = make_identity(2);

    result.matrix(0,0) =  std::cos(angle);
    result.matrix(1,0) = -std::sin(angle);
    result.matrix(0,1) =  std::sin(angle);
    result.matrix(1,1) =  std::cos(angle);

    return result;
  }

  static Transformation make_rotation_3z(double angle)
  {
    Transformation result = make_identity(3);

    result.matrix(0,0) =  std::cos(angle);
    result.matrix(1,0) = -std::sin(angle);
    result.matrix(0,1) =  std::sin(angle);
    result.matrix(1,1) =  std::cos(angle);

    return result;
  }

  static std::pair<Transformation, Transformation> make_facet_projection(viennagrid_plc plc, viennagrid_element_id facet_id, viennagrid_numeric eps);

  PointType operator()(PointType const & pt) const
  {
    PointType result = prod(matrix, pt);
    return result+translate;
  }

  void operator()(viennagrid_numeric const * pt, viennagrid_numeric * result) const
  {
    prod(matrix, pt, result);
    viennagrid_point_add(translate.size(), result, &translate[0], result);
  }

  MeshType operator()(MeshType const & mesh) const
  {
    MeshType result;
    viennagrid::copy(mesh, result);
    viennagrid::affine_transform(result, &matrix.values[0], translate);
    return result;
  }


  bool valid() const { return matrix.valid(); }

  void print() const
  {
    if (!valid())
    {
      std::cout << "Invalid transformation" << std::endl;
    }
    else
    {
      for (int i = 0; i != matrix.row_count; ++i)
      {
        std::cout << "[";
        for (int j = 0; j != matrix.column_count; ++j)
        {
          std::cout << matrix(i,j);
          if (j != matrix.column_count-1)
            std::cout << ",";
        }
        std::cout << "]   (" << translate[i] << ")" << std::endl;
      }
    }
  }





  bool test_transform(std::vector<PointType> const & pts_from,
                      std::vector<PointType> const & pts_to,
                      viennagrid_numeric eps)
  {
    if (pts_from.size() != pts_to.size())
      return false;

    std::vector<bool> used(pts_from.size(), false);

    for (std::size_t i = 0; i != pts_from.size(); ++i)
    {
      PointType transformed = (*this)( pts_from[i] );

      bool found = false;
      for (std::size_t j = 0; j != pts_to.size(); ++j)
      {
        if (used[j])
          continue;

        if ( viennagrid::detail::is_equal(eps, transformed, pts_to[j]) )
        {
          used[j] = true;
          found = true;
          break;
        }
      }

      if (!found)
        return false;
    }

    return true;
  }


  static Transformation find_2d_rigid(std::vector<PointType> pts_from,
                                      std::vector<PointType> const & pts_to,
                                      viennagrid_numeric eps,
                                      bool swap_orientation = false)
  {
    Transformation result;

//     std::vector<PointType> tmp_from = pts_from;
    Transformation base = Transformation::make_identity(2);
    if (swap_orientation)
    {
      base = Transformation::make_reflection_x(2);

      for (std::size_t i = 0; i != pts_from.size(); ++i)
        pts_from[i][0] = -pts_from[i][0];
    }


    PointType centroid_from = centroid(pts_from.begin(), pts_from.end());
    PointType centroid_to = centroid(pts_to.begin(), pts_to.end());

    viennagrid_numeric radius0 = viennagrid::norm_2( pts_from[0]-centroid_from );

    std::vector<viennagrid_numeric> angles;
    for (std::size_t i = 0; i != pts_to.size(); ++i)
    {
      viennagrid_numeric tmp = viennagrid::norm_2( pts_to[i]-centroid_to );
      if ( std::abs(radius0-tmp) < eps )
      {
        PointType v_from = pts_from[0]-centroid_from;
        PointType v_to = pts_to[i]-centroid_to;

        angles.push_back( angle(v_from) - angle(v_to) );
      }
    }

    for (std::size_t i = 0; i != angles.size(); ++i)
    {
      result = composition(
                   composition(
                       Transformation::make_translation(centroid_to),
                       composition(
                           Transformation::make_rotation_2(angles[i]),
                           Transformation::make_translation(-centroid_from)
                       )
                   ),
                   base
               );

      if ( result.test_transform(pts_from, pts_to, eps) )
      {
  //       std::cout << "Transformation found!!" << std::endl;
        return result;
      }
    }

  //   std::cout << "No transformation found!!" << std::endl;
  //   assert(false);
    return result;
  }



  MatrixType matrix;
  PointType translate;
};



inline Transformation composition(Transformation const & T1, Transformation const & T2)
{
  Transformation result;
  result.matrix = prod(T1.matrix, T2.matrix);
  result.translate = prod(T1.matrix, T2.translate) + T1.translate;
  return result;
}


inline void rotate_2(viennagrid_numeric * pt, double angle)
{
  double ca = std::cos(angle);
  double sa = std::sin(angle);

  viennagrid_numeric tmp[2];
  tmp[0] = ca*pt[0] - sa*pt[1];
  tmp[1] = sa*pt[0] + ca*pt[1];

  pt[0] = tmp[0];
  pt[1] = tmp[1];
}

inline PointType rotate(PointType const & pt, double angle)
{
  PointType tmp = pt;
  rotate_2(&tmp[0], angle);
  return tmp;
}


inline PointType rotate_z(PointType const & pt, double angle)
{
  PointType tmp = pt;
  rotate_2(&tmp[0], angle);
  return tmp;
}






inline void print_point(viennagrid_dimension geo_dim, viennagrid_numeric const * pt)
{
  std::cout << "(";
  for (viennagrid_dimension i = 0; i != geo_dim; ++i)
  {
    std::cout << pt[+i];
    if (i != geo_dim-1)
      std::cout << ",";
  }
  std::cout << ")";
}

inline PointType get_point(viennagrid_plc plc, viennagrid_element_id vid)
{
  viennagrid_dimension geo_dim;
  viennagrid_plc_geometric_dimension_get(plc, &geo_dim);

  viennagrid_numeric * tmp;
  viennagrid_plc_vertex_coords_get(plc, vid, &tmp);

  PointType result(geo_dim);
  std::copy(tmp, tmp+geo_dim, result.begin());
  return result;
}

inline viennagrid_element_id make_vertex(viennagrid_plc plc, PointType const & pt)
{
  viennagrid_element_id vtx;
  viennagrid_plc_vertex_create(plc, &pt[0], &vtx);
  return vtx;
}

inline viennagrid_element_id make_line(viennagrid_plc plc, viennagrid_element_id vtx0, viennagrid_element_id vtx1)
{
  viennagrid_element_id line;
  viennagrid_plc_line_create(plc, vtx0, vtx1, &line);
  return line;
}

inline viennagrid_plc make_plc(viennagrid_dimension geo_dim)
{
  viennagrid_plc plc;
  viennagrid_plc_create(&plc);
  viennagrid_plc_geometric_dimension_set(plc, geo_dim);
  return plc;
}



inline void add_seed_points(viennagrid_plc plc, viennamesh::seed_point_container const & seed_points)
{
  for (std::size_t i = 0; i != seed_points.size(); ++i)
    viennagrid_plc_seed_point_add(plc, &seed_points[i].first[0], seed_points[i].second);
}




MeshType plc_to_linemesh(viennagrid_plc plc);
inline void write_mesh(MeshType const & mesh, std::string const & filename)
{
  viennagrid::io::vtk_writer<MeshType> writer;
  writer(mesh, filename);
}




StatisticType boundary_min_angle_statistic(MeshType const & mesh);

QuantityField min_angles(MeshType const & mesh);
QuantityField radius_edge_ratio(MeshType const & mesh);
StatisticType statistics_from_quantity_field(MeshType const & mesh, QuantityField const & qf,
                                             viennagrid_numeric min, viennagrid_numeric max, std::size_t count);




int mesh_size(MeshType const & mesh);
int mr_mesh_size(MeshType const & mesh);



static const int sizeof_int             = 4;
static const int sizeof_ptr             = 8;
static const int sizeof_numeric         = 8;
static const int sizeof_element_type    = 1;
static const int sizeof_index_type      = 4;
static const int sizeof_region_id_type  = 2;


#endif
