#include "matrix.hpp"


std::vector<double> randomVector(int size)
{
  std::vector<double> result(size);
  for (int i = 0; i != size; ++i)
    result[i] = (double)rand()/(double)(RAND_MAX);
  return result;
}


void prod(CSR const & matrix, std::vector<double> const & x, std::vector<double> & result)
{
  result.resize(matrix.size());
  std::fill(result.begin(), result.end(), 0.0);
  for (std::size_t i = 0; i != matrix.offsets.size()-1; ++i)
  {
    double tmp = 0;
    int j = matrix.offsets[i];
    int ei = matrix.offsets[i+1];

    for (; j != ei; ++j)
      tmp += x[matrix.indices[j]] * matrix.values[j];

    result[i] = tmp;
  }
}

CSR fromMesh(MeshType const & mesh)
{
  CSR result;
  srand (time(NULL));

  result.offsets.push_back(0);

  ElementRangeType vertices(mesh, 0);
  for (ElementRangeIterator vit = vertices.begin(); vit != vertices.end(); ++vit)
  {
    result.offsets.push_back( result.offsets.back() );

    NeighborRangeType neighbor_vertices(mesh, *vit, viennagrid::cell_dimension(mesh), 0);
    for (NeighborRangeIterator nvit = neighbor_vertices.begin(); nvit != neighbor_vertices.end(); ++nvit)
    {
      result.values.push_back( (double)rand()/(double)(RAND_MAX) );
      result.indices.push_back( (*nvit).id().index() );
      result.offsets.back()++;
    }
  }

  return result;
}


void prod(TCSR const & matrix, std::vector<double> const & x, std::vector<double> & result)
{
  result.resize(matrix.size);
  std::fill(result.begin(), result.end(), 0.0);
  for (std::size_t i = 0; i != matrix.templates.size(); ++i)
  {
    for (std::size_t j = 0; j != matrix.templates[i].instances.size(); ++j)
    {
      for (int row = 0; row != matrix.templates[i].matrix.size(); ++row)
      {
        int column_start = matrix.templates[i].matrix.offsets[row];
        int column_end = matrix.templates[i].matrix.offsets[row+1];
        for (int k = column_start; k != column_end; ++k)
        {
          int column = matrix.templates[i].matrix.indices[k];
          double value = matrix.templates[i].matrix.values[k];

          int global_row = matrix.templates[i].instances[j].local_to_global_mapping[row];
          int global_column = matrix.templates[i].instances[j].local_to_global_mapping[column];

          if ((global_row != -1) && (global_column != -1))
            result[global_row] += value * x[global_column];
        }
      }
    }
  }
}

TCSR fromMesh(boost::shared_ptr<TemplatedMesh> const & tm)
{
  int global_vertex_count = viennagrid::elements(tm->structure_instance, 0).size();

  TCSR result;

  result.templates.resize( tm->templates.size() );
  for (std::size_t i = 0; i != result.templates.size(); ++i)
  {
    result.templates[i].matrix = fromMesh(tm->templates[i]);
  }

  for (std::size_t i = 0; i != tm->instances.size(); ++i)
  {
    int template_index = tm->instances[i].mesh_template_index;
    result.templates[template_index].instances.push_back( TCSRInstance(result.templates[template_index].matrix.size(), global_vertex_count) );
  }

  result.size = global_vertex_count;

  return result;
}


