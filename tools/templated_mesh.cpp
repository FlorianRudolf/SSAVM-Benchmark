#include "templated_mesh.hpp"

MeshType const & MeshInstance::mesh_template() const
{
  return templated_mesh->templates[mesh_template_index];
}

Transformation const & MeshInstance::transformation() const
{
  return templated_mesh->transformations[transformation_index];
}

int MeshInstance::global_vertex_count() const
{
  int count = 0;

  for (std::map< ElementType, std::vector<ElementType> >::iterator it = templated_mesh->global_to_local_mapping.begin();
       it != templated_mesh->global_to_local_mapping.end(); ++it)
  {
    std::vector<ElementType> const & tmp = (*it).second;
    if (tmp.size() > 1)
    {
      for (std::size_t i = 0; i != tmp.size(); ++i)
      {
        if (tmp[i].get_mesh() == instance_mesh)
        {
          ++count;
        }
      }
    }
  }

  return count;
}

