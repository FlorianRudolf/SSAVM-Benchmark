#ifndef TEMPLATED_MESH_HPP_
#define TEMPLATED_MESH_HPP_

#include "common.hpp"
#include <boost/enable_shared_from_this.hpp>


struct TemplatedMesh;

struct MeshInstance
{
  MeshInstance(shared_ptr<TemplatedMesh> const & templated_mesh_) : templated_mesh(templated_mesh_) {}

  MeshType instance_mesh;
  int mesh_template_index;
  int transformation_index;
  int region_id;

  shared_ptr<TemplatedMesh> templated_mesh;

  MeshType const & mesh_template() const;
  Transformation const & transformation() const;

  void setup()
  {
    instance_mesh = transformation()(mesh_template());
  }

  int cell_count() const
  {
    return viennagrid::cells(mesh_template()).size();
  }

  int global_vertex_count() const;
};

struct TemplatedMesh : public boost::enable_shared_from_this<TemplatedMesh>
{
  TemplatedMesh(PointType const & min_, PointType const & max_) : min(min_), max(max_) {}

  std::vector<MeshType> templates;
  std::vector<Transformation> transformations;
  std::vector<MeshInstance> instances;

  PointType min;
  PointType max;

  int geo_dim() const { return min.size(); }

  MeshType structure_instance;
  std::map< ElementType, std::vector<ElementType> > global_to_local_mapping; // maps a vertex in the structure instance to the instance vertex
  int global_shared_vertex_count;

  QuantityField template_index_field;
  QuantityField instance_index_field;


  int instance_count(int mesh_template_index) const
  {
    int count = 0;
    for (std::size_t i = 0; i != instances.size(); ++i)
    {
      if (instances[i].mesh_template_index == mesh_template_index)
        ++count;
    }
    return count;
  }


  int add_template(MeshType const & mesh_template)
  {
    templates.push_back(mesh_template);
    return templates.size()-1;
  }


  void add_instance(int mesh_template_index, int transformation_index, int region_id)
  {
    MeshInstance mi( shared_from_this() );
    mi.mesh_template_index = mesh_template_index;
    mi.transformation_index = transformation_index;
    mi.region_id = region_id;
    instances.push_back(mi);
  }

  void add_instance(int mesh_template_index, Transformation const & transformation, int region_id)
  {
    transformations.push_back(transformation);

    MeshInstance mi( shared_from_this() );
    mi.mesh_template_index = mesh_template_index;
    mi.transformation_index = transformations.size()-1;
    mi.region_id = region_id;
    instances.push_back(mi);
  }


  void add_trivial_instance(MeshType const & mesh_template,
                            int identity_index = -1)
  {
    viennagrid_dimension geo_dim = viennagrid::geometric_dimension(mesh_template);

    if (identity_index <= 0)
    {
      identity_index = transformations.size();
      transformations.push_back( Transformation::make_identity(geo_dim) );
    }

    RegionRangeType regions(mesh_template);
    if (regions.empty())
    {
      int template_index = add_template(mesh_template);
      add_instance(template_index, identity_index, 0);
    }
    else
    {
      int index = 0;
      for (RegionRangeIterator rit = regions.begin(); rit != regions.end(); ++rit, ++index)
      {
        MeshType tmp;
        viennagrid::copy(*rit, tmp);

        int template_index = add_template(tmp);
        add_instance(template_index, identity_index, (*rit).id());
      }
    }
  }



  void add_instances_rotatation_z(int rotational_symmetry_order,
                                  MeshType const & mesh_template,
                                  std::vector<PointType> const & centers,
                                  int rotation_index_start = -1)
  {
    double angle = 2*M_PI / rotational_symmetry_order;

    if (rotation_index_start < 0)
    {
      rotation_index_start = transformations.size();
      for (std::size_t i = 0; i != centers.size(); ++i)
      {
        for (int j = 0; j != rotational_symmetry_order; ++j)
          transformations.push_back( composition(Transformation::make_translation(centers[i]),
                                                 Transformation::make_rotation_3z(angle*j)) );
      }
    }

    RegionRangeType regions(mesh_template);
    if (regions.empty() || (regions.size() == 1))
    {
      int template_index = add_template(mesh_template);
      for (std::size_t i = 0; i != centers.size(); ++i)
        for (int j = 0; j != rotational_symmetry_order; ++j)
          add_instance(template_index, rotation_index_start+rotational_symmetry_order*i+j, 0);
    }
    else
    {
      int index = 0;
      for (RegionRangeIterator rit = regions.begin(); rit != regions.end(); ++rit, ++index)
      {
        MeshType tmp;
        viennagrid::copy(*rit, tmp);
        int template_index = add_template(tmp);

        for (std::size_t i = 0; i != centers.size(); ++i)
          for (int j = 0; j != rotational_symmetry_order; ++j)
            add_instance(template_index, rotation_index_start+rotational_symmetry_order*i+j, (*rit).id());
      }
    }
  }



  void add_instances_rotatation_z(int rotational_symmetry_order,
                                  MeshType const & mesh_template,
                                  int rotation_index_start = -1)
  {
    double angle = 2*M_PI / rotational_symmetry_order;

    if (rotation_index_start < 0)
    {
      rotation_index_start = transformations.size();
      for (int i = 0; i != rotational_symmetry_order; ++i)
        transformations.push_back( Transformation::make_rotation_3z(angle*i) );
    }

    RegionRangeType regions(mesh_template);
    if (regions.empty() || (regions.size() == 1))
    {
      int template_index = add_template(mesh_template);
      for (int i = 0; i != rotational_symmetry_order; ++i)
        add_instance(template_index, rotation_index_start+i, 0);
    }
    else
    {
      int index = 0;
      for (RegionRangeIterator rit = regions.begin(); rit != regions.end(); ++rit, ++index)
      {
        MeshType tmp;
        viennagrid::copy(*rit, tmp);
        int template_index = add_template(tmp);

        for (int i = 0; i != rotational_symmetry_order; ++i)
          add_instance(template_index, rotation_index_start+i, (*rit).id());
      }
    }
  }


  void add_instances_rotatation(int rotational_symmetry_order,
                                MeshType const & mesh_template,
                                int rotation_index_start = -1)
  {
    double angle = 2*M_PI / rotational_symmetry_order;

    if (rotation_index_start < 0)
    {
      rotation_index_start = transformations.size();
      for (int i = 0; i != rotational_symmetry_order; ++i)
        transformations.push_back( Transformation::make_rotation_2(angle*i) );
    }

    RegionRangeType regions(mesh_template);
    if (regions.empty() || (regions.size() == 1))
    {
      int template_index = add_template(mesh_template);
      for (int i = 0; i != rotational_symmetry_order; ++i)
        add_instance(template_index, rotation_index_start+i, 0);
    }
    else
    {
      int index = 0;
      for (RegionRangeIterator rit = regions.begin(); rit != regions.end(); ++rit, ++index)
      {
        MeshType tmp;
        viennagrid::copy(*rit, tmp);
        int template_index = add_template(tmp);

        for (int i = 0; i != rotational_symmetry_order; ++i)
          add_instance(template_index, rotation_index_start+i, (*rit).id());
      }
    }
  }



  void setup(viennagrid_numeric eps)
  {
    viennagrid::clear(structure_instance);
    CopyMapType cm(structure_instance, eps, min, max);

    for (std::size_t i = 0; i != instances.size(); ++i)
    {
      instances[i].setup();
    }

    for (std::size_t i = 0; i != instances.size(); ++i)
    {
      RegionType region = structure_instance.get_or_create_region(instances[i].region_id);

      CellRangeType cells(instances[i].instance_mesh);
      for (CellRangeIterator cit = cells.begin(); cit != cells.end(); ++cit)
      {
        ElementType nc = cm(*cit);
        viennagrid::add(region, nc);

        template_index_field.set(nc, instances[i].mesh_template_index);
        instance_index_field.set(nc, i);
//         std::cout << i << " " << instances[i].mesh_template_index << std::endl;
      }
    }

    for (std::size_t i = 0; i != instances.size(); ++i)
    {
      VertexRangeType vertices(instances[i].instance_mesh);
      for (VertexRangeIterator vit = vertices.begin(); vit != vertices.end(); ++vit)
      {
        ElementType sivtx = cm.find_vertex(*vit);
        global_to_local_mapping[sivtx].push_back(*vit);
      }
    }

    global_shared_vertex_count = 0;
    for (std::map< ElementType, std::vector<ElementType> >::iterator it = global_to_local_mapping.begin();
                                                                     it != global_to_local_mapping.end();
                                                                   ++it)
    {
      if ((*it).second.size() != 1)
      {
        ++global_shared_vertex_count;
      }
    }
  }



  int structure_instance_cell_count() const
  {
    int cell_count = 0;

    for (std::size_t i = 0; i != instances.size(); ++i)
    {
      cell_count += instances[i].cell_count();
    }

    return cell_count;
  }

  int size_transformation_functions() const
  {
    int tmp = 0;

    tmp += sizeof_ptr + sizeof_numeric * geo_dim() * geo_dim();     // matrix
    tmp += sizeof_ptr + sizeof_numeric * geo_dim();                 // translation

    return tmp * transformations.size();
  }




  int structure_instance_size() const
  {
    return mesh_size(structure_instance);
  }

  int structure_instance_mr_size() const
  {
    return mr_mesh_size(structure_instance);
  }



  int size() const
  {
    int tmp = 0;

    tmp += sizeof_int + sizeof_ptr + size_transformation_functions(); // transformation_count and transformations

    tmp += sizeof_int + sizeof_ptr;    // template_count and templates
    for (std::size_t i = 0; i != templates.size(); ++i)
    {
      tmp += mesh_size( templates[i] );                               // mesh
      tmp += sizeof_int;                                                // instance_count
      tmp += sizeof_ptr + instance_count(i) * sizeof_int;                  // transformation_indices
      tmp += sizeof_ptr + instance_count(i) * sizeof_region_id_type;       // region_ids
    }

    return tmp;
  }


  int size_SVB() const
  {
    int tmp = 0;

    tmp += sizeof_int + sizeof_ptr + size_transformation_functions(); // transformation_count and transformations

    tmp += sizeof_int + sizeof_ptr;    // template_count and templates
    for (std::size_t i = 0; i != templates.size(); ++i)
    {
      tmp += mesh_size( templates[i] );   // mesh

      tmp += sizeof_int + sizeof_ptr; // instance_count and instances
      for (std::size_t j = 0; j != instances.size(); ++j)
      {
        if (instances[j].mesh_template_index != i)
          continue;

        tmp += sizeof_region_id_type + sizeof_int; // region_id and transformation_index
        tmp += sizeof_int + sizeof_ptr; // local_to_global_vertex_mapping_count and local_to_global_vertex_mapping
        tmp += instances[i].global_vertex_count() * 2 * sizeof_int;
      }
    }

    tmp += sizeof_int + sizeof_ptr; // global_vertex_count and global_vertices
    tmp += sizeof_numeric * geo_dim() * global_shared_vertex_count;

    return tmp;
  }


};

#endif
