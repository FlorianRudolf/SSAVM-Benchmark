#ifndef DISS_MULTIMESH_ELEMENT_COPY_MAP_HPP_
#define DISS_MULTIMESH_ELEMENT_COPY_MAP_HPP_

#include "viennagrid/viennagrid.hpp"

namespace viennagrid
{


  inline bool bin(int val, int index)
  {
    return ((val >> index) & 0b1) == 0b1;
  }

  template<typename PointType>
  PointType binary_coord(PointType const & pt, int mask)
  {
    PointType result = pt;
    for (std::size_t i = 0; i != pt.size(); ++i)
    {
      if (!bin(mask, i))
        result[i] = 0;
    }

    return result;
  }


  class nd_vertex_tree_node
  {
  public:
    typedef viennagrid::mesh MeshType;
    typedef viennagrid::result_of::point<MeshType>::type PointType;
    typedef viennagrid::result_of::element<MeshType>::type ElementType;

    nd_vertex_tree_node(PointType const & min_, PointType const & max_, viennagrid_numeric eps_) : min(min_), max(max_), eps(eps_)
    {
      assert( min.size() == max.size() );
      assert( min < max );
    }

    std::size_t geo_dim() const { return min.size(); }

    nd_vertex_tree_node * find_nodes(PointType const & pt,
                                     viennagrid_numeric eps,
                                     std::vector<nd_vertex_tree_node*> & nodes)
    {
      if (does_intersect(pt, eps))
      {
        if (children.empty())
        {
          nodes.push_back(this);
          return is_inside(pt) ? this : NULL;
        }
        else
        {
          nd_vertex_tree_node * result = NULL;

          for (std::size_t i = 0; i != children.size(); ++i)
          {
            nd_vertex_tree_node * tmp = children[i]->find_nodes(pt, eps, nodes);
            if (tmp)
              result = tmp;
          }

          return result;
        }
      }
      else
        return NULL;
    }

    nd_vertex_tree_node * find_node(PointType const & pt)
    {
      if (is_inside(pt))
      {
        if (children.empty())
          return this;
        else
        {
          for (std::size_t i = 0; i != children.size(); ++i)
          {
            nd_vertex_tree_node * tmpnode = children[i]->find_node(pt);
            if (tmpnode)
              return tmpnode;
          }

          std::cout << "!!!!!!!" << std::endl;
          assert(false);
          return NULL;
        }
      }
      else
        return NULL;
    }

    ElementType find_vertex(PointType const & pt, viennagrid_numeric eps) const
    {
      for (std::size_t i = 0; i != vertices.size(); ++i)
      {
        if ( viennagrid::detail::is_equal(eps, pt, viennagrid::get_point(vertices[i])) )
          return vertices[i];
      }

      return ElementType();
    }

    void add_vertex(ElementType const & element,
                    int max_vertices_per_node)
    {
      PointType pt = viennagrid::get_point(element);
      nd_vertex_tree_node * node = find_node(pt);
      node->add_vertex_impl(pt, element, max_vertices_per_node);
    }

    bool is_inside(PointType const & pt) const
    {
      return (min <= pt) && (pt < max);
    }

    bool does_intersect(PointType const & pt, viennagrid_numeric radius) const
    {
      viennagrid_numeric radius_squared = radius * radius;

      for (std::size_t i = 0; i != pt.size(); ++i)
      {
        if (pt[i] < min[i])
          radius_squared -= (pt[i]-min[i])*(pt[i]-min[i]);
        else if (pt[i] > max[i])
          radius_squared -= (pt[i]-max[i])*(pt[i]-max[i]);
      }

      return radius_squared > 0;
    }

    std::size_t vertex_count() const
    {
      if (children.empty())
        return vertices.size();

      std::size_t count = 0;
      for (std::size_t i = 0; i != children.size(); ++i)
        count += children[i]->vertex_count();

      return count;
    }

  private:

    void add_vertex_impl(PointType const & pt,
                         ElementType const & element,
                         int max_vertices_per_node)
    {
      if (!is_inside(pt))
        return;

      if (vertices.size() < max_vertices_per_node)
      {
        vertices.push_back(element);
      }
      else
      {
//         std::cout << "min = " << min << "  max = " << max << "   (count = " << vertices.size() << ")" << std::endl;
        for (int i = 0; i != (1 << geo_dim()); ++i)
        {
          boost::shared_ptr<nd_vertex_tree_node> child( new nd_vertex_tree_node(lower(i), upper(i), eps) );

          for (std::size_t j = 0; j != vertices.size(); ++j)
            child->add_vertex_impl( viennagrid::get_point(vertices[j]), vertices[j], max_vertices_per_node );

//           std::cout << "  new min = " << lower(i) << "  max = " << upper(i) << "   (count = " << child->vertices.size() << ")" << std::endl;

          children.push_back( child );
        }

        vertices.clear();

        add_vertex(element, max_vertices_per_node);
      }
    }

    PointType lower(int mask) const
    {
      return min + binary_coord(size(), mask) * 0.5;
    }

    PointType upper(int mask) const
    {
      return (min+max) * 0.5 + binary_coord(size(), mask) * 0.5;
    }

    PointType size() const
    {
      return max-min;
    }

    std::vector< boost::shared_ptr<nd_vertex_tree_node> > children;
    PointType min;
    PointType max;
    viennagrid_numeric eps;

    std::vector<ElementType> vertices;
  };






  /** @brief A helper class for element copy operation between two differen meshes.
    *
    * @tparam SrcMeshT      The mesh type of the source mesh
    * @tparam DstMeshT      The mesh type of the destination mesh
    */
  class multimesh_element_copy_map
  {
  public:

    typedef viennagrid::const_mesh SrcMeshType;
    typedef viennagrid::mesh DstMeshType;

    typedef viennagrid::result_of::point<SrcMeshType>::type PointType;

    typedef viennagrid::result_of::coord<DstMeshType>::type DstNumericType;
    typedef viennagrid::result_of::element<SrcMeshType>::type SrcElementType;
    typedef viennagrid::result_of::element<DstMeshType>::type DstElementType;
    typedef viennagrid::result_of::element_id<SrcMeshType>::type SrcElementIDType;

    /** @brief The constructor, requires the destination mesh where the elements are copied to.
      *
      * @param  dst_mesh_                The destination mesh
      */
//     multimesh_element_copy_map(DstMeshType const & dst_mesh_in, PointType const & min, PointType const & max) : dst_mesh_(dst_mesh_in), eps(-1.0), ndtree(new nd_vertex_tree_node(min,max)) {}
    multimesh_element_copy_map(DstMeshType const & dst_mesh_in, viennagrid_numeric eps_in, PointType const & min, PointType const & max) : dst_mesh_(dst_mesh_in), eps(eps_in), ndtree(new nd_vertex_tree_node(min,max,eps_in)) {}

    /** @brief Copies one vertex to the destination mesh. If the vertex is already present in the destination mesh, the vertex handle of this vertex is return, otherwise a new vertex is created in the destination mesh.
      *
      * @param  src_vertex              The vertex to be copied
      */
    template<bool element_is_const>
    DstElementType operator()( base_element<element_is_const> const & src )
    {
      assert( src.get_mesh() != dst_mesh() );

      if (src.is_vertex())
        return copy_vertex(src);
      return copy_element(src);
    }


    void operator()(SrcMeshType const & src_mesh)
    {
      typedef viennagrid::result_of::const_cell_range<SrcMeshType>::type ConstCellRangeType;
      typedef viennagrid::result_of::iterator<ConstCellRangeType>::type ConstCellRangeIterator;

      ConstCellRangeType cells(src_mesh);
      for (ConstCellRangeIterator cit = cells.begin(); cit != cells.end(); ++cit)
        (*this)(*cit);
    }




    template<bool element_is_const>
    DstElementType find_vertex(base_element<element_is_const> const & src)
    {
      assert( src.is_vertex() );

      typename std::map<SrcElementType, DstElementType>::iterator vit = vertex_map.find( src );
      if (vit != vertex_map.end())
        return vit->second;
      else
      {
        std::cout << "ERROR!!" << std::endl;
      }

      return DstElementType();
    }


    /** @brief Copies a whole element including its vertices to the destination mesh.
      *
      * @param  el                      The element to be copied
      */
    template<bool element_is_const>
    DstElementType copy_vertex( base_element<element_is_const> const & src )
    {
      typename std::map<SrcElementType, DstElementType>::iterator vit = vertex_map.find( src );
      if (vit != vertex_map.end())
        return vit->second;
      else
      {
        PointType src_pt = viennagrid::get_point(src);
        std::vector<nd_vertex_tree_node*> nodes;
        nd_vertex_tree_node * node = ndtree->find_nodes(src_pt, eps, nodes);

        DstElementType vtx;
        for (std::size_t i = 0; i != nodes.size(); ++i)
        {
          vtx = nodes[i]->find_vertex(src_pt, eps);
          if (vtx.valid())
            break;
        }

        if (!vtx.valid())
        {
          vtx = viennagrid::make_vertex(dst_mesh(), src_pt);
          node->add_vertex(vtx, 10);
        }


//         PointType src_pt = viennagrid::get_point(src);
//         nd_vertex_tree_node * node = ndtree->find_node(src_pt);
//         DstElementType vtx = node->find_vertex(src_pt, eps);
//         if (!vtx.valid())
//         {
//           vtx = viennagrid::make_vertex(dst_mesh(), src_pt);
//           node->add_vertex(vtx, 100);
//         }
//         else
//         {
//           DstElementType tmp = viennagrid::make_unique_vertex(dst_mesh(), viennagrid::get_point(src), eps);
//           if (tmp != vtx)
//           {
//             std::cout << "ERROR !!!!!!!!!!!!!" << std::endl;
//             std::cout << "     " << vtx << std::endl;
//             std::cout << "     " << tmp << std::endl;
//           }
//         }


//         DstElementType vtx = viennagrid::make_unique_vertex(dst_mesh(), viennagrid::get_point(src), eps);

        vertex_map[src] = vtx;
        return vtx;
      }
    }


    /** @brief Copies a whole element including its vertices to the destination mesh.
      *
      * @param  el                      The element to be copied
      */
    template<bool element_is_const>
    DstElementType copy_element( base_element<element_is_const> const & src )
    {
      typedef base_element<element_is_const> SrcElementType;
      typedef typename viennagrid::result_of::const_element_range<SrcElementType>::type ConstVerticesOnElementRangeType;
      typedef typename viennagrid::result_of::iterator<ConstVerticesOnElementRangeType>::type ConstVerticesOnElementIteratorType;

      std::vector<DstElementType> dst_vertices;
      ConstVerticesOnElementRangeType vertices(src, 0);
      for (ConstVerticesOnElementIteratorType vit = vertices.begin(); vit != vertices.end(); ++vit)
        dst_vertices.push_back( (*this)(*vit) );

      DstElementType dst = viennagrid::make_element( dst_mesh(), src.tag(), dst_vertices.begin(), dst_vertices.end() );
      return dst;
    }

    DstMeshType const & dst_mesh() const { return dst_mesh_; }

  private:

    DstMeshType const & dst_mesh_;
    std::map<SrcElementType, DstElementType> vertex_map;

    boost::shared_ptr<nd_vertex_tree_node> ndtree;

    viennagrid_numeric eps;
  };
}

#endif
