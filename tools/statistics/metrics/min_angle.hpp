#ifndef VIENNAMESH_STATISTICS_METRICS_MIN_ANGLE_HPP
#define VIENNAMESH_STATISTICS_METRICS_MIN_ANGLE_HPP

/* ============================================================================
   Copyright (c) 2011-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
                ViennaMesh - The Vienna Meshing Framework
                            -----------------

                    http://viennamesh.sourceforge.net/

   License:         MIT (X11), see file LICENSE in the base directory
=============================================================================== */

#include "forwards.hpp"
#include "viennagrid/algorithm/angle.hpp"

namespace viennamesh
{
  namespace metrics
  {

    template<typename ElementT, typename NumericLimitsT>
    viennagrid_numeric min_angle_impl(ElementT const & element, NumericLimitsT, viennagrid::triangle_tag)
    {
      typedef typename viennagrid::result_of::point<ElementT>::type PointType;

//       typedef typename PointAccessorT::value_type            PointType;
//       typedef typename viennagrid::result_of::coord<PointType>::type            NumericType;

      PointType a = viennagrid::get_point( viennagrid::vertices(element)[0] );
      PointType b = viennagrid::get_point( viennagrid::vertices(element)[1] );
      PointType c = viennagrid::get_point( viennagrid::vertices(element)[2] );

      viennagrid_numeric alpha = viennagrid::angle( b, c, a );
      viennagrid_numeric beta = viennagrid::angle( a, c, b );
      viennagrid_numeric gamma = M_PI - alpha - beta;

      return std::min(std::min( alpha, beta ), gamma);
    }


    template<typename ElementT, typename NumericLimitsT>
    viennagrid_numeric min_angle_impl(ElementT const & element, NumericLimitsT, viennagrid::tetrahedron_tag)
    {
      typedef typename viennagrid::result_of::point<ElementT>::type PointType;

//       typedef typename PointAccessorT::value_type            PointType;
//       typedef typename viennagrid::result_of::coord<PointType>::type            NumericType;

      PointType a = viennagrid::get_point( viennagrid::vertices(element)[0] );
      PointType b = viennagrid::get_point( viennagrid::vertices(element)[1] );
      PointType c = viennagrid::get_point( viennagrid::vertices(element)[2] );
      PointType d = viennagrid::get_point( viennagrid::vertices(element)[3] );

      viennagrid_numeric alpha = viennagrid::solid_angle( b, c, d, a );
      viennagrid_numeric beta = viennagrid::solid_angle( a, c, d, b );
      viennagrid_numeric gamma = viennagrid::solid_angle( a, b, d, c );
      viennagrid_numeric delta = viennagrid::solid_angle( a, b, c, d );

      return std::min( std::min( alpha, beta ), std::min(gamma, delta) );
    }
  }


  struct min_angle_tag {};

  template<typename ElementT, typename NumericLimitsT>
  viennagrid_numeric min_angle( ElementT const & element, NumericLimitsT numeric_limits )
  {
    if (element.tag().is_triangle())
    {
      return metrics::min_angle_impl(element, numeric_limits, viennagrid::triangle_tag());
    }
    else if (element.tag().is_triangle())
    {
      return metrics::min_angle_impl(element, numeric_limits, viennagrid::tetrahedron_tag());
    }

    return -1;
  }

  template<typename ElementT>
  viennagrid_numeric min_angle(ElementT const & element )
  {
    return min_angle(element, std::numeric_limits<viennagrid_numeric>() );
  }

//   template<typename ElementT>
//   typename viennagrid::result_of::coord< ElementT >::type min_angle(ElementT const & element)
//   {
//     return min_angle( viennagrid::root_mesh_point_accessor(element), element);
//   }


  namespace detail
  {
    template<typename PointAccessorT, typename ElementT, typename NumericLimitsT>
    typename viennagrid::result_of::coord<typename PointAccessorT::value_type>::type metric( PointAccessorT const point_accessor, ElementT const & element, NumericLimitsT numeric_limits, min_angle_tag)
    {
      return min_angle(point_accessor, element, numeric_limits);
    }

    template<typename PointAccessorT, typename ElementT>
    typename viennagrid::result_of::coord<typename PointAccessorT::value_type>::type metric( PointAccessorT const point_accessor, ElementT const & element, min_angle_tag)
    {
      return min_angle(point_accessor, element);
    }

    template<typename ElementT>
    typename viennagrid::result_of::coord< ElementT >::type metric( ElementT const & element, min_angle_tag)
    {
      return min_angle(element);
    }
  }

}


#endif
