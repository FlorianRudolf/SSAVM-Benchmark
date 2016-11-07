#include <iostream>
#include <vector>
#include <algorithm>
#include <set>

typedef std::vector< std::pair<int,int> > GraphType;
typedef std::vector< std::pair<int,int> > NeighborsType;
typedef std::vector<int> CircleType;

NeighborsType graph_neighbors(GraphType const & graph, int node)
{
  NeighborsType result;

  for (std::size_t i = 0; i != graph.size(); ++i)
  {
    if (graph[i].first == node)
      result.push_back( std::make_pair(i, graph[i].second) );

    if (graph[i].second == node)
      result.push_back( std::make_pair(i, graph[i].first) );
  }

  return result;
}

int node_count(GraphType const & graph)
{
  int result = 0;
  for (std::size_t i = 0; i != graph.size(); ++i)
  {
    result = std::max(result, std::max(graph[i].first, graph[i].second));
  }
  return result+1;
}


void find_circle(GraphType const & graph,
                 std::vector<NeighborsType> const & neighbors,
                 int starting_node,
                 int current_node,
                 std::vector<bool> edge_flags,
                 CircleType circle,
                 std::set<CircleType> & circles)
{
  circle.push_back(current_node);

  for (std::size_t i = 0; i != neighbors[current_node].size(); ++i)
  {
    if (edge_flags[neighbors[current_node][i].first])
      continue;

    if (neighbors[current_node][i].second == starting_node)
    {
      CircleType::iterator it = std::min_element(circle.begin(), circle.end());
      std::rotate(circle.begin(), it, circle.end());
      if (circle[1] > circle.back())
      {
        std::reverse(circle.begin(), circle.end());
        std::rotate(circle.rbegin(), circle.rbegin()+1, circle.rend());
      }

      circles.insert(circle);
      return;
    }

    edge_flags[neighbors[current_node][i].first] = true;
    find_circle(graph, neighbors, starting_node, neighbors[current_node][i].second, edge_flags, circle, circles);
    edge_flags[neighbors[current_node][i].first] = false;
  }
}


std::set<CircleType> find_all_circles(GraphType const & graph)
{
  std::size_t nc = node_count(graph);

  std::set<CircleType> circles;
  std::vector<NeighborsType> neighbors(nc);
  for (int i = 0; i != nc; ++i)
    neighbors[i] = graph_neighbors(graph, i);

  for (int i = 0; i != nc; ++i)
  {
    std::vector<bool> edge_flags(graph.size(), false);
    std::vector<int> circle;

    find_circle(graph, neighbors, i, i, edge_flags, circle, circles);
  }

  return circles;
}

int main()
{
  GraphType graph;

  graph.push_back( std::make_pair(0,1) );
  graph.push_back( std::make_pair(0,2) );
  graph.push_back( std::make_pair(1,3) );
  graph.push_back( std::make_pair(1,4) );
  graph.push_back( std::make_pair(2,6) );
  graph.push_back( std::make_pair(3,6) );
  graph.push_back( std::make_pair(3,7) );
  graph.push_back( std::make_pair(4,7) );
  graph.push_back( std::make_pair(5,6) );
  graph.push_back( std::make_pair(6,7) );
  graph.push_back( std::make_pair(7,8) );





  std::set<CircleType> circles = find_all_circles(graph);

  for (std::set<CircleType>::iterator it = circles.begin(); it != circles.end(); ++it)
  {
    std::cout << "  found circle: ";
    for (std::size_t i = 0; i != (*it).size(); ++i)
      std::cout << (*it)[i] << ",";
    std::cout << std::endl;
  }



  return 0;
}