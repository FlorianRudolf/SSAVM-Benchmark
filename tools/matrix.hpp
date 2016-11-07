#ifndef DISS_MATRIX_HPP_
#define DISS_MATRIX_HPP_

#include <stdlib.h>
#include <time.h>
#include "common.hpp"
#include "templated_mesh.hpp"


std::vector<double> randomVector(int size);

struct CSR
{
  std::vector<double> values;
  std::vector<int> indices;
  std::vector<int> offsets;

  int size() const { return offsets.size()-1; }

  int memory_usage() const
  {
    return sizeof(double)*values.size() +
           sizeof(int)*indices.size() +
           sizeof(int)*offsets.size();
  }
};



void prod(CSR const & matrix, std::vector<double> const & x, std::vector<double> & result);
CSR fromMesh(MeshType const & mesh);




struct TCSRInstance
{
  TCSRInstance() {}
  TCSRInstance(int local_size, int global_size)
  {
    local_to_global_mapping.resize(local_size);
    for (std::size_t i = 0; i != local_to_global_mapping.size(); ++i)
      local_to_global_mapping[i] = rand() % global_size;
  }

  std::vector<int> local_to_global_mapping;
  int memory_usage() const { return sizeof(int)*local_to_global_mapping.size(); }
};

struct TSCRTemplate
{
  CSR matrix;
  std::vector<TCSRInstance> instances;

  int memory_usage() const
  {
    int mem = matrix.memory_usage();
    for (std::size_t i = 0; i != instances.size(); ++i)
      mem += instances[i].memory_usage();
    return mem;
  }
};

struct TCSR
{
  std::vector<TSCRTemplate> templates;
  int size;

  int memory_usage() const
  {
    int mem = 0;
    for (std::size_t i = 0; i != templates.size(); ++i)
      mem += templates[i].memory_usage();
    return mem;
  }
};


void prod(TCSR const & matrix, std::vector<double> const & x, std::vector<double> & result);
TCSR fromMesh(boost::shared_ptr<TemplatedMesh> const & tm);



#endif
