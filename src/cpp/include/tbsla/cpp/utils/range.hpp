#ifndef TBSLA_CPP_RANGE
#define TBSLA_CPP_RANGE

#include <cstddef>

namespace tbsla { namespace utils { namespace range {

  size_t lnv(size_t size, long long int local_pos, long long int number_pos);
  size_t pflv(size_t size, long long int local_pos, long long int number_pos);
}}}

#endif
