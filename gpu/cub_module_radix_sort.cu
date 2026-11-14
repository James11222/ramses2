#include "cub_module_radix_sort.h"
#include <cub/cub.cuh>
#include <cstdint>

/// working space for device scan
void *g_cub_module_radix_sort_workspace = nullptr;
size_t g_cub_module_radix_sort_work_size = 0;
int g_cub_module_radix_sort_verbose = 0;


#define CUDA_CHECK(call) \
{ \
  cudaError_t ierr = call; \
  if (cudaSuccess != ierr) \
  { \
    std::cerr << "ERROR: " << __FILE__ << ":" << __LINE__ << std::endl \
              << "ERROR: " << cudaGetErrorString(ierr) << std::endl; \
    abort(); \
  } \
}

///  checks the alighnment of the passed ptr and aborts if it is not correctly aligned
#if !defined(NDEBUG)
// @returns true if the pointer is aligned to it's type's size
template<typename T>
bool is_aligned(T *ptr)
{
  if (0 == uintptr_t(ptr) % sizeof(T))
    return true;
  return false;
}

#define CHECK_ALIGN(var) \
if (!is_aligned(var)) \
{ \
  std::cerr << "ERROR: " #var " is not " << sizeof(decltype(var)) \
            << " byte aligned." << std::endl; \
  abort(); \
}

#else
#define CHECK_ALIGN(var)
#endif

// --------------------------------------------------------------------------
void cub_module_radix_sort_allocate_workspace(int32_t num_elem)
{
  // check if we have enough work space
  size_t work_size = 0;
  CUDA_CHECK(cub::DeviceRadixSort::SortPairs(nullptr, work_size,
                                             (int64_t*)nullptr, (int64_t*)nullptr,
                                             (int32_t*)nullptr, (int32_t*)nullptr,
                                             num_elem))
  if (g_cub_module_radix_sort_work_size < work_size)
  {
    // free old work space
    if (g_cub_module_radix_sort_workspace)
    {
      cudaFree(g_cub_module_radix_sort_workspace);
      g_cub_module_radix_sort_workspace = nullptr;
      g_cub_module_radix_sort_work_size = 0;
    }

    // allocate new work space
    g_cub_module_radix_sort_work_size = work_size;
    CUDA_CHECK(cudaMalloc(&g_cub_module_radix_sort_workspace,
                          g_cub_module_radix_sort_work_size))
    if (0 < g_cub_module_radix_sort_verbose)
      std::cerr << "allocated " << work_size << " bytes on the device" << std::endl;
  }
}

/// grow the retained work space to at least work_size bytes
static void cub_module_radix_sort_reserve(size_t work_size)
{
  if (g_cub_module_radix_sort_work_size < work_size)
  {
    if (g_cub_module_radix_sort_workspace)
    {
      cudaFree(g_cub_module_radix_sort_workspace);
      g_cub_module_radix_sort_workspace = nullptr;
      g_cub_module_radix_sort_work_size = 0;
    }
    CUDA_CHECK(cudaMalloc(&g_cub_module_radix_sort_workspace, work_size))
    g_cub_module_radix_sort_work_size = work_size;
    if (0 < g_cub_module_radix_sort_verbose)
      std::cerr << "allocated " << work_size << " bytes on the device" << std::endl;
  }
}

/// sort pairs in place over two key buffers and two value buffers.
/// CUB ping-pongs between the buffers of each pair, so it needs no third
/// buffer of num_elem keys and values in the work space (only histograms).
/// Both buffers of each pair are overwritten.
/// @returns 0 if the sorted values are in vals, 1 if they are in vals_alt
template<typename KeyT>
static int32_t radix_sort_inplace(KeyT *keys, KeyT *keys_alt,
                                  int32_t *vals, int32_t *vals_alt,
                                  int32_t num_elem,
                                  int32_t start_bit, int32_t end_bit)
{
  CHECK_ALIGN(keys)
  CHECK_ALIGN(keys_alt)
  CHECK_ALIGN(vals)
  CHECK_ALIGN(vals_alt)

  cub::DoubleBuffer<KeyT> d_keys(keys, keys_alt);
  cub::DoubleBuffer<int32_t> d_vals(vals, vals_alt);

  size_t work_size = 0;
  CUDA_CHECK(cub::DeviceRadixSort::SortPairs(nullptr, work_size, d_keys, d_vals,
                                             num_elem, start_bit, end_bit))
  cub_module_radix_sort_reserve(work_size);

  CUDA_CHECK(cub::DeviceRadixSort::SortPairs(g_cub_module_radix_sort_workspace,
                                             work_size, d_keys, d_vals,
                                             num_elem, start_bit, end_bit))

  return d_vals.Current() == vals ? 0 : 1;
}

int32_t cub_module_radix_sort_inplace(int64_t *keys, int64_t *keys_alt,
                                      int32_t *vals, int32_t *vals_alt,
                                      int32_t num_elem,
                                      int32_t start_bit, int32_t end_bit)
{
  return radix_sort_inplace(keys, keys_alt, vals, vals_alt, num_elem, start_bit, end_bit);
}

int32_t cub_module_radix_sort_inplace32(uint32_t *keys, uint32_t *keys_alt,
                                        int32_t *vals, int32_t *vals_alt,
                                        int32_t num_elem,
                                        int32_t start_bit, int32_t end_bit)
{
  return radix_sort_inplace(keys, keys_alt, vals, vals_alt, num_elem, start_bit, end_bit);
}

/// input an array of keys to sort
void cub_module_radix_sort(int64_t *keys_in, int64_t *keys_out,
                           int32_t *idx_in, int32_t *idx_out,
                           int32_t idx0, int32_t num_elem,
                           int32_t start_bit, int32_t end_bit)
{
  if (1 < g_cub_module_radix_sort_verbose)
    std::cerr << "cub_module_radix_sort" << std::endl;

  // check for misaligned buffers.
  CHECK_ALIGN(keys_in)
  CHECK_ALIGN(keys_out)
  CHECK_ALIGN(idx_in)
  CHECK_ALIGN(idx_out)

  // allocate workspace
  cub_module_radix_sort_allocate_workspace(num_elem);

  CUDA_CHECK(cub::DeviceRadixSort::SortPairs(g_cub_module_radix_sort_workspace,
                                             g_cub_module_radix_sort_work_size,
                                             keys_in, keys_out, idx_in, idx_out,
                                             num_elem, start_bit, end_bit))
}
