#ifndef cub_module_radix_sort_h
#define cub_module_radix_sort_h

#include <cstdint>

extern "C"
{
/** sorts keys and generates a permutation array.
 * @param[in] keys_in keys to sort
 * @param[out] keys_out array for the sorted keys
 * @param[in] idx_in an empty array of the same length as keys_in
 * @param[out] idx_out array for the sorted indices. this is the result
 * @param[in] head_idx the 1 based index into the keys
 * @param[in] num_elem the number of keys to sort
 * @param[in] start_bit the first bit to use in the radix sort
 * @param[in] end_bit the last bit to use in the radix sort
 */
void cub_module_radix_sort(int64_t *keys_in, int64_t *keys_out,
                           int32_t *idx_in, int32_t *idx_out,
                           int32_t head_idx,  int32_t num_elem,
                           int32_t start_bit, int32_t end_bit);

/** sorts pairs in place, ping-ponging between two key buffers and two value
 * buffers. Unlike cub_module_radix_sort, the work space holds no copy of the
 * keys and values, but both buffers of each pair are overwritten.
 * @param[in,out] keys keys to sort
 * @param[in,out] keys_alt second key buffer of the same length
 * @param[in,out] vals values carried with the keys
 * @param[in,out] vals_alt second value buffer of the same length
 * @param[in] num_elem the number of keys to sort
 * @param[in] start_bit the first bit to use in the radix sort
 * @param[in] end_bit the last bit to use in the radix sort (exclusive)
 * @returns 0 if the sorted values are in vals, 1 if they are in vals_alt
 */
int32_t cub_module_radix_sort_inplace(int64_t *keys, int64_t *keys_alt,
                                      int32_t *vals, int32_t *vals_alt,
                                      int32_t num_elem,
                                      int32_t start_bit, int32_t end_bit);

/** same as cub_module_radix_sort_inplace for 32 bit keys. */
int32_t cub_module_radix_sort_inplace32(uint32_t *keys, uint32_t *keys_alt,
                                        int32_t *vals, int32_t *vals_alt,
                                        int32_t num_elem,
                                        int32_t start_bit, int32_t end_bit);

/** allocate the workspace for a sort of num_elem.
 * note: the sort call will always check if there is enough space and allocate
 * as needed.
 */
void cub_module_radix_sort_allocate_workspace(int32_t num_elem);
}

#endif
