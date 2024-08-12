//
// Created by cepheid on 5/16/24.
//

#ifndef MEMORY_H
#define MEMORY_H

#include <cstdint>
#include <cassert>

#define DEFAULT_ALIGN (2 * sizeof(void*))

namespace memory {
  inline bool is_power_of_two(const uintptr_t x) {
    return (x & (x-1)) == 0;
  }

  inline uintptr_t align_forward(const uintptr_t ptr, const size_t align)
  {
    assert(is_power_of_two(align));

    auto p = ptr;
    const auto a = static_cast<uintptr_t>(align);
    // Same as (p % a) but faster as 'a' is a power of two
    const auto modulo = p & (a - 1);

    if ((p & (a - 1)) != 0) {
      // If 'p' address is not aligned, push the address to the
      // next value which is aligned
      p += a - modulo;
    }
    return p;
  }

  // struct Arena {
  //   unsigned char* buf;
  //   size_t         buf_len;
  //   size_t         prev_offset; // This will be useful for later on
  //   size_t         curr_offset;
  //
  //   void* alloc_align(size_t size, size_t align);
  //   void* alloc(size_t size);
  //   void* resize_align(void* old_memory, size_t old_size, size_t new_size, size_t align=DEFAULT_ALIGN);
  //   void* resize(void* old_memory, size_t old_size, size_t new_size);
  //   void free_all();
  //   void init(void* backing_buffer, size_t backing_buffer_length);
  //
  //   static void free([[maybe_unused]] void* ptr) {} // Do nothing
  // };

  struct Pool_Free_Node {
    Pool_Free_Node *next;
  };

  struct Pool {
    unsigned char *buf;
    size_t buf_len;
    size_t chunk_size;

    Pool_Free_Node *head; // Free List Head

    void* alloc();

    void init(void* buffer, size_t buffer_length, size_t block_size, size_t block_align=DEFAULT_ALIGN);
    void free(void *ptr);
    void free_all();
  };


}
#endif //MEMORY_H
