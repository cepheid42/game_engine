#ifndef MEMORY_H
#define MEMORY_H

#include <cstdint>
#include <cstring>
#include <cassert>

#define DEFAULT_ALIGN (2 * sizeof(void*))

namespace memory {
  inline bool is_power_of_two(const uintptr_t x) {
    return (x & (x-1)) == 0;
  }

  inline uintptr_t align_forward(const uintptr_t ptr, const std::size_t align)
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

  struct Pool_Free_Node {
    Pool_Free_Node *next;
  };

  template<typename T>
  struct Pool {
    using pointer = T*;

    Pool();
    Pool(const std::size_t size);

    T* buf;
    std::size_t buf_len;
    std::size_t chunk_size;

    Pool_Free_Node* head; // Free List Head

    pointer allocate(std::size_t n) {
      auto* node = head;
      if (!node) {
        assert(false && "Pool allocator has no free memory");
        return nullptr;
      }

      head = head->next;
      return std::memset(node, 0, chunk_size);
    }

    void deallocate(pointer p, std::size_t n) {

    }


    void free(void *ptr);
    void free_all();
  };


}
#endif //MEMORY_H
