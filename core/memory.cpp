//
// Created by cepheid on 5/17/24.
//

#include <cstdint>
#include <cstring>

#include "memory.h"

namespace memory {
  void Pool::init(void* buffer_, size_t buffer_length_, size_t block_size_, size_t block_align_)
  {
    // Align backing buffer to the specified chunk alignment
    auto initial_start = reinterpret_cast<uintptr_t>(buffer_);
    auto start = align_forward(initial_start, block_align_);
    buffer_length_ -= start - initial_start;

    // Align chunk size up to the required chunk_alignment
    chunk_size = align_forward(chunk_size, block_align_);

    // Assert that the parameters passed are valid
    assert(chunk_size >= sizeof(Pool_Free_Node) and "Chunk size is too small");
    assert(buffer_length_ >= chunk_size and "Backing buffer length is smaller than the chunk size");

    // Store the adjusted parameters
    buf = static_cast<unsigned char*>(buffer_);
    buf_len = buffer_length_;
    chunk_size = block_size_;
    head = nullptr; // Free List Head

    // Set up the free list for free chunks
    free_all();
  }

  void* Pool::alloc() {
    // Get latest free node
    auto* node = head;

    if (node == nullptr) {
      assert(false && "Pool allocator has no free memory");
      return nullptr;
    }

    // Pop free node
    head = head->next;

    // Zero memory by default
    return std::memset(node, 0, chunk_size);
  }

  void Pool::free(void *ptr) {
    void* start = buf;
    void* end = &buf[buf_len];

    if (ptr == nullptr) {
      // Ignore nullptr pointers
      return;
    }

    if (!(start <= ptr and ptr < end)) {
      assert(false && "Memory is out of bounds of the buffer in this pool");
      return;
    }

    // Push free node
    auto* node = static_cast<Pool_Free_Node*>(ptr);
    node->next = head;
    head = node;
  }

  void Pool::free_all() {
    size_t chunk_count = buf_len / chunk_size;

    // Set all chunks to be free
    for (size_t i = 0; i < chunk_count; i++) {
      void* ptr = &buf[i * chunk_size];
      auto* node = static_cast<Pool_Free_Node*>(ptr);
      // Push free node onto thte free list
      node->next = head;
      head = node;
    }
  }
}
