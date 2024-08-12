//
// Created by cepheid on 5/17/24.
//

#include <cstdint>
#include <cstring>

#include "memory.h"

namespace memory {
  // void* Arena::alloc_align(size_t size, size_t align=DEFAULT_ALIGN)
  // {
  //   // Align 'curr_offset' forward to the specified alignment
  //   auto curr_ptr = reinterpret_cast<uintptr_t>(buf + curr_offset);
  //   auto offset = align_forward(curr_ptr, align);
  //   offset -= reinterpret_cast<uintptr_t>(buf); // Change to relative offset
  //
  //   // Check to see if the backing memory has space left
  //   if (offset + size <= buf_len) {
  //     void *ptr = &buf[offset];
  //     prev_offset = offset;
  //     curr_offset = offset + size;
  //
  //     // Zero new memory by default
  //     std::memset(ptr, 0, size);
  //     return ptr;
  //   }
  //   // Return nullptr if the arena is out of memory (or handle differently)
  //   return nullptr;
  // }
  //
  // void* Arena::alloc(size_t size) {
  //   return alloc_align(size);
  // }
  //
  // void Arena::init(void* backing_buffer, size_t backing_buffer_length)
  // {
  //   buf = static_cast<unsigned char*>(backing_buffer);
  //   buf_len = backing_buffer_length;
  //   curr_offset = 0;
  //   prev_offset = 0;
  // }
  //
  // void* Arena::resize_align(void* old_memory, size_t old_size, size_t new_size, size_t align)
  // {
  //   assert(is_power_of_two(align));
  //   auto* old_mem = static_cast<unsigned char*>(old_memory);
  //
  //   if (old_mem == nullptr or old_size == 0) {
  //     return alloc_align(new_size, align);
  //   }
  //
  //   if (buf <= old_mem and old_mem < buf + buf_len) {
  //     if (buf + prev_offset == old_mem) {
  //       curr_offset = prev_offset + new_size;
  //       if (new_size > old_size) {
  //         // Zero the new memory by default
  //         std::memset(&buf[curr_offset], 0, new_size - old_size);
  //       }
  //       return old_memory;
  //     }
  //     void* new_memory = alloc_align(new_size, align);
  //     size_t copy_size = old_size < new_size ? old_size : new_size;
  //     // Copy across old memory to the new memory
  //     std::memmove(new_memory, old_memory, copy_size);
  //     return new_memory;
  //   }
  //
  //   assert(false && "Memory is out of bounds of the buffer in this arena");
  //   return nullptr;
  // }
  //
  // // Because C doesn't have default parameters
  // void* Arena::resize(void* old_memory, size_t old_size, size_t new_size) {
  //   return resize_align(old_memory, old_size, new_size);
  // }
  //
  // void Arena::free_all() {
  //   curr_offset = 0;
  //   prev_offset = 0;
  // }

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
