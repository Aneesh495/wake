#include "wake/core/memory.hpp"
#include <cstdlib>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

namespace wake {

void* allocate_aligned(std::size_t bytes, std::size_t alignment) {
    if (bytes == 0) return nullptr;

#if defined(_MSC_VER)
    return _aligned_malloc(bytes, alignment);
#else
    void* ptr = nullptr;
    // posix_memalign requires alignment to be a power of two and a multiple of sizeof(void*)
    std::size_t effective_align = alignment;
    if (effective_align < sizeof(void*)) {
        effective_align = sizeof(void*);
    }
    if (posix_memalign(&ptr, effective_align, bytes) != 0) {
        return nullptr;
    }
    return ptr;
#endif
}

void free_aligned(void* ptr) {
    if (!ptr) return;

#if defined(_MSC_VER)
    _aligned_free(ptr);
#else
    std::free(ptr);
#endif
}

} // namespace wake
