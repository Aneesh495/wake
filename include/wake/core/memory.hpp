#pragma once

#include "wake/core/constants.hpp"
#include <cstdlib>
#include <new>
#include <memory>
#include <cstring>
#include <utility>

namespace wake {

// Low level aligned allocation and deallocation
void* allocate_aligned(std::size_t bytes, std::size_t alignment = SIMD_ALIGNMENT);
void free_aligned(void* ptr);

// Standard-compliant aligned STL allocator
template <typename T, std::size_t Alignment = SIMD_ALIGNMENT>
class AlignedAllocator {
public:
    using value_type = T;
    using pointer = T*;
    using const_pointer = const T*;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    template <typename U>
    struct rebind {
        using other = AlignedAllocator<U, Alignment>;
    };

    constexpr AlignedAllocator() noexcept = default;
    constexpr AlignedAllocator(const AlignedAllocator&) noexcept = default;
    template <typename U>
    constexpr AlignedAllocator(const AlignedAllocator<U, Alignment>&) noexcept {}

    pointer allocate(size_type n) {
        if (n == 0) return nullptr;
        void* ptr = allocate_aligned(n * sizeof(T), Alignment);
        if (!ptr) throw std::bad_alloc();
        return static_cast<pointer>(ptr);
    }

    void deallocate(pointer p, size_type) noexcept {
        free_aligned(p);
    }

    bool operator==(const AlignedAllocator&) const noexcept { return true; }
    bool operator!=(const AlignedAllocator&) const noexcept { return false; }
};

// RAII aligned buffer for SIMD vector distribution arrays
template <typename T>
class AlignedBuffer {
public:
    AlignedBuffer() : data_(nullptr), size_(0), capacity_(0) {}

    explicit AlignedBuffer(std::size_t size, std::size_t alignment = SIMD_ALIGNMENT)
        : data_(nullptr), size_(size), capacity_(size), alignment_(alignment) {
        if (size_ > 0) {
            data_ = static_cast<T*>(allocate_aligned(size_ * sizeof(T), alignment_));
            if (!data_) throw std::bad_alloc();
            std::memset(data_, 0, size_ * sizeof(T));
        }
    }

    ~AlignedBuffer() {
        if (data_) {
            free_aligned(data_);
            data_ = nullptr;
        }
    }

    AlignedBuffer(const AlignedBuffer& other)
        : size_(other.size_), capacity_(other.size_), alignment_(other.alignment_) {
        if (size_ > 0) {
            data_ = static_cast<T*>(allocate_aligned(size_ * sizeof(T), alignment_));
            if (!data_) throw std::bad_alloc();
            std::memcpy(data_, other.data_, size_ * sizeof(T));
        }
    }

    AlignedBuffer& operator=(const AlignedBuffer& other) {
        if (this != &other) {
            if (data_) free_aligned(data_);
            size_ = other.size_;
            capacity_ = other.size_;
            alignment_ = other.alignment_;
            if (size_ > 0) {
                data_ = static_cast<T*>(allocate_aligned(size_ * sizeof(T), alignment_));
                if (!data_) throw std::bad_alloc();
                std::memcpy(data_, other.data_, size_ * sizeof(T));
            } else {
                data_ = nullptr;
            }
        }
        return *this;
    }

    AlignedBuffer(AlignedBuffer&& other) noexcept
        : data_(other.data_), size_(other.size_), capacity_(other.capacity_), alignment_(other.alignment_) {
        other.data_ = nullptr;
        other.size_ = 0;
        other.capacity_ = 0;
    }

    AlignedBuffer& operator=(AlignedBuffer&& other) noexcept {
        if (this != &other) {
            if (data_) free_aligned(data_);
            data_ = other.data_;
            size_ = other.size_;
            capacity_ = other.capacity_;
            alignment_ = other.alignment_;
            other.data_ = nullptr;
            other.size_ = 0;
            other.capacity_ = 0;
        }
        return *this;
    }

    void resize(std::size_t new_size, std::size_t alignment = SIMD_ALIGNMENT) {
        if (new_size == size_ && alignment_ == alignment) return;
        T* new_data = nullptr;
        if (new_size > 0) {
            new_data = static_cast<T*>(allocate_aligned(new_size * sizeof(T), alignment));
            if (!new_data) throw std::bad_alloc();
            std::memset(new_data, 0, new_size * sizeof(T));
            if (data_ && size_ > 0) {
                std::size_t copy_bytes = (std::min(size_, new_size)) * sizeof(T);
                std::memcpy(new_data, data_, copy_bytes);
            }
        }
        if (data_) free_aligned(data_);
        data_ = new_data;
        size_ = new_size;
        capacity_ = new_size;
        alignment_ = alignment;
    }

    void zero() {
        if (data_ && size_ > 0) {
            std::memset(data_, 0, size_ * sizeof(T));
        }
    }

    T* data() noexcept { return data_; }
    const T* data() const noexcept { return data_; }

    T& operator[](std::size_t index) noexcept { return data_[index]; }
    const T& operator[](std::size_t index) const noexcept { return data_[index]; }

    std::size_t size() const noexcept { return size_; }
    std::size_t bytes() const noexcept { return size_ * sizeof(T); }

private:
    T* data_{nullptr};
    std::size_t size_{0};
    std::size_t capacity_{0};
    std::size_t alignment_{SIMD_ALIGNMENT};
};

} // namespace wake
