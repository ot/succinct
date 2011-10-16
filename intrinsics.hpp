#pragma once

#include <stdint.h>

#ifndef SUCCINCT_USE_INTRINSICS
#define SUCCINCT_USE_INTRINSICS 1
#endif

#if SUCCINCT_USE_INTRINSICS

#include <xmmintrin.h>

#if defined(__GNUC__) || defined(__clang__)
#    define __INTRIN_INLINE inline __attribute__((__always_inline__))
#elif defined(_MSC_VER)
#    define __INTRIN_INLINE inline __forceinline
#else
#    error Unsupported platform
#endif

namespace succinct {
namespace intrinsics {

    __INTRIN_INLINE uint64_t byteswap64(uint64_t value)
    {
#if defined(__GNUC__) || defined(__clang__)
        return __builtin_bswap64(value);
#elif defined(_MSC_VER)
        return _byteswap_uint64(value);
#else
#     error Unsupported platform
#endif
    }

    __INTRIN_INLINE bool bsf64(unsigned long* const index, const uint64_t mask)
    {
#if defined(__GNUC__) || defined(__clang__)
	if (mask) {
	    __asm__("bsf %[mask], %[index]" : [index] "=r" (*index) : [mask] "mr" (mask));
	    return true;
	} else {
	    return false;
	}
#elif defined(_MSC_VER)
        return _BitScanForward64(index, mask) != 0;
#else
#     error Unsupported platform
#endif
    }

    __INTRIN_INLINE bool bsr64(unsigned long* const index, const uint64_t mask)
    {
#if defined(__GNUC__) || defined(__clang__)
	if (mask) {
	    __asm__("bsr %[mask], %[index]" : [index] "=r" (*index) : [mask] "mr" (mask));
	    return true;
	} else {
	    return false;
	}
#elif defined(_MSC_VER)
        return _BitScanReverse64(index, mask) != 0;
#else
#     error Unsupported platform
#endif
    }

}
}

#endif


namespace succinct {
namespace intrinsics {

    template <typename T>
    __INTRIN_INLINE void prefetch(T const* ptr)
    {
#if SUCCINCT_USE_INTRINSICS
        _mm_prefetch((const char*)ptr, _MM_HINT_T0);
#endif
    }
}
}

