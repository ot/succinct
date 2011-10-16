#pragma once

#include <stdint.h>
#include "intrinsics.hpp"

namespace succinct {
namespace broadword {

    static const uint64_t ones_step_4  = 0x1111111111111111ULL;
    static const uint64_t ones_step_8  = 0x0101010101010101ULL;
    static const uint64_t ones_step_9  = 1ULL << 0 | 1ULL << 9 | 1ULL << 18 | 1ULL << 27 | 1ULL << 36 | 1ULL << 45 | 1ULL << 54;
    static const uint64_t msbs_step_8  = 0x80ULL * ones_step_8;
    static const uint64_t msbs_step_9  = 0x100ULL * ones_step_9;
    static const uint64_t incr_step_8  = 0x80ULL << 56 | 0x40ULL << 48 | 0x20ULL << 40 | 0x10ULL << 32 | 0x8ULL << 24 | 0x4ULL << 16 | 0x2ULL << 8 | 0x1;
    static const uint64_t inv_count_step_9 = 1ULL << 54 | 2ULL << 45 | 3ULL << 36 | 4ULL << 27 | 5ULL << 18 | 6ULL << 9 | 7ULL;

    static const uint64_t magic_mask_1 = 0x5555555555555555ULL;
    static const uint64_t magic_mask_2 = 0x3333333333333333ULL;
    static const uint64_t magic_mask_3 = 0x0F0F0F0F0F0F0F0FULL;
    static const uint64_t magic_mask_4 = 0x00FF00FF00FF00FFULL;
    static const uint64_t magic_mask_5 = 0x0000FFFF0000FFFFULL;
    static const uint64_t magic_mask_6 = 0x00000000FFFFFFFFULL;

    inline uint64_t leq_step_8(uint64_t x, uint64_t y)
    {
        return ((((y | msbs_step_8) - (x & ~msbs_step_8)) ^ (x ^ y)) & msbs_step_8) >> 7;
    }

    inline uint64_t uleq_step_8(uint64_t x, uint64_t y)
    {
        return (((((y | msbs_step_8) - (x & ~msbs_step_8)) ^ (x ^ y)) ^ (x & ~y)) & msbs_step_8) >> 7;
    }

    inline uint64_t zcompare_step_8(uint64_t x)
    {
        return ((x | ((x | msbs_step_8) - ones_step_8)) & msbs_step_8) >> 7;
    }
    
    inline uint64_t uleq_step_9(uint64_t x, uint64_t y) 
    {
        return (((((y | msbs_step_9) - (x & ~msbs_step_9)) | (x ^ y)) ^ (x & ~y)) & msbs_step_9 ) >> 8;
    }

    inline uint64_t byte_counts(uint64_t x) 
    {
        x = x - ((x & 0xa * ones_step_4) >> 1);
        x = (x & 3 * ones_step_4) + ((x >> 2) & 3 * ones_step_4);
        x = (x + (x >> 4)) & 0x0f * ones_step_8;
        return x;
    }

    inline uint64_t bytes_sum(uint64_t x)
    {
        return x * ones_step_8 >> 56;
    }

    inline uint64_t popcount(uint64_t x) 
    {
	// TODO(ot): may be worth to try the popcnt instruction
	// on platforms that support it
        return bytes_sum(byte_counts(x));
    }
        
    inline uint64_t reverse_bytes(uint64_t x) 
    {
#if SUCCINCT_USE_INTRINSICS
        return intrinsics::byteswap64(x);
#else
        x = ((x >> 8)  & magic_mask_4) | ((x & magic_mask_4) << 8);
        x = ((x >> 16) & magic_mask_5) | ((x & magic_mask_5) << 16);
        x = ((x >> 32)               ) | ((x               ) << 32);
        return x;
#endif
    }

    inline uint64_t select_in_word(const uint64_t x, const uint64_t k) 
    {
        assert(k < popcount(x));

        uint64_t byte_sums = byte_counts(x) * ones_step_8;

        const uint64_t k_step_8 = k * ones_step_8;
        const uint64_t place = (leq_step_8(byte_sums, k_step_8) * ones_step_8 >> 53) & ~0x7;

        const uint64_t byte_rank = k - (((byte_sums << 8) >> place) & 0xFF);
        const uint64_t spread_bits = (x >> place & 0xFF) * ones_step_8 & incr_step_8;
        const uint64_t bit_sums = zcompare_step_8(spread_bits) * ones_step_8;

        const uint64_t byte_rank_step_8 = byte_rank * ones_step_8;

        return place + (leq_step_8(bit_sums, byte_rank_step_8) * ones_step_8 >> 56);
    }
    
    inline uint64_t same_msb(uint64_t x, uint64_t y)
    {
        return (x ^ y) <= (x & y);
    }

    inline uint8_t msb(uint64_t x, unsigned long& ret) {
#if SUCCINCT_USE_INTRINSICS
        return intrinsics::bsr64(&ret, x);
#else
        if (!x) {
            return false;
        }
        // From Knuth
        ret = (unsigned long)
            ((same_msb(x, x & ~magic_mask_1)     ) +
             (same_msb(x, x & ~magic_mask_2) << 1) +
             (same_msb(x, x & ~magic_mask_3) << 2) +
             (same_msb(x, x & ~magic_mask_4) << 3) +
             (same_msb(x, x & ~magic_mask_5) << 4) +
             (same_msb(x, x & ~magic_mask_6) << 5) );
        return true;
#endif
    }

    inline uint8_t msb(uint64_t x) {
        assert(x);
        unsigned long ret;
        msb(x, ret);
        return (uint8_t)ret;
    }

    namespace detail {
        const uint8_t lsb_bit_table[64] = {
            63,  0, 58,  1, 59, 47, 53,  2,
            60, 39, 48, 27, 54, 33, 42,  3,
            61, 51, 37, 40, 49, 18, 28, 20,
            55, 30, 34, 11, 43, 14, 22,  4,
            62, 57, 46, 52, 38, 26, 32, 41,
            50, 36, 17, 19, 29, 10, 13, 21,
            56, 45, 25, 31, 35, 16,  9, 12,
            44, 24, 15,  8, 23,  7,  6,  5
        };
    }
 
    inline uint8_t lsb(uint64_t x, unsigned long& ret) {
#if SUCCINCT_USE_INTRINSICS
        return intrinsics::bsf64(&ret, x);
#else
        if (!x) {
            return false;
        }
        // Adapted from LSB of Chess Programming Wiki
        const uint64_t debruijn64 = 0x07EDD5E59A4E28C2ULL;
        ret = detail::lsb_bit_table[((x & -int64_t(x)) * debruijn64) >> 58];
        return true;
#endif
    }

    inline uint8_t lsb(uint64_t x) {
        assert(x);
        unsigned long ret;
        lsb(x, ret);
        return (uint8_t)ret;
    }
}
}
