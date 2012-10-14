#include "bp_vector.hpp"
#include "util.hpp"

namespace succinct {
    
    namespace {
        class excess_tables
        {
        public:
            excess_tables() {
                for (int c = 0; c < 256; ++c) {
                    for (unsigned char i = 0; i < 9; ++i) {
                        m_fwd_pos[c][i] = 0;
                        m_bwd_pos[c][i] = 0;
                    }
                    // populate m_fwd_pos
                    int excess = 0;
		    int fwd_min = 0;
                    for (char i = 0; i < 8; ++i) {
                        if ((c >> i) & 1) { // opening
                            ++excess;
                        } else { // closing
                            --excess;
                            if (excess < 0 && 
                                m_fwd_pos[c][-excess] == 0) { // not already found
                                m_fwd_pos[c][-excess] = i + 1;
                            }
                        }
			fwd_min = std::max(-excess, fwd_min);
                    }
		    m_fwd_min[c] = fwd_min;

                    // populate m_bwd_pos 
                    excess = 0;
                    int bwd_min = 0;
                    for (unsigned char i = 0; i < 8; ++i) {
                        if ((c << i) & 128) { // opening
                            ++excess;
                            if (excess > 0 && 
                                m_bwd_pos[c][(unsigned char)excess] == 0) { // not already found
                                m_bwd_pos[c][(unsigned char)excess] = i + 1;
                            }
                        } else { // closing
                            --excess;
                        }
			bwd_min = std::max(excess, bwd_min);
                    }
		    m_bwd_min[c] = bwd_min;
                }
            }

            uint8_t m_fwd_pos[256][9];
            uint8_t m_bwd_pos[256][9];
            
	    uint8_t m_bwd_min[256];
	    uint8_t m_fwd_min[256];
        };
        
        const static excess_tables tables;

        inline bool find_close_in_word(uint64_t word, uint64_t byte_counts, bp_vector::excess_t cur_exc, uint64_t& ret) 
        {
            assert(cur_exc > 0 && cur_exc <= 64);
            const uint64_t cum_exc_step_8 = (cur_exc + ((2 * byte_counts - 8 * broadword::ones_step_8) << 8)) * broadword::ones_step_8;
	    
	    uint64_t min_exc_step_8 = 0;
	    for (int i = 0; i < 8; ++i) {
		size_t shift = i * 8;
		min_exc_step_8 |= ((uint64_t)(tables.m_fwd_min[(word >> shift) & 0xFF])) << shift;
	    }

	    const uint64_t has_result = broadword::leq_step_8(cum_exc_step_8, min_exc_step_8);
	    
	    unsigned long shift;
	    if (broadword::lsb(has_result, shift)) {
                uint8_t bit_pos = tables.m_fwd_pos[(word >> shift) & 0xFF][(cum_exc_step_8 >> shift) & 0xFF];
                assert(bit_pos > 0);
		ret = shift + bit_pos - 1;
		return true;
	    }
	    return false;
        }

        inline bool find_open_in_word(uint64_t word, uint64_t byte_counts, bp_vector::excess_t cur_exc, uint64_t& ret) {
            assert(cur_exc > 0 && cur_exc <= 64);
            const uint64_t rev_byte_counts = broadword::reverse_bytes(byte_counts);
            const uint64_t cum_exc_step_8 = (cur_exc - ((2 * rev_byte_counts - 8 * broadword::ones_step_8) << 8)) * broadword::ones_step_8;

	    uint64_t max_exc_step_8 = 0;
	    for (int i = 0; i < 8; ++i) {
		size_t shift = i * 8;
		max_exc_step_8 |= ((uint64_t)(tables.m_bwd_min[(word >> (64 - shift - 8)) & 0xFF])) << shift;
	    }

	    const uint64_t has_result = broadword::leq_step_8(cum_exc_step_8, max_exc_step_8);
	    
	    unsigned long shift;
	    if (broadword::lsb(has_result, shift)) {
                uint8_t bit_pos = tables.m_bwd_pos[(word >> (64 - shift - 8)) & 0xFF][(cum_exc_step_8 >> shift) & 0xFF];
                assert(bit_pos > 0);
		ret = 64 - (shift + bit_pos);
		return true;
	    }
	    return false;
        }
    }

    inline bool bp_vector::find_close_in_block(uint64_t block_offset, bp_vector::excess_t excess, uint64_t start, uint64_t& ret) const {
        if (excess > excess_t((bp_block_size - start) * 64)) {
            return false;
        }
        assert(excess > 0);
        for (uint64_t sub_block_offset = start; sub_block_offset < bp_block_size; ++sub_block_offset) {
            uint64_t sub_block = block_offset + sub_block_offset;
            uint64_t word = m_bits[sub_block];
            uint64_t byte_counts = broadword::byte_counts(word);
            assert(excess > 0);
            if (excess <= 64) {
                if (find_close_in_word(word, byte_counts, excess, ret)) {
                    ret += sub_block * 64;
                    return true;
                }
            }
            excess += static_cast<excess_t>(2 * broadword::bytes_sum(byte_counts) - 64);
        }
        return false;
    }
    
    uint64_t bp_vector::find_close(uint64_t pos) const
    {
        assert((*this)[pos]); // check there is an opening parenthesis in pos
        uint64_t ret = -1;
        // Search in current word
        uint64_t word_pos = (pos + 1) / 64;
        uint64_t shift = (pos + 1) % 64;
        uint64_t shifted_word = m_bits[word_pos] >> shift;
        // Pad with "open"
        uint64_t padded_word = shifted_word | (-!!shift & (~0ULL << (64 - shift)));
        uint64_t byte_counts = broadword::byte_counts(padded_word);

        excess_t word_exc = 1;
        if (find_close_in_word(padded_word, byte_counts, word_exc, ret)) {
            ret += pos + 1;
            return ret;
        }
        
        // Otherwise search in the local block
        uint64_t block = word_pos / bp_block_size;
        uint64_t block_offset = block * bp_block_size;
        uint64_t sub_block = word_pos % bp_block_size;
        uint64_t local_rank = broadword::bytes_sum(byte_counts) - shift; // subtract back the padding
        excess_t local_excess = static_cast<excess_t>((2 * local_rank) - (64 - shift));
        if (find_close_in_block(block_offset, local_excess + 1, sub_block + 1, ret)) {
            return ret;
        }

        // Otherwise, find the first appropriate block
        excess_t pos_excess = static_cast<excess_t>(2 * rank(pos) - pos);
        uint64_t found_block = search_min_tree<1>(block + 1, pos_excess);
        uint64_t found_block_offset = found_block * bp_block_size;
        excess_t found_block_excess = get_block_excess(found_block);

        // Search in the found block
        bool found = find_close_in_block(found_block_offset, found_block_excess - pos_excess, 0, ret);
        assert(found); (void)found;
        return ret;
    }

    inline bool bp_vector::find_open_in_block(uint64_t block_offset, bp_vector::excess_t excess, uint64_t start, uint64_t& ret) const {
        if (excess > excess_t(start * 64)) {
            return false;
        }
        assert(excess >= 0);
        uint64_t block = block_offset / bp_block_size;

        for (uint64_t sub_block_offset = start - 1; sub_block_offset + 1 > 0; --sub_block_offset) {
            assert(excess > 0);
            uint64_t sub_block = block_offset + sub_block_offset;
            uint64_t word = m_bits[sub_block];
            uint64_t byte_counts = broadword::byte_counts(word);
            if (excess <= 64) {
                if (find_open_in_word(word, byte_counts, excess, ret)) {
                    ret += sub_block * 64;
                    return true;
                }
            }
            excess -= static_cast<excess_t>(2 * broadword::bytes_sum(byte_counts) - 64);
        }
        return false;
    }

    uint64_t bp_vector::find_open(uint64_t pos) const
    {
        assert(pos);
        uint64_t ret = -1;
        // Search in current word
        uint64_t word_pos = (pos / 64);
        uint64_t len = pos % 64;
        // Rest is padded with "close"
        uint64_t shifted_word = -!!len & (m_bits[word_pos] << (64 - len));
        uint64_t byte_counts = broadword::byte_counts(shifted_word);

        excess_t word_exc = 1;
        if (find_open_in_word(shifted_word, byte_counts, word_exc, ret)) {
            ret += pos - 64;
            return ret;
        }

        // Otherwise search in the local block
        uint64_t block = word_pos / bp_block_size;
        uint64_t block_offset = block * bp_block_size;
        uint64_t sub_block = word_pos % bp_block_size;
        uint64_t local_rank = broadword::bytes_sum(byte_counts); // no need to subtract the padding
        excess_t local_excess = -static_cast<excess_t>((2 * local_rank) - len);
        if (find_open_in_block(block_offset, local_excess + 1, sub_block, ret)) {
            return ret;
        }

        // Otherwise, find the first appropriate block 
	excess_t pos_excess = static_cast<excess_t>(2 * rank(pos) - pos) - 1;
        uint64_t found_block = search_min_tree<0>(block - 1, pos_excess);
        uint64_t found_block_offset = found_block * bp_block_size;
        // Since search is backwards, have to add the current block
        excess_t found_block_excess = get_block_excess(found_block + 1);

        // Search in the found block
        bool found = find_open_in_block(found_block_offset, found_block_excess - pos_excess, bp_block_size, ret);
        assert(found); (void)found;
        return ret;
    }

    template <int direction>
    inline bool bp_vector::search_block_in_superblock(uint64_t block, excess_t excess, size_t& found_block) const 
    {
	size_t superblock = block / superblock_size;
        excess_t superblock_excess = get_block_excess(superblock * superblock_size);
	if (direction) {
	    for (size_t cur_block = block; 
		 cur_block < std::min((superblock + 1) * superblock_size, (size_t)m_block_excess_min.size()); 
		 ++cur_block) {
		if (excess >= superblock_excess + m_block_excess_min[cur_block]) {
		    found_block = cur_block;
		    return true;
		}
	    }
	} else {
	    for (size_t cur_block = block; 
		 cur_block + 1 >= (superblock * superblock_size) + 1;
		 --cur_block) {
		if (excess >= superblock_excess + m_block_excess_min[cur_block]) {
		    found_block = cur_block;
		    return true;
		}
	    }
	}

	return false;
    }

    inline bp_vector::excess_t bp_vector::get_block_excess(uint64_t block) const {
        uint64_t sub_block_idx = block * bp_block_size;
        uint64_t block_pos = sub_block_idx * 64;
        excess_t excess = static_cast<excess_t>(2 * sub_block_rank(sub_block_idx) - block_pos);
        assert(excess >= 0);
        return excess;
    }

    inline bool bp_vector::in_node_range(uint64_t node, excess_t excess) const {
        assert(m_superblock_excess_min[node] != size());
        return excess >= m_superblock_excess_min[node];
    }

    template <int direction>
    inline uint64_t bp_vector::search_min_tree(uint64_t block, excess_t excess) const 
    {
	size_t found_block;
	if (search_block_in_superblock<direction>(block, excess, found_block)) {
		return found_block;
	}
        
        int direction_sign = 2 * direction - 1;
        size_t cur_superblock = block / superblock_size;
	size_t cur_node = m_internal_nodes + cur_superblock;
        while (true) {
            assert(cur_node);
            bool going_back = (cur_node & 1) == direction;
            if (!going_back) {
                size_t next_node = cur_node + direction_sign;
                if (in_node_range(next_node, excess)) {
                    cur_node = next_node;
                    break;
                }
            }
            cur_node /= 2;
        }

	assert(cur_node);

        while (cur_node < m_internal_nodes) {
            uint64_t next_node = cur_node * 2 + (1 - direction);
            if (in_node_range(next_node, excess)) {
                cur_node = next_node;
                continue;
            }
            next_node += direction_sign;
            // if it is not one child, it must be the other
            assert(in_node_range(next_node, excess));
            cur_node = next_node;
        }
	
	size_t next_superblock = cur_node - m_internal_nodes;
	bool ret = search_block_in_superblock<direction>(next_superblock * superblock_size + (1 - direction) * (superblock_size - 1),
                                                         excess, found_block);
	assert(ret); (void)ret;

        return found_block;
    }

    
    uint64_t bp_vector::excess_rmq(uint64_t a, uint64_t b) const
    {
        return 0;
    }

    
    void bp_vector::build_min_tree() 
    {
	if (!size()) return;

        std::vector<block_min_excess_t> block_excess_min;
        excess_t cur_block_min = 0, cur_superblock_excess = 0;
        for (uint64_t sub_block = 0; sub_block < m_bits.size(); ++sub_block) {
            if (sub_block % bp_block_size == 0) {
                if (sub_block % (bp_block_size * superblock_size) == 0) {
                    cur_superblock_excess = 0;
                }
                if (sub_block) {
                    assert(cur_block_min >= std::numeric_limits<block_min_excess_t>::min());
                    assert(cur_block_min <= std::numeric_limits<block_min_excess_t>::max());
                    block_excess_min.push_back((block_min_excess_t)cur_block_min);
                    cur_block_min = cur_superblock_excess;
                }
            }
            uint64_t word = m_bits[sub_block];
            uint64_t mask = 1ULL;
            // for last block stop at bit boundary
            uint64_t n_bits = 
                (sub_block == m_bits.size() - 1 && size() % 64) 
                ? size() % 64 
                : 64;
	    // XXX(ot) use tables.m_fwd_{min,max}
            for (uint64_t i = 0; i < n_bits; ++i) {
                cur_superblock_excess += (word & mask) ? 1 : -1;
                cur_block_min = std::min(cur_block_min, cur_superblock_excess);
                mask <<= 1;
            }
        }
        // Flush last block mins
        assert(cur_block_min >= std::numeric_limits<block_min_excess_t>::min());
        assert(cur_block_min <= std::numeric_limits<block_min_excess_t>::max());
        block_excess_min.push_back((block_min_excess_t)cur_block_min);
        
        size_t n_blocks = util::ceil_div(data().size(), bp_block_size);
        assert(n_blocks == block_excess_min.size());

	size_t n_superblocks = (n_blocks + superblock_size - 1) / superblock_size;
	
        size_t n_complete_leaves = 1;
        while (n_complete_leaves < n_superblocks) n_complete_leaves <<= 1; // XXX(ot): I'm sure this can be done with broadword::msb...
	// n_complete_leaves is the smallest power of 2 >= n_superblocks
        m_internal_nodes = n_complete_leaves;
 	size_t treesize = m_internal_nodes + n_superblocks;
        
        std::vector<excess_t> superblock_excess_min(treesize);

        // Fill in the leaves of the tree
        for (size_t superblock = 0; superblock < n_superblocks; ++superblock) {
	    excess_t cur_super_min = static_cast<excess_t>(size());
            excess_t superblock_excess = get_block_excess(superblock * superblock_size);

	    for (size_t block = superblock * superblock_size; 
		 block < std::min((superblock + 1) * superblock_size, n_blocks);
		 ++block) {
                cur_super_min = std::min(cur_super_min, superblock_excess + block_excess_min[block]);
            }
            assert(cur_super_min >= 0 && cur_super_min < size());
             
            superblock_excess_min[m_internal_nodes + superblock] = cur_super_min;
        }
	
	// fill in the internal nodes with past-the-boundary values
	// (they will also serve as sentinels in debug)
	for (size_t node = 0; node < m_internal_nodes; ++node) {
	    superblock_excess_min[node] = static_cast<excess_t>(size());
	}

        // Fill bottom-up the other layers: each node updates the parent
        for (size_t node = treesize - 1; node > 1; --node) {
	    size_t parent = node / 2;
            superblock_excess_min[parent] = std::min(superblock_excess_min[parent], // same node
						     superblock_excess_min[node]);
        }

        m_block_excess_min.steal(block_excess_min);
        m_superblock_excess_min.steal(superblock_excess_min);
    }
}
