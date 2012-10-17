#pragma once

#include <vector>

#include <boost/range.hpp>

#include "bp_vector.hpp"
#include "util.hpp"

namespace succinct {

    // This class implements a cartesian-tree-based RMQ data
    // structure, using the 2d-Min-Heap DFUDS representation described
    // in "Space-Efficient Preprocessing Schemes for Range Minimum
    // Queries on Static Arrays", Johannes Fischer and Volker Heun,
    // SIAM J. Comput., 40(2), 465–492.

    // We made a few variations:
    //
    // - The rmq() operation in the paper checks whether x is parent
    //   of w - 1, which can be written as select0(x - 1) <
    //   find_open(w - 1). We use instead the fact that the excess
    //   between x and w (both excluded) is strictly greater than the
    //   excess of w, so the formula above holds iff excess(select0(x
    //   - 1) + 1) <= excess(w). This is faster because a select0 is
    //   faster than find_open+rank0.
    //
    // - The construction is done in reverse order so that the input
    //   array can be traversed left-to-right. This involves
    //   re-mapping all the indices at query time
    //     
    // - Our data structures have 0-based indices, so the operations
    //   are slightly different from those in the paper

    class cartesian_tree : boost::noncopyable {
    public:

        cartesian_tree() {}
        
        template <typename Range>
        cartesian_tree(Range const& v)
        {
            typedef typename 
                boost::range_value<Range>::type value_type;
            typedef typename
                boost::range_iterator<const Range>::type iter_type;

            typedef std::pair<value_type, size_t> value_pos_type;
            std::vector<value_pos_type> s;
            size_t idx = 0;

	    uint64_t n = 2 * boost::size(v) + 2;
            bit_vector_builder bp(n);
	    uint64_t i = n - 1;
            
            for (iter_type it = boost::begin(v); it != boost::end(v); ++it) {
                value_pos_type cur(*it, idx++);
		i--; // prepend 0
                
                while (!s.empty() && s.back() > cur) {
                    s.pop_back();
                    bp.set(i--, 1); // prepend 1
                }
                
                s.push_back(cur);
            }
	    
	    i--; // fake root
            
            while (!s.empty()) {
                s.pop_back();
		bp.set(i--, 1); 
            }

	    bp.set(i--, 1); 

            bp_vector(&bp, false, true).swap(m_bp);
        }
        
        // NOTE: this is RMQ in the interval [a, b], b inclusive
        // XXX(ot): maybe change this to [a, b), for consistency with
        // the rest of the library?
        uint64_t rmq(uint64_t a, uint64_t b) const
        {
            typedef bp_vector::excess_t excess_t;

	    assert(a <= b);
            if (a == b) return a;
         
	    uint64_t n = size();

            uint64_t t = m_bp.select0(n - b - 1);
            excess_t exc_t = excess_t(t) - 2 * (n - b - 1);
            assert(exc_t - 1 == m_bp.excess(t + 1));

            uint64_t x = m_bp.select0(n - b);
            uint64_t y = m_bp.select0(n - a);

            excess_t exc_w;
            uint64_t w = m_bp.excess_rmq(x, y, exc_w);
            assert(m_bp[w - 1] == 0);

            uint64_t ret;
            if (exc_w >= exc_t - 1) {
                ret = b;
            } else {
                ret = n - m_bp.rank0(w);
            }

            assert(ret >= a);
            assert(ret <= b);
            return ret;
        }

	bp_vector const& get_bp() const
	{
	    return m_bp;
	}

	uint64_t size() const
	{
	    return m_bp.size() / 2 - 1;
	}
	
        template <typename Visitor>
        void map(Visitor& visit) 
        {
            visit
                (m_bp, "m_bp");
        }
	
        void swap(cartesian_tree& other) 
        {
            other.m_bp.swap(m_bp);
        }
        
    protected:
        bp_vector m_bp;
    };

}
