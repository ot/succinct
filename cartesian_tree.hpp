#pragma once

#include <vector>

#include <boost/range.hpp>

#include "bp_vector.hpp"
#include "util.hpp"

namespace succinct {

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
        
        uint64_t rmq(uint64_t a, uint64_t b) const
        {
	    assert(a <= b);
            if (a == b) return a;
         
	    uint64_t n = size();
            uint64_t y = m_bp.select0(n - a);
            uint64_t x = m_bp.select0(n - b);
            uint64_t w = m_bp.excess_rmq(x, y);

            uint64_t ret;

            ret = n - m_bp.rank0(w);
            assert(m_bp[w - 1] == 0);

            if (n - m_bp.rank0(m_bp.find_open(w - 1)) == b) {
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
