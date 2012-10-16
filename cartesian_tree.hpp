#pragma once

#include <stack>

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
                boost::range_reverse_iterator<const Range>::type riter_type;

            // XXX int
            typedef std::pair<value_type, int> value_pos_type;
            std::stack<value_pos_type> s;

            std::vector<bool> bp;

            int p = 0;
            
            for (riter_type it = boost::rbegin(v); it != boost::rend(v); ++it) {
                value_pos_type cur(*it, p--);
                bp.push_back(0);
                
                while (!s.empty() && s.top() > cur) {
                    s.pop();
                    bp.push_back(1);
                }
                
                s.push(cur);
            }
	    
	    bp.push_back(0); // fake root	    
            
            while (!s.empty()) {
                s.pop();
                bp.push_back(1);
            }

	    bp.push_back(1); // fake root	    

            std::reverse(bp.begin(), bp.end());
               
            // XXX reverse select support?
            bp_vector(bp, false, true).swap(m_bp);
        }
        
        uint64_t rmq(uint64_t a, uint64_t b) const
        {
            assert(a <= b);
            if (a == b) return a;
            uint64_t x = m_bp.select0(a + 1);
            uint64_t y = m_bp.select0(b + 1);
            uint64_t w = m_bp.excess_rmq(x, y);

            uint64_t ret;

            ret = m_bp.rank0(w);
            assert(m_bp[w - 1] == 0);

            if (m_bp.rank0(m_bp.find_open(w - 1)) == a + 1) {
                ret = a;
            } else {
                ret = m_bp.rank0(w) - 1;
            }
            
            assert(ret >= a);
            assert(ret <= b);
            return ret;
        }

	bp_vector const& get_bp() const
	{
	    return m_bp;
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
