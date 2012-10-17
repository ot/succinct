#pragma once

#include <vector>
#include <queue>

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "cartesian_tree.hpp"

namespace succinct {

    // XXX(ot): write a version that enumerates results one by one
    // XXX(ot): implement arbitrary comparator
    template <typename Vector>
    class topk_vector : boost::noncopyable {
    public:
	typedef Vector vector_type;
	typedef typename vector_type::value_type value_type;
	typedef boost::tuple<value_type, uint64_t> entry_type;
	
        template <typename Range>
        topk_vector(Range const& v)
	{
	    cartesian_tree(v, std::greater<typename boost::range_value<Range>::type>())
		.swap(m_cartesian_tree);
	    vector_type(v).swap(m_v);
	}
	
	value_type const
	operator[](uint64_t idx) const
	{
	    return m_v[idx];
	}

	uint64_t size() const
	{
	    return m_v.size();
	}

	std::vector<entry_type> 
	topk(uint64_t a, uint64_t b, size_t k) const
	{
	    assert(a <= b);

	    using boost::get;
	    using boost::tie;

	    std::vector<entry_type> ret(std::min(size_t(b - a + 1), k));
	    std::priority_queue<queue_element_type, 
				std::vector<queue_element_type>, 
				value_index_comparator> q;
	    
	    uint64_t m = m_cartesian_tree.rmq(a, b);
	    q.push(queue_element_type(m_v[m], m, a, b));

	    for (size_t i = 0; i < ret.size() - 1; ++i) {
		value_type cur_mid_val;
		uint64_t cur_mid, cur_a, cur_b;
		tie(cur_mid_val, cur_mid, cur_a, cur_b) = q.top(); 
		q.pop();
		
		ret[i] = entry_type(cur_mid_val, cur_mid);
		
		if (cur_mid != cur_a) {
		    m = m_cartesian_tree.rmq(cur_a, cur_mid - 1);
		    q.push(queue_element_type(m_v[m], m, cur_a, cur_mid - 1));
		}

		if (cur_mid != cur_b) {
		    m = m_cartesian_tree.rmq(cur_mid + 1, cur_b);
		    q.push(queue_element_type(m_v[m], m, cur_mid + 1, cur_b));
		}
	    }

	    ret.back() = entry_type(get<0>(q.top()), get<1>(q.top()));
	    
	    return ret;
	}


        template <typename Visitor>
        void map(Visitor& visit) 
        {
            visit
                (m_v, "m_v");
                (m_cartesian_tree, "m_cartesian_tree");
        }
	
        void swap(topk_vector& other) 
        {
            other.m_v.swap(m_v);
            other.m_cartesian_tree.swap(m_cartesian_tree);
        }

    private:
	
	typedef boost::tuple<value_type, uint64_t, uint64_t, uint64_t> queue_element_type;

	struct value_index_comparator {
	    template <typename Tuple>
	    bool operator()(Tuple const& a, Tuple const& b) const 
	    {
		using boost::get;
		// lexicographic, increasing on value and decreasing
		// on index
		return (get<0>(a) < get<0>(b) ||
			(get<0>(a) == get<0>(b) &&
			 get<1>(a) > get<1>(b)));
	    }
	};
	
	vector_type m_v;
	cartesian_tree m_cartesian_tree;
    };

}
