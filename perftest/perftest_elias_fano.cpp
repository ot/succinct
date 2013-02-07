#include <iostream>
#include <vector>

#include <boost/lexical_cast.hpp>

#include "util.hpp"

#include "elias_fano.hpp"
#include "mapper.hpp"

#include "perftest_common.hpp"

void ef_enumeration_benchmark(uint32_t m, uint32_t max_value)
{
    std::vector<uint32_t> v(m);
    for (size_t i = 0; i < m; ++i) { 
	v[i] = rand() % max_value;
    }

    std::sort(v.begin(), v.end());
    succinct::elias_fano::elias_fano_builder bvb(max_value, m);
    for (size_t i = 0; i < m; ++i) { 
	bvb.push_back(v[i]);
    }

    succinct::elias_fano ef(&bvb);
    
    double elapsed;
    uint32_t foo;
    SUCCINCT_TIMEIT(elapsed) {
	succinct::elias_fano::select_enumerator it(ef, 0);
	for (size_t i = 0; i < m; ++i) {
	    foo += it.next();
	}
    }
    volatile uint32_t vfoo = foo;
    std::cerr << "Elapsed: " << elapsed / 1000 << " msec\n"
	      << m / elapsed << " Mcodes/s" << std::endl;
}

int main(int argc, char** argv)
{
    size_t m = boost::lexical_cast<uint32_t>(argv[1]);
    size_t max_value = boost::lexical_cast<uint32_t>(argv[2]);
        
    ef_enumeration_benchmark(m, max_value);
}
