#include "ex.h"
#include "optional.hpp"

#include <vector>
#include <map>
#include <iostream>


namespace GiNaC {

struct CMatcher;
using opt_exmap = nonstd::optional<exmap>;
using opt_bool = nonstd::optional<bool>;
using opt_CMatcher = nonstd::optional<CMatcher>;
using uvec = std::vector<size_t>;
using nonstd::nullopt;

struct CMatcher_state {
};

struct CMatcher {
        CMatcher(const ex &source_, const ex & pattern_, exmap& map_)
         : source(source_), pattern(pattern_), map(map_)
        {
                ret_val = init();
                if (ret_val and not ret_val.value()) {
                        finished = true;
                        ret_map.reset();
                }
        }
        opt_bool init();
        void run();
        void noncomm_run();
        void no_global_wild();
        void with_global_wild();
        void perm_run(const exvector&, const exvector&);
        void comb_run(const exvector&, const exvector&);
        opt_exmap get()
        {
                if (ret_val) {
                        if (not ret_val.value())
                                return nullopt;
                        ret_val.reset();
                        return ret_map;
                }
                ret_map.reset();
                ++level;
                run();    // guarantees to set ret, and if true, map
                --level;
                return ret_map;
        }
        bool get_alt(size_t);
        void clear_ret()
        {
                ret_val.reset();
                ret_map.reset();
        }

        static int level;
        ex source, pattern;
        opt_bool ret_val;
        opt_exmap ret_map;
        exmap map;
        size_t N{0}, P{0}, wi{0};
        exvector ops, pat, gws, gwp;
        std::vector<opt_CMatcher> cms;
        std::vector<exmap> map_repo;
        std::vector<bool> m_cmatch;
        bool finished{false};
        uvec perm, comb, wild_ind;
        enum Type { unset, comm, noncomm, comm_plus };
        Type type {unset};
};

}
