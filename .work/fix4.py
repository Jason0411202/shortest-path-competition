s = open('hy.cpp', 'rb').read().decode()
def rep(old, new, cnt=1):
    global s
    assert s.count(old) == cnt, old
    s = s.replace(old, new)

rep("            auto& target = fwd ? ci : co;", "            auto& target = (fwd && directed) ? ci : co;")
rep('''        auto& cin = directed ? ci[i] : co[i];
        LI.start[h] = LI.hub.size();
        for (auto& e : cin) { LI.hub.push_back(e.first); LI.dist.push_back(e.second); }
        LI.len[h] = (uint32_t)cin.size();''', '''        if (directed) {
            LI.start[h] = LI.hub.size();
            for (auto& e : ci[i]) { LI.hub.push_back(e.first); LI.dist.push_back(e.second); }
            LI.len[h] = (uint32_t)ci[i].size();
        }''')
rep('''            if (!directed && !fwd) {
                LI.start[hv] = LO.start[hv]; LI.len[hv] = LO.len[hv];
                break;
            }
            Labels& L = fwd ? LO : LI;
            const Labels& other = fwd ? LI : LO;''', '''            if (!directed && !fwd) break;
            Labels& L = fwd ? LO : LI;
            const Labels& other = (fwd && directed) ? LI : LO;''')
rep('''    const uint32_t* hb = LI.hub.data() + LI.start[b];
    const uint32_t* eb = hb + LI.len[b];
    const W* db = LI.dist.data() + LI.start[b];''', '''    const Labels& R = directed ? LI : LO;
    const uint32_t* hb = R.hub.data() + R.start[b];
    const uint32_t* eb = hb + R.len[b];
    const W* db = R.dist.data() + R.start[b];''')
open('hy.cpp', 'wb').write(s.encode())
