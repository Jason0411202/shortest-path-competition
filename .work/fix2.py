s = open('sp.cpp', 'rb').read().decode()
old = "            bool stalled = false;\n            for (uint32_t k = H[u]; k < H[u + 1]; ++k) {"
assert s.count(old) == 2
s = s.replace(old, "            bool stalled = false;\n            if (STALL) for (uint32_t k = H[u]; k < H[u + 1]; ++k) {")
s = s.replace('static int USE_PLL = 0;', 'static int USE_PLL = 0;\nstatic int STALL = 1;')
s = s.replace('    if (std::getenv("SHUF")) SHUF = 1;', '    if (std::getenv("SHUF")) SHUF = 1;\n    if (const char* e = std::getenv("STALL")) STALL = std::atoi(e);')
open('sp.cpp', 'wb').write(s.encode())
