#include "common.inc"
int main(int argc, char** argv) {
    g_debug = std::getenv("SPC_DEBUG") != nullptr;
    g_t0 = std::chrono::steady_clock::now();
    read_graph(argv[1]); tlog("graph");
    read_queries(argv[2]); tlog("queries");
    write_answers(argv[3]); tlog("write");
    std::_Exit(0);
}
