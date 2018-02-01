#include <chrono>
#include <getopt.h>
#include <memory>
#include "snv.hpp"
#include "del.hpp"
#include "read_parsers.hh"

using namespace oxli;
using namespace read_parsers;
typedef Read Sequence;
typedef std::chrono::time_point<std::chrono::system_clock> timepoint;

typedef struct
{
    uint delsize;
    ulong interval;
    uint ksize;
    ulong limit;
    uint histmax;
    double sampling_rate;
    int seed;
    std::string muttype;
    std::string seqfile;
    std::string refrcounts;
} ProgramArgs;

void print_usage(std::ostream& stream = std::cerr)
{
    stream << "Usage: snv-hist [options] genome.fa kmercounts.counttable\n";
    stream << "  options:\n";
    stream << "    -h    print this help message and exit\n";
    stream << "    -i    debug output interval (default: 1000000)\n";
    stream << "    -k    k-mer length (default: 31)\n";
    stream << "    -l    limit number of positions to process\n";
    stream << "    -m    max histogram value (default: 16)\n";
    stream << "    -r    sampling rate (default: 1.0)\n";
    stream << "    -s    seed for random number generator (default: 42)\n";
    stream << "    -t    mutation type: snv (default), del\n";
    stream << "    -z    deletion size (default: 5)\n";
}

void parse_args(int argc, const char **argv, ProgramArgs *args)
{
    char c;
    while ((c = getopt (argc, (char *const *)argv, "d:hi:k:l:m:r:s:t:y:z:")) != -1) {
        if (c == 'd') {
            args->delsize = atoi(optarg);
        }
        else if (c == 'h') {
            print_usage(std::cout);
            exit(0);
        }
        else if (c == 'i') {
            args->interval = atoi(optarg);
        }
        else if (c == 'k') {
            args->ksize = atoi(optarg);
        }
        else if (c == 'l') {
            args->limit = atoi(optarg);
        }
        else if (c == 'm') {
            args->histmax = atoi(optarg);
        }
        else if (c == 'r') {
            args->sampling_rate = atof(optarg);
        }
        else if (c == 's') {
            args->seed = atoi(optarg);
        }
        else if (c == 't') {
            args->muttype = optarg;
        }
        else if (c == 'z') {
            args->delsize = atoi(optarg);
        }
        else {
            std::cerr << "Unknown option '" << c << "'\n";
            print_usage();
            exit(1);
        }
    }

    args->seqfile = argv[optind];
    args->refrcounts = argv[optind + 1];
}

int main(int argc, const char **argv)
{
    if (argc == 1) {
        print_usage(std::cout);
        return 0;
    }

    ProgramArgs args = {5, 1000000, 31, 0, 16, 1.0, 42, "snv", "", ""};
    parse_args(argc, argv, &args);

    timepoint alloc_start = std::chrono::system_clock::now();
    std::cerr << "# loading reference genome k-mers counts from " << args.refrcounts << "...";
    std::vector<uint64_t> tablesizes = {1};
    Counttable counttable(1, tablesizes);
    counttable.load(args.refrcounts);
    timepoint alloc_end = std::chrono::system_clock::now();
    std::chrono::duration<double> alloc_elapsed = alloc_end - alloc_start;
    std::cerr << "done! (in " << alloc_elapsed.count() << " seconds)\n";

    timepoint query_start = std::chrono::system_clock::now();
    std::cerr << "# iterating through " << args.seqfile << ", querying k-mer abundances...";
    Logger logger(args.interval, std::cerr);
    std::unique_ptr<Mutator> mut = NULL;
    if(args.muttype == "snv") {
        mut.reset(new MutatorSNV(args.ksize, logger, args.histmax, args.limit));
    }
    else if(args.muttype == "del") {
        mut.reset(new MutatorDel(args.ksize, args.delsize, logger, args.histmax, args.limit));
    }
    else {
        std::cerr << "Error: unknown mutation type '" << args.muttype << "'\n";
        exit(1);
    }
    mut->set_sampling_rate(args.sampling_rate, args.seed);

    FastxParserPtr parser = get_parser<FastxReader>(args.seqfile);
    Sequence seq;
    while (!parser->is_complete()) {
        try {
            seq = parser->get_next_read();
        } catch (NoMoreReadsAvailable &e) {
            break;
        }
        mut->process(seq.sequence, counttable);
    }
    timepoint query_end = std::chrono::system_clock::now();
    std::chrono::duration<double> query_elapsed = query_end - query_start;

    std::cerr << "done! processed " << mut->get_mut_count()
              << " mutations (in " << query_elapsed.count() << " seconds)\n\n";
    std::cout << *mut;

    return 0;
}
