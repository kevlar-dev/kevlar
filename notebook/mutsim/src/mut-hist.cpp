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
    ulong memory;
    std::string seqfile;
    std::string refrfile;
    std::string refrcounts;
} ProgramArgs;

void print_usage(std::ostream& stream = std::cerr)
{
    stream << "Usage: mut-hist [options] genome.fa\n";
    stream << "  options:\n";
    stream << "    -c    load genome k-mer counts from the specified counttable\n";
    stream << "    -g    load genome k-mer counts from the specified seqfile\n";
    stream << "    -h    print this help message and exit\n";
    stream << "    -i    debug output interval (default: 1000000)\n";
    stream << "    -k    k-mer length (default: 31)\n";
    stream << "    -l    limit number of positions to process\n";
    stream << "    -m    max histogram value (default: 16)\n";
    stream << "    -r    sampling rate (default: 1.0)\n";
    stream << "    -s    seed for random number generator (default: 42)\n";
    stream << "    -t    mutation type: snv (default), del\n";
    stream << "    -y    memory consumption (in bytes; default: 2000000000)\n";
    stream << "    -z    deletion size (default: 5)\n";
}

void parse_args(int argc, const char **argv, ProgramArgs *args)
{
    char c;
    while ((c = getopt (argc, (char *const *)argv, "c:g:hi:k:l:m:r:s:t:y:z:")) != -1) {
        if (c == 'c') {
            args->refrcounts = optarg;
        }
        else if (c == 'g') {
            args->refrfile = optarg;
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
        else if (c == 'y') {
            args->memory = strtoul(optarg, NULL, 10);
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
    if (args->refrfile.compare("") == 0) {
        args->refrfile = args->seqfile;
    }
}

Counttable load_from_seqfile(std::string &seqfile, uint ksize, ulong memory)
{
    timepoint alloc_start = std::chrono::system_clock::now();
    std::cerr << "# allocating counttable...";
    std::vector<uint64_t> tablesizes = get_n_primes_near_x(4, memory / 4);
    Counttable counttable(ksize, tablesizes);
    timepoint alloc_end = std::chrono::system_clock::now();
    std::chrono::duration<double> alloc_elapsed = alloc_end - alloc_start;
    std::cerr << "done! (in " << alloc_elapsed.count() << " seconds)\n";

    std::cerr << "# consuming k-mers from seqfile " << seqfile << "...";
    timepoint consume_start = std::chrono::system_clock::now();
    unsigned int seqs_consumed = 0;
    unsigned long long kmers_consumed = 0;
    counttable.consume_seqfile<FastxReader>(seqfile, seqs_consumed, kmers_consumed);
    timepoint consume_end = std::chrono::system_clock::now();
    std::chrono::duration<double> consume_elapsed = consume_end - consume_start;
    std::cerr << "done! consumed " << seqs_consumed << " sequence(s) and "
              << kmers_consumed << " " << ksize << "-mers (in "
              << consume_elapsed.count() << " seconds)\n";

    return counttable;
}

Counttable load_from_counttable(std::string &counttablefile)
{
    timepoint alloc_start = std::chrono::system_clock::now();
    std::cerr << "# loading reference genome k-mers counts from " << counttablefile << "...";
    std::vector<uint64_t> tablesizes = {1};
    Counttable counttable(1, tablesizes);
    counttable.load(counttablefile);
    timepoint alloc_end = std::chrono::system_clock::now();
    std::chrono::duration<double> alloc_elapsed = alloc_end - alloc_start;
    std::cerr << "done! (in " << alloc_elapsed.count() << " seconds)\n";    
    return counttable;
}

int main(int argc, const char **argv)
{
    if (argc == 1) {
        print_usage(std::cout);
        return 0;
    }

    ProgramArgs args = {5, 1000000, 31, 0, 16, 1.0, 42, "snv", 2000000000, "", "", ""};
    parse_args(argc, argv, &args);

    Counttable counttable = args.refrcounts.compare("") == 0 
                              ? load_from_counttable(args.refrcounts)
                              : load_from_seqfile(args.refrfile, args.ksize, args.memory);

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
