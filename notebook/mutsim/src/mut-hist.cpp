#include <getopt.h>
#include <memory>
#include "primes.hpp"
#include "snv.hpp"
#include "del.hpp"
#include "read_parsers.hh"

using namespace khmer;
using namespace read_parsers;
typedef Read Sequence;

typedef struct
{
    uint delsize;
    ulong interval;
    uint ksize;
    ulong limit;
    uint histmax;
    uint numtables;
    std::string muttype;
    ulong targetsize;
    std::string infile;
    std::string refrfile;
} ProgramArgs;

void print_usage(std::ostream& stream = std::cerr)
{
    stream << "Usage: snv-hist [options] seqs.fa [genome.fa]\n";
    stream << "  options:\n";
    stream << "    -h    print this help message and exit\n";
    stream << "    -i    debug output interval (default: 1000000)\n";
    stream << "    -k    k-mer length (default: 31)\n";
    stream << "    -l    limit number of positions to process\n";
    stream << "    -m    max histogram value (default: 16)\n";
    stream << "    -n    num tables (default: 4)\n";
    stream << "    -t    mutation type: snv (default), del\n";
    stream << "    -x    approx table size (default: 5e8)\n";
    stream << "    -z    deletion size (default: 5)\n";
}

void parse_args(int argc, const char **argv, ProgramArgs *args)
{
    char c;
    while ((c = getopt (argc, (char *const *)argv, "d:hi:k:l:m:n:t:x:z:")) != -1) {
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
        else if (c == 'n') {
            args->numtables = atoi(optarg);
        }
        else if (c == 't') {
            args->muttype = optarg;
        }
        else if (c == 'x') {
            args->targetsize = atoi(optarg);
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

    args->infile = argv[optind];
    args->refrfile = argv[optind];
    if (argc > optind + 1) {
        args->refrfile = argv[optind + 1];
    }
}

int main(int argc, const char **argv)
{
    if (argc == 1) {
        print_usage(std::cout);
        return 0;
    }

    ProgramArgs args = {5, 1000000, 31, 0, 16, 4, "snv", 500000000, "", ""};
    parse_args(argc, argv, &args);

    std::cerr << "# allocating countgraph\n";
    std::vector<HashIntoType> tablesizes = get_n_primes_near_x(args.targetsize, args.numtables);
    Countgraph countgraph(args.ksize, tablesizes);

    std::cerr << "# consuming reference...";
    unsigned int seqs_consumed = 0;
    unsigned long long kmers_consumed = 0;
    countgraph.consume_fasta(args.refrfile, seqs_consumed, kmers_consumed);
    std::cerr << "consumed " << seqs_consumed << " sequence(s) and "
              << kmers_consumed << " " << args.ksize << "-mers\n";

    std::cerr << "# querying k-mer abundance...\n\n";
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

    IParser *parser = IParser::get_parser(args.infile);
    Sequence seq;
    while (!parser->is_complete()) {
        try {
            seq = parser->get_next_read();
        } catch (NoMoreReadsAvailable &e) {
            break;
        }
        mut->process(seq.sequence, countgraph);
    }
    std::cout << *mut;
    delete parser;

    return 0;
}
