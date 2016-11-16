#include <assert.h>
#include <iomanip>
#include "del.hpp"

MutatorDel::MutatorDel(uint ksize, uint delsize, Logger& l, uint maxabund, ulong lim)
    : Mutator(ksize, l, maxabund, lim), delsize(delsize)
{

}

ulong MutatorDel::process(std::string& sequence, Countgraph& countgraph)
{
    ulong kmercount = 0;
    for (ulong i = k - 1; i + k + delsize <= sequence.length(); i++) {
        if (limit > 0 && nuclcount > limit) {
            break;
        }
        nuclcount++;

        ulong min1 = i - k + 1;
        ulong length1 = k - 1;
        ulong min2 = i + delsize;
        ulong length2 = k;
std::cerr << "DEBUG " << min1 << ' ' << length1 << ' ' << min2 << ' ' << length2 << '\n';
        std::string delseq = sequence.substr(min1, length1) +
                             sequence.substr(min2, length2);
        assert(delseq.size() == (2*k) - 1);

        Deletion del(delseq, *this, countgraph);
        bool doprint = logger.increment();
        if (doprint) {
            float mb = (float)logger.counter / (float)1000000.0;
            logger.stream << std::setiosflags(std::ios::fixed)
                          << std::setprecision(1) << "# ...processed " << mb
                          << " Mb of sequence\n";
            del.print(logger.stream);
            logger.stream << "# " << abund_hist << "\n# " << unique_hist << '\n';
        }
        kmercount += k;
    }
    return kmercount;
}

MutatorDel::Deletion::Deletion(std::string& seq, MutatorDel& m, Countgraph& countgraph)
    : sequence(seq), mut(m)
{
    int unique_count = 0;
    std::vector<std::string> kmers;
    countgraph.get_kmers(sequence, kmers);
    assert(kmers.size() == mut.k);
    for (auto& kmer : kmers) {
        uint kmer_freq = countgraph.get_count(kmer.c_str());
        mut.abund_hist.increment(kmer_freq);
        abunds.push_back(kmer_freq);
        if (kmer_freq == 0) {
            unique_count++;
        }
    }
    mut.unique_hist.increment(unique_count);
}

void MutatorDel::Deletion::print(std::ostream& stream)
{
    stream << '>' << sequence << '\n';
    for (uint i = 0; i < mut.k; i++) {
        stream << ' ';
    }
    stream << '|' << '\n';
    for (uint i = 0; i < mut.k; i++) {
        uint abund = abunds[i];
        if (i > 0) {
            stream << ' ';
        }
        stream << abund;
    }
    stream << '\n';
}
