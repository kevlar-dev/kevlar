#include <assert.h>
#include <iomanip>
#include "snv.hpp"

MutatorSNV::MutatorSNV(uint ksize, Logger& l, uint maxabund, ulong lim)
    : Mutator(ksize, l, maxabund, lim), nuclcount(0), nucl("ACGT")
{

}

ulong MutatorSNV::process(std::string& sequence, Countgraph& countgraph)
{
    ulong kmercount = 0;
    for (ulong i = k - 1; i + k <= sequence.length(); i++) {
        if (limit > 0 && nuclcount > limit) {
            break;
        }
        if (skip_nucl()) {
            continue;
        }
        nuclcount++;

        ulong minindex = i - k + 1;
        uint subseqlength = (2*k) - 1;
        std::string snv_interval = sequence.substr(minindex, subseqlength);

        SingleNucleotideVariant snv(snv_interval, *this, countgraph);
        bool doprint = logger.increment();
        if (doprint) {
            float mb = (float)logger.counter / (float)1000000.0;
            logger.stream << std::setiosflags(std::ios::fixed)
                          << std::setprecision(1) << "# ...processed " << mb
                          << " Mb of sequence\n";
            snv.print(logger.stream);
            logger.stream << "# " << abund_hist << "\n# " << unique_hist << '\n';
        }
        kmercount += k * 3;
    }
    return kmercount;
}

MutatorSNV::SingleNucleotideVariant::SingleNucleotideVariant(std::string& seq, MutatorSNV& m, Countgraph& countgraph)
    : sequence(seq), mut(m)
{
    for (auto bp : mut.nucl) {
        if (bp == sequence[mut.k - 1]) {
            continue;
        }

        std::string mutseq = sequence;
        mutseq[mut.k - 1] = bp;

        int unique_count = 0;
        std::vector<std::string> kmers;
        countgraph.get_kmers(mutseq, kmers);
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
}

void MutatorSNV::SingleNucleotideVariant::print(std::ostream& stream)
{
    stream << '>' << sequence << '\n';
    for (uint i = 0; i < mut.k; i++) {
        stream << ' ';
    }
    stream << '|' << '\n';
    for (uint i = 0; i < mut.k * 3; i++) {
        uint abund = abunds[i];
        if (i % mut.k != 0) {
            stream << ' ';
        }
        stream << abund;
        if ((i + 1) % mut.k == 0) {
            stream << '\n';
        }
    }
}

ulong MutatorSNV::get_mut_count()
{
    return nuclcount * 3;
}
