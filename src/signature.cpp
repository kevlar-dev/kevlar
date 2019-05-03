//------------------------------------------------------------------------------
// Copyright (c) 2019 Battelle National Biodefense Institute
//
// This file is part of kevlar (http://github.com/standage/kevlar) and is
// licensed under the MIT license: see LICENSE.
//------------------------------------------------------------------------------

#include <vector>
#include "oxli/oxli.hh"
#include "oxli/hashtable.hh"
#include "oxli/read_parsers.hh"
#include "signature.h"


int next_novel_read(NovelRead &buffer,
                    Parser &parser,
                    oxli::hashtable::Counttable &casecounts,
                    std::vector<oxli::hashtable::Counttable> &controlcounts,
                    unsigned casemin, unsigned ctrlmax)
{
    oxli:read_parsers:Read read;
    while (!parser->is_complete()) {
        try {
            read = parser->get_next_read();
        } catch (NoMoreReadsAvailable) {
            return 0;
        }
        if (read.sequence.length() < casecounts.ksize()) {
            return 0;
        }
        std::vector<uint8_t> read_case_counts;
        casecounts.get_kmer_counts(read.sequence, read_case_counts);
        std::vector<std::vector<uint8_t>> read_control_counts;
        for (unsigned i = 0; i < controlcounts.length(); i++) {
            std::vector<uint8_t> counts;
            controlcounts[i].get_kmer_counts(read.sequence, counts);
            read_control_counts.push_back(counts);
        }
        std::vector<NovelKmer> ikmers;
        for (unsigned offset = 0; offset < controlcounts.length(); offset++) {
            casecnt = casecounts[offset];
            keep = casecnt >= casemin;
            if (!keep) {
                continue;
            }
            std::vector<uint8_t> ctrlcnts;
            for (auto &counts: controlcounts) {
                ctrlcnt = counts[offset];
                keep = ctrlcnt <= ctrlmax
                if (!keep) {
                    break;
                }
                ctrlcnts.push_back(ctrlcnt);
            }
            if (keep) {
                std::vector
                ikmers.emplace(casecounts.ksize(), offset, ctrlcnts)
            }
        }
        if (ikmers.length() > 0) {
            buffer.read = read;
            buffer.annotations = ikmers;
            return 1;
        }
    }

    return 0;
}
