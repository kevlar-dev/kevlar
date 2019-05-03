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


typedef oxli::read_parsers::ReadParserPtr<oxli.read_parsers:FastxReader> Parser;

typedef struct
{
    unsigned ksize;
    unsigned offset;
    std::vector<uint8_t> abunds;
} NovelKmer;


typedef struct
{
    oxli::read_parsers::Read read;
    std::vector<NovelKmer> annotations;
} NovelRead;


int next_novel_read(NovelRead &buffer,
                    Parser &parser,
                    oxli::hashtable::Counttable &casecounts,
                    std::vector<oxli::hashtable::Counttable> &controlcounts,
                    unsigned casemin, unsigned ctrlmax);
