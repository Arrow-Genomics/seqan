// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Code for writing SAM.
// ==========================================================================

#ifndef INCLUDE_SEQAN_BAM_IO_WRITE_SAM_H_
#define INCLUDE_SEQAN_BAM_IO_WRITE_SAM_H_

namespace seqan {

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function write()                                             BamHeaderRecord
// ----------------------------------------------------------------------------

template <typename TTarget, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void write(TTarget & target,
                  BamHeaderRecord const & header,
                  BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> const & /*context*/,
                  Sam const & /*tag*/)
{
    char const * headerTypes[] = {"@HD", "@SQ", "@RG", "@PG", "@CO"};
    write(target, headerTypes[header.type]);

    if (header.type == BAM_HEADER_COMMENT && !empty(header.tags))
    {
        writeValue(target, '\t');
        write(target, header.tags[0].i2);
    }
    else
    {
        for (unsigned i = 0; i < length(header.tags); ++i)
        {
            writeValue(target, '\t');
            write(target, header.tags[i].i1);
            writeValue(target, ':');
            write(target, header.tags[i].i2);
        }
    }

    writeValue(target, '\n');
}

// ----------------------------------------------------------------------------
// Function write()                                                   BamHeader
// ----------------------------------------------------------------------------

template <typename TTarget, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void write(TTarget & target,
                  BamHeader const & header,
                  BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> const & context,
                  Sam const & tag)
{
    String<bool> writtenSeqInfos;
    resize(writtenSeqInfos, length(contigNames(context)), false);

    size_t globalRefId = 0;
    for (unsigned i = 0; i < length(header); ++i)
    {
        BamHeaderRecord const & record = header[i];
        if (record.type == BAM_HEADER_REFERENCE)
            for (unsigned j = 0; j < length(record.tags); ++j)
                if (record.tags[j].i1 == "SN")
                {
                    if (getIdByName(globalRefId, contigNamesCache(context), record.tags[j].i2))
                        writtenSeqInfos[globalRefId] = true;
                    break;
                }

        write(target, record, context, tag);
    }

    // Write missing @SQ header records.
    SEQAN_ASSERT_LEQ(length(contigLengths(context)), length(contigNames(context)));
    for (unsigned i = 0; i < length(contigLengths(context)); ++i)
    {
        if (writtenSeqInfos[i])
            continue;
        write(target, "@SQ\tSN:");
        write(target, contigNames(context)[i]);
        write(target, "\tLN:");
        appendNumber(target, contigLengths(context)[i]);
        writeValue(target, '\n');
    }

}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_BAM_IO_WRITE_SAM_H_
