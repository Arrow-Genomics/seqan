// ==========================================================================
//                           arg_parse_exceptions.h
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: Stephan Aiche <stephan.aiche@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_EXCEPTIONS_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_EXCEPTIONS_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class ParseError
// ----------------------------------------------------------------------------

//TODO(singer): now in stream/tokinization.h
/*
.Internal.Class.ParseError
..cat:Miscellaneous
..summary:General ParseError.
*/

/*class ParseError : public RuntimeError
{
public:
    ParseError(std::string const & option) :
        RuntimeError(option)
    {}
};*/

// ----------------------------------------------------------------------------
// Class InvalidOption
// ----------------------------------------------------------------------------

/*
.Internal.Class.InvalidOption
..cat:Miscellaneous
..summary:Thrown if an unknown option was set on the command line.
*/

class InvalidOption : public ParseError
{
public:
    InvalidOption(std::string const & option) :
        ParseError("illegal option -- " + option)
    {}
};

// ----------------------------------------------------------------------------
// Class MissingArgument
// ----------------------------------------------------------------------------

/*
.Internal.Class.MissingArgument
..cat:Miscellaneous
..summary:Thrown if an option was set on the command line but without giving the
required arguments for this option.
*/

class MissingArgument : public ParseError
{
public:
    MissingArgument(std::string const & option) :
        ParseError("option requires an argument -- " + option)
    {}
};

// ----------------------------------------------------------------------------
// Class NotEnoughArguments
// ----------------------------------------------------------------------------

/*
.Internal.Class.NotEnoughArguments
..cat:Miscellaneous
..summary:Thrown if an option was set on the command line but not enough arguments for
this option were provided.
*/

class NotEnoughArguments : public ParseError
{
public:
    NotEnoughArguments(std::string const & option) :
        ParseError("option requires more arguments -- " + option)
    {}
};

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_EXCEPTIONS_H_
