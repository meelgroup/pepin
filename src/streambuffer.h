/*******************************************************************************
MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*******************************************************************************/

#ifndef STREAMBUFFER_H
#define STREAMBUFFER_H

static const unsigned chunk_limit = 148576;

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>
#include <memory>
#include <cmath>

#ifdef USE_ZLIB
#include <zlib.h>
struct GZ {
    static inline int read(void* buf, size_t num, size_t count, gzFile f)
    {
        return gzread(f, buf, num*count);
    }
};
#endif

struct FN {
    static inline int read(void* buf, size_t num, size_t count, FILE* f)
    {
        return fread(buf, num, count, f);
    }
};

struct CH {
    static inline int read(void* buf, size_t num, size_t count, const char*& f)
    {
        int toread = num*count;
        char* mybuf = (char*)buf;

        int read = 0;
        while(*f != 0 && read < toread) {
            *mybuf = *f;
            mybuf++;
            f++;
            read++;
        }
        return read;
    }
};

template<typename A, typename B>
class StreamBuffer
{
    A  in;
    void assureLookahead() {
        if (pos >= size) {
            pos  = 0;
            size = B::read(buf.get(), 1, chunk_limit, in);
        }
    }
    int     pos;
    int     size;
    std::unique_ptr<char[]> buf;

    void advance()
    {
        operator++();
    }
    int value()
    {
        return operator*();
    }

public:
    StreamBuffer(A i) :
        in(i)
        , pos(0)
        , size(0)
        , buf(new char[chunk_limit]())
    {
        assureLookahead();
    }

    int  operator *  () {
        return (pos >= size) ? EOF : buf[pos];
    }
    void operator ++ () {
        pos++;
        assureLookahead();
    }

    void skipWhitespace()
    {
        char c = value();
        while (c == '\t' || c == '\r' || c == ' ') {
            advance();
            c = value();
        }
    }

    void skipLine()
    {
        for (;;) {
            if (value() == EOF || value() == '\0') return;
            if (value() == '\n') {
                advance();
                return;
            }
            advance();
        }
    }

    bool skipEOL(const size_t lineNum)
    {
        for (;;) {
            if (value() == EOF || value() == '\0') return true;
            if (value() == '\n') {
                advance();
                return true;
            }
            if (value() != '\r') {
                std::cerr
                << "PARSE ERROR! Unexpected char (hex: " << std::hex
                << std::setw(2)
                << std::setfill('0')
                << "0x" << value()
                << std::setfill(' ')
                << std::dec
                << ")"
                << " At line " << lineNum+1
                << " we expected an end of line character (\\n or \\r + \\n)"
                << std::endl;
                return false;
            }
            advance();
        }
        exit(-1);
    }

    inline bool parseRatio(uint32_t& dividend, uint32_t& divisor, size_t lineNum)
    {
        int32_t head;
        bool rc = parseInt(head, lineNum);
        if (!rc) {
            return false;
        }
        if (head < 0) {
            std::cerr << "ERROR: dividend is negative!" << std::endl;
            exit(-1);
        }
        dividend = head;
        if (value() == '/') {
            advance();
            int64_t tail;
            rc = parseInt<int64_t>(tail, lineNum);
            if (!rc) {
                return false;
            }
            if (tail < 0) {
                std::cerr << "ERROR: divisor is negative!" << std::endl;
                exit(-1);
            }
            divisor = tail;
        } else {
            std::cerr << "ERROR: Only dividend was provided!" << std::endl;
            exit(-1);
        }
        return true;
    }

    template<class T=int32_t>
    inline bool parseInt(T& ret, size_t lineNum, bool allow_eol = false)
    {
        T val = 0;
        T mult = 1;
        skipWhitespace();
        if (value() == '-') {
            mult = -1;
            advance();
        } else if (value() == '+') {
            advance();
        }

        char c = value();
        if (allow_eol && c == '\n') {
            ret = std::numeric_limits<T>::max();
            return true;
        }
        if (c < '0' || c > '9') {
            std::cerr
            << "PARSE ERROR! Unexpected char (dec: '" << c << ")"
            << " At line " << lineNum
            << " we expected a number"
            << std::endl;
            return false;
        }

        while (c >= '0' && c <= '9') {
            T val2 = val*10 + (c - '0');
            if (val2 < val) {
                std::cerr << "PARSE ERROR! At line " << lineNum
                << " the variable number is to high"
                << std::endl;
                return false;
            }
            val = val2;
            advance();
            c = value();
        }
        ret = mult*val;
        return true;
    }

    void parseString(std::string& str)
    {
        str.clear();
        skipWhitespace();
        while (value() != ' ' && value() != '\n' && value() != EOF) {
            str.push_back(value());
            advance();
        }
    }
};

#endif //STREAMBUFFER_H
