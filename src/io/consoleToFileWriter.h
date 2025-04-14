#ifndef CONSOLETOFILEWRITER_H
#define CONSOLETOFILEWRITER_H

#include <iostream>
#include <fstream>
#include <streambuf>

namespace olb {

class DoubleBuffer : public std::streambuf {
public:
    DoubleBuffer(std::streambuf* consoleBuf, std::streambuf* fileBuf) : consoleBuffer(consoleBuf), fileBuffer(fileBuf) {}

protected:
    // This function is called for every character written to the stream
    virtual int overflow(int c) override {
        if (c != EOF) {
            // Write the character to both buffers
            if (consoleBuffer->sputc(c) == EOF || fileBuffer->sputc(c) == EOF) {
                return EOF;
            }
        }
        return c;
    }

    // Synchronize both buffers
    virtual int sync() override {
        if (consoleBuffer->pubsync() == 0 && fileBuffer->pubsync() == 0) {
            return 0;
        }
        return -1;
    }

private:
    std::streambuf* consoleBuffer; // Original console buffer
    std::streambuf* fileBuffer;    // File buffer
};

}

#endif
