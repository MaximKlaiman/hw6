#ifndef HASH_H
#define HASH_H

#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
#include <string>
#include <cctype>

typedef std::size_t HASH_INDEX_T;

struct MyStringHash {
    HASH_INDEX_T rValues[5] { 983132572, 1468777056, 552714139, 984953261, 261934300 };
    MyStringHash(bool debug = true)
    {
        if(false == debug){
            generateRValues();
        }
    }
    // hash function entry point (i.e. this is h(k))
    HASH_INDEX_T operator()(const std::string& k) const
    {
        unsigned long long w[5] = {0ULL, 0ULL, 0ULL, 0ULL, 0ULL};

        int n = (int)k.size();
        int wi = 4;

        for(int end = n; end > 0 && wi >= 0; end -= 6)
        {
            int start = end - 6;
            if(start < 0) start = 0;

            unsigned long long val = 0ULL;

            // base-36 conversion (avoid pow): val = val*36 + digit
            for(int j = start; j < end; ++j)
            {
                char c = (char)std::tolower((unsigned char)k[j]);
                val = val * 36ULL + (unsigned long long)letterDigitToNumber(c);
            }

            w[wi] = val;   // w[4] is last chunk, then w[3], ...
            wi--;
        }

        unsigned long long h = 0ULL;
        for(int i = 0; i < 5; ++i)
        {
            h += w[i] * (unsigned long long)rValues[i];
        }
        return (HASH_INDEX_T)h;
    }

    // A likely helper function is to convert a-z,0-9 to an integral value 0-35
    HASH_INDEX_T letterDigitToNumber(char letter) const
    {
        if(letter >= 'a' && letter <= 'z') {
            return (HASH_INDEX_T)(letter - 'a');
        }
        if(letter >= '0' && letter <= '9') {
            return (HASH_INDEX_T)(26 + (letter - '0'));
        }
        return 0;
    }

    // Code to generate the random R values
    void generateRValues()
    {
        // obtain a seed from the system clock:
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::mt19937 generator (seed);  // mt19937 is a standard random number generator

        // Simply call generator() [it has an operator()] to get another random number
        for(int i{ 0 }; i < 5; ++i)
        {
            rValues[i] = generator();
        }
    }
};

#endif
