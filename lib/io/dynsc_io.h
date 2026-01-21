#ifndef DYNSC_IO_H
#define DYNSC_IO_H

#include <string>
#include <vector>
#include <tuple>

class dynsc_io {
public:
    dynsc_io();
    
    // Read update sequence from file
    // Returns tuple (k, m, n) where k is number of updates, m is number of sets, n is max number of elements
    std::tuple<int, int, int> readUpdateSequence(const std::string& filename, std::vector<std::vector<int>>& update_sequence);
};

#endif // DYNSC_IO_H

