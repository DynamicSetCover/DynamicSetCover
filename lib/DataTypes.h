#ifndef DYNSC_DATATYPES_H
#define DYNSC_DATATYPES_H

#include <vector>

// Structure to store results for each update
struct UpdateResult {
    long long time_nanoseconds;  // Time taken for the update
    int recourse;                // Recourse for this update
    int set_cover_size;          // Set cover size after the update
    std::vector<int> set_cover;   // Actual set cover (0-indexed)
};

#endif // DYNSC_DATATYPES_H