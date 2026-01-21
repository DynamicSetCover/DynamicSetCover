#ifndef DYNSC_CTRL_H
#define DYNSC_CTRL_H

#include <vector>
#include "../DataTypes.h"
#include <chrono>
#include "alg1.h"
#include "alg2.h"
#include "alg3.h"
#include "alg4.h"


//   - Constructor: AlgorithmType(int n, int m, double eps)
//   - Methods: updateElement(const std::vector<int>&), getRecourse(), getSetCoverSize(), getSetCoverRef()
template<typename AlgorithmType>
inline std::vector<UpdateResult> run_algorithm_loop(const std::vector<std::vector<int>>& update_sequence,
                                                    int k,
                                                    int n,
                                                    int m,
                                                    double eps,
                                                    bool store_set_covers = true) {
    std::vector<UpdateResult> results;
    results.reserve(k); // Pre-allocate memory
    
    AlgorithmType alg(n, m, eps);
    
    // Iterate over all k updates
    for (int t = 0; t < k && t < static_cast<int>(update_sequence.size()); ++t) {
        // Start chrono timer for this update step
        auto start = std::chrono::steady_clock::now();
        
        // Extract update information from the sequence
        const auto& update = update_sequence[t];
        if (update.empty()) {
            // For empty updates, still record a result with 0s
            if (store_set_covers) {
                results.push_back({0, 0, alg.getSetCoverSize(), alg.getSetCoverRef()});
            } else {
                results.push_back({0, 0, alg.getSetCoverSize(), std::vector<int>()});
            }
            continue;
        }
        
        // Call element update on the algorithm
        alg.updateElement(update);

        // Stop timer
        auto stop = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
        
        // Fetch recourse, set cover size, and set cover for this update (O(1))
        int update_recourse = alg.getRecourse();
        int set_cover_size = alg.getSetCoverSize();
        
        // Store results for this update
        if (store_set_covers) {
            const std::vector<int>& set_cover_ref = alg.getSetCoverRef(); // O(1), no copy, 0-indexed
            results.push_back({duration.count(), update_recourse, set_cover_size, set_cover_ref});
        } else {
            results.push_back({duration.count(), update_recourse, set_cover_size, std::vector<int>()});
        }
    }
    
    return results;
}

#endif // DYNSC_CTRL_H