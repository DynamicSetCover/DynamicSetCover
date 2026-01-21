#include "alg3.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <sstream>
#include <string>
#include <functional>

// Structure to represent an update operation
struct Update {
    int operation;  // 0 = insert, 1 = delete
    int element_id; // element identifier
    std::vector<int> sets; // sets containing the element (only for insertions)
};

// Anonymous namespace for helper functions (to avoid multiple definition errors)
namespace {
    // Helper: floor log base beta for n >= 1
    inline int floor_log_beta(int n, double beta) {
        if (n <= 1) return 0;
        return static_cast<int>(std::log(n) / std::log(beta)); // safe for n>=1
    }

    // Helper: Calculate which bucket a set belongs to based on uncovered count
    inline int getBucket(int uncovered_count, const std::vector<int>& bucket_boundaries) {
        if (uncovered_count == 0) return -1; // No bucket for sets with 0 uncovered
        // Linear search from highest bucket 
        for (int j = bucket_boundaries.size() - 1; j >= 0; j--) {
            if (uncovered_count >= bucket_boundaries[j]) {
                return j;
            }
        }
        return 0; // Should not reach here, but return 0 as fallback
    }
}

// Constructor: initialize data structures
Algorithm3::Algorithm3(int n, int m, double eps)
    : n(n), m(m), eps(eps), beta(1 + eps), sc_size(0), last_update_recourse(0), last_set_cover_size(0) {
    // Calculate number of levels
    num_levels = floor_log_beta(n, beta) + 1;
    
    // Initialize data structures
    element_levels.resize(n, -1);
    set_levels.resize(m, -1);
    cov_sets.resize(m);
    cov_sets_size.resize(m, 0);
    NjS.resize(m, std::vector<std::vector<int>>(num_levels));
    NjS_size.resize(m, std::vector<int>(num_levels, 0));
    level_to_elements.resize(num_levels + 1);
    S.resize(num_levels, 0);
    D.resize(num_levels, 0.0);
    set_cover.clear(); // Initially empty set cover
}

// Process a single update from update sequence
void Algorithm3::updateElement(const std::vector<int>& update_row) {
    if (update_row.empty()) {
        return; // Skip empty updates
    }
    
    // Convert update_row to Update struct
    Update update;
    update.operation = update_row[0];
    update.element_id = update_row[1];
    
    // For insertions, extract sets (convert from 1-indexed to 0-indexed)
    if (update.operation == 0 && update_row.size() > 2) {
        for (size_t i = 2; i < update_row.size(); ++i) {
            update.sets.push_back(update_row[i] - 1); // Convert to 0-indexed
        }
    }
    
    // Call applyUpdate with the Update struct
    applyUpdate(update);
}

// Add set s to set_cover vector (maintains O(1) access)
void Algorithm3::addSetToCover(int s) {
    // Check if set is already in the cover vector
    for (int set_id : set_cover) {
        if (set_id == s) {
            return; // Already in cover
        }
    }
    set_cover.push_back(s);
}

// Remove set s from set_cover vector (maintains O(1) access)
void Algorithm3::removeSetFromCover(int s) {
    // Find and remove set s from set_cover using swap-and-pop
    for (size_t i = 0; i < set_cover.size(); ++i) {
        if (set_cover[i] == s) {
            // Swap with last element
            set_cover[i] = set_cover.back();
            set_cover.pop_back();
            return;
        }
    }
}


// Anonymous namespace for helper functions (to avoid multiple definition errors)
namespace {
    // Helper: Encode (s, j, e) into a single integer key for the hash map
    // Assumes: s < m, j < num_levels, e < n
    inline long long encodeKey(int s, int j, int e, int num_levels, int n) {
        return (long long)s * num_levels * n + (long long)j * n + e;
    }

    // Helper: Delete element e from NjS[s][j] in O(1) time using swap-and-pop
    // Requires: NjS_index must contain the key for (s, j, e)
    // Returns: true if element was found and deleted, false otherwise
    bool deleteFromNjS(int s, int j, int e,
                        int num_levels, int n,
                        std::vector<std::vector<std::vector<int>>>& NjS,
                        std::vector<std::vector<int>>& NjS_size,
                    std::unordered_map<long long, int>& NjS_index) {
        // Bounds checking
        if (s < 0 || s >= (int)NjS.size()) {
            return false; // Invalid set index
        }
        if (j < 0 || j >= (int)NjS[s].size()) {
            return false; // Invalid level index
        }
        
        long long key = encodeKey(s, j, e, num_levels, n);
        
        // Check if element exists in the index map
        if (NjS_index.find(key) == NjS_index.end()) {
            return false; // Element not found
        }
        
        int idx = NjS_index[key]; // Get the index of element e
        
        // Check if vector is empty or index is out of bounds
        if (NjS[s][j].empty() || idx < 0 || idx >= (int)NjS[s][j].size()) {
            // Inconsistent state - clean up and return
            NjS_index.erase(key);
            return false;
        }
        
        int last_idx = NjS[s][j].size() - 1;
        
        if (idx != last_idx) {
            // Swap with last element
            int last_elem = NjS[s][j][last_idx];
            NjS[s][j][idx] = last_elem;
            // Update index map for the swapped element
            long long swapped_key = encodeKey(s, j, last_elem, num_levels, n);
            NjS_index[swapped_key] = idx;
        }
        
        // Remove from vector (pop last element)
        NjS[s][j].pop_back();
        // Remove from index map
        NjS_index.erase(key);
        // Decrement size
        if (s < (int)NjS_size.size() && j < (int)NjS_size[s].size()) {
            NjS_size[s][j]--;
        }
        
        return true;
    }

    // Helper: Add element e to NjS[s][j] and maintain index map
    void addToNjS(int s, int j, int e,
               int num_levels, int n,
               std::vector<std::vector<std::vector<int>>>& NjS,
               std::vector<std::vector<int>>& NjS_size,
                    std::unordered_map<long long, int>& NjS_index) {
        // Bounds checking
        if (s < 0 || s >= (int)NjS.size()) {
            return; // Invalid set index
        }
        if (j < 0 || j >= num_levels) {
            return; // Invalid level index
        }
        
        // Ensure NjS[s] has enough levels
        if (j >= (int)NjS[s].size()) {
            NjS[s].resize(j + 1);
        }
        // Ensure NjS_size[s] has enough levels
        if (s >= (int)NjS_size.size()) {
            NjS_size.resize(s + 1);
        }
        if (j >= (int)NjS_size[s].size()) {
            NjS_size[s].resize(j + 1, 0);
        }
        
        int new_idx = NjS[s][j].size();
        NjS[s][j].push_back(e);
        long long key = encodeKey(s, j, e, num_levels, n);
        NjS_index[key] = new_idx;
        NjS_size[s][j]++;
    }

    // Helper: Delete element e from cov_sets[s] in O(1) time using swap-and-pop
    // Returns: true if element was found and deleted, false otherwise
    bool deleteFromCovSet(int e,
                       std::vector<std::vector<int>>& cov_sets,
                       std::vector<int>& cov_sets_size,
                       std::unordered_map<int, std::pair<int, int>>& cov_sets_index) {
        // Check if element exists in the index map
        if (cov_sets_index.find(e) == cov_sets_index.end()) {
            return false; // Element not found
        }
        
        auto [s, idx] = cov_sets_index[e]; // Get the set and index of element e
        int last_idx = cov_sets[s].size() - 1;
        
        if (idx != last_idx) {
            // Swap with last element
            int last_elem = cov_sets[s][last_idx];
            cov_sets[s][idx] = last_elem;
            // Update index map for the swapped element
            cov_sets_index[last_elem] = {s, idx};
        }
        
        // Remove from vector (pop last element)
        cov_sets[s].pop_back();
        // Remove from index map
        cov_sets_index.erase(e);
        // Decrement size
        cov_sets_size[s]--;
        
        return true;
    }

    // Helper: Add element e to cov_sets[s] and maintain index map
    void addToCovSet(int s, int e,
                 std::vector<std::vector<int>>& cov_sets,
                 std::vector<int>& cov_sets_size,
                    std::unordered_map<int, std::pair<int, int>>& cov_sets_index) {
        int new_idx = cov_sets[s].size();
        cov_sets[s].push_back(e);
        cov_sets_index[e] = {s, new_idx};
        cov_sets_size[s]++;
    }

    // Helper: Remove set s from cover if its cover set is empty
    // Sets set_levels[s] = -1 and decrements sc_size
    void removeSetFromCoverIfEmpty(int s,
                                std::vector<int>& S,
                                std::vector<int>& set_levels,
                                std::vector<int>& cov_sets_size,
                                int& sc_size,
                                std::function<void(int)>& trackSetChange) {
        // Bounds check
        if (s < 0 || s >= (int)set_levels.size() || s >= (int)cov_sets_size.size()) {
            return;
        }
        
        if (set_levels[s] >= 0 && cov_sets_size[s] == 0) {
            int old_level = set_levels[s];
            trackSetChange(s); // Track this set change
            set_levels[s] = -1;
            sc_size --;
            if (old_level >= 0 && old_level < (int)S.size()) {
                S[old_level]--;
            }
        }
    }

    // Helper: Get the set and index for element e in O(1) time
    // Returns: pair (s, index) if found, (-1, -1) otherwise
    std::pair<int, int> getCovSetPosition(int e,
                                  const std::unordered_map<int, std::pair<int, int>>& cov_sets_index) {
        auto it = cov_sets_index.find(e);
        if (it == cov_sets_index.end()) {
            return {-1, -1}; // Element not found
        }
        return it->second;
    }

    // Helper: Remove element e from level_to_elements at its current level
    void removeFromLevelToElements(int e,
                                int current_level,
                                int num_levels,
                                std::vector<std::vector<int>>& level_to_elements) {
        if (current_level < 0) {
            // Level -1 stored at index num_levels
            std::vector<int>& vec = level_to_elements[num_levels];
            vec.erase(std::remove(vec.begin(), vec.end(), e), vec.end());
        } else if (current_level < num_levels) {
            std::vector<int>& vec = level_to_elements[current_level];
            vec.erase(std::remove(vec.begin(), vec.end(), e), vec.end());
        }
    }

    // Helper: Add element e to level_to_elements at its new level
    void addToLevelToElements(int e,
                           int new_level,
                           int num_levels,
                           std::vector<std::vector<int>>& level_to_elements) {
        if (new_level < 0) {
            // Level -1 stored at index num_levels
            level_to_elements[num_levels].push_back(e);
        } else if (new_level < num_levels) {
            level_to_elements[new_level].push_back(e);
        }
    }

    // Helper: Update level_to_elements when element moves from old_level to new_level
    void updateLevelToElements(int e,
                            int old_level,
                            int new_level,
                            int num_levels,
                            std::vector<std::vector<int>>& level_to_elements) {
        if (old_level != new_level) {
            removeFromLevelToElements(e, old_level, num_levels, level_to_elements);
            addToLevelToElements(e, new_level, num_levels, level_to_elements);
        }
    }


    // Check if any set containing element e is PD (Positive Dirty)
    // A set s is PD with respect to level j if:
    //   sum(NjS_size[s][i] for i=0 to j-1) >= beta^(j+1)
    // Checks all sets containing e, from high level to low level
    // For each level j from num_levels down to min_level:
    //   - Checks all sets s containing e if they are PD w.r.t. level j
    // Returns: pair (set, level) if PD set found, (-1, -1) otherwise
    // Time Complexity: O(f * (num_levels - min_level)) where f is max frequency
    std::pair<int, int> isSetPD(int e,
                       int min_level, int num_levels,
                       const std::vector<std::vector<int>>& adj,
                       const std::vector<std::vector<int>>& NjS_size,
                       double beta) {
        // Precompute prefix sums for each set containing e
        // For each set s, compute sum for j = num_levels first, then decrement
        std::vector<int> prefix_sums(adj[e].size());
        
        // Initialize prefix sums for j = num_levels (sum of i=0 to num_levels-1)
        for (size_t idx = 0; idx < adj[e].size(); ++idx) {
            int s = adj[e][idx];
            int sum = 0;
            for (int i = 0; i < num_levels; ++i) {
                if (i < (int)NjS_size[s].size()) {
                    sum += NjS_size[s][i];
                }
            }
            prefix_sums[idx] = sum;
        }
        
        // Check levels from high to low
        for (int j = num_levels; j >= min_level; --j) {
            // For each set s containing element e
            for (size_t idx = 0; idx < adj[e].size(); ++idx) {
                int s = adj[e][idx];
                
                // Update prefix sum: subtract NjS_size[s][j] if j < num_levels
                // (for j = num_levels, we already have the full sum)
                if (j < num_levels && j < (int)NjS_size[s].size()) {
                    prefix_sums[idx] -= NjS_size[s][j];
                }
                
                // Calculate beta^(j+1)
                double threshold = std::pow(beta, j + 1);
                
                // Check if set s is PD with respect to level j
                if (prefix_sums[idx] >= threshold) {
                    return {s, j}; // Return the first PD set and level found
                }
            }
        }
        
        return {-1, -1}; // No set containing e is PD
    }


    // Function to perform greedy MSC on participating elements and sets
    // Called by reset(k) to run static greedy MSC algorithm on the subset
    void msc(int n,
         int m,
         double beta,
         int num_levels,
         const std::vector<std::vector<int>>& adj,
         std::vector<int>& element_levels,
         std::vector<int>& set_levels,
         int& sc_size,
         std::vector<std::vector<int>>& cov_sets,
         std::vector<int>& cov_sets_size,
         std::unordered_map<int, std::pair<int, int>>& cov_sets_index,
         int max_truncate_level,
         const std::unordered_set<int>& participating_elements,
         const std::unordered_set<int>& participating_sets,
         std::vector<std::vector<int>>* level_to_elements)
    {
        // Reset participating elements and sets
        for (int e : participating_elements) {
            if (e < (int)element_levels.size()) {
                int old_level = element_levels[e];
                element_levels[e] = -1;
                if (level_to_elements != nullptr) {
                    removeFromLevelToElements(e, old_level, num_levels, *level_to_elements);
                    addToLevelToElements(e, -1, num_levels, *level_to_elements);
                }
            }
        }
        
        for (int s : participating_sets) {
            set_levels[s] = -1;
        }
        
        // Remove participating elements from cov_sets_index and clear their cov_sets
        for (int e : participating_elements) {
            cov_sets_index.erase(e);
        }
        
        for (int s : participating_sets) {
            cov_sets[s].clear();
            cov_sets_size[s] = 0;
        }

        // Build set_lst and num_uncvrd for participating elements/sets
        std::unordered_map<int, int> num_uncvrd;         // number of currently uncovered elements in each set
        std::unordered_map<int, std::vector<int>> set_lst;     // elements in each set (static)
        std::unordered_set<int> cvrd;                     // elements that are covered

        for (int e : participating_elements) {
            if (e < (int)adj.size()) {
                for (int s : adj[e]) {
                    if (participating_sets.find(s) != participating_sets.end()) {
                        set_lst[s].push_back(e);
                        num_uncvrd[s]++;   // initially all elements are uncovered
                    }
                }
            }
        }

        // Find maximum uncovered count to determine number of buckets needed
        int max_uncvrd = 0;
        for (const auto& [s, count] : num_uncvrd) {
            if (participating_sets.find(s) != participating_sets.end() && count > max_uncvrd) {
                max_uncvrd = count;
            }
        }
        
        // Calculate number of buckets: log_{beta}(max_uncvrd) = log(max_uncvrd) / log(beta)
        int num_buckets = (max_uncvrd > 0) ? (int)std::ceil(std::log(max_uncvrd) / std::log(beta)) + 1 : 1;
        
        // Precompute bucket boundaries: bucket_boundaries[j] = beta^j
        // Use iterative multiplication to avoid expensive pow() calls
        std::vector<int> bucket_boundaries(num_buckets);
        bucket_boundaries[0] = 1;
        double current = 1.0;
        for (int j = 1; j < num_buckets; j++) {
            current *= beta;
            bucket_boundaries[j] = (int)std::floor(current);
            if (bucket_boundaries[j] <= bucket_boundaries[j-1]) {
                bucket_boundaries[j] = bucket_boundaries[j-1] + 1; // Ensure strictly increasing
            }
        }
        
        // Initialize buckets: buckets[j] contains set indices in bucket j
        std::vector<std::vector<int>> buckets(num_buckets);
        for (const auto& [s, count] : num_uncvrd) {
            if (participating_sets.find(s) != participating_sets.end() && count > 0) {
                int bucket_idx = getBucket(count, bucket_boundaries);
                if (bucket_idx >= 0 && bucket_idx < num_buckets) {
                    buckets[bucket_idx].push_back(s);
                }
            }
        }

        // Find highest non-empty bucket
        int highest_bucket = -1;
        for (int j = num_buckets - 1; j >= 0; j--) {
            if (!buckets[j].empty()) {
                highest_bucket = j;
                break;
            }
        }

        int uncovered = participating_elements.size(); // All participating elements start uncovered

        // Main greedy loop with bucket-based lazy update
        while (uncovered > 0) {
            // Find highest non-empty bucket with a valid set (lazy update)
            int s = -1;
            int cnt = 0;
            while (highest_bucket >= 0) {
                // Clean up invalid sets from this bucket (sets that moved to different buckets)
                while (!buckets[highest_bucket].empty()) {
                    int candidate = buckets[highest_bucket].back();
                    buckets[highest_bucket].pop_back();
                    
                    // Check if this set is still valid (has uncovered elements, in correct bucket)
                    auto it = num_uncvrd.find(candidate);
                    if (it != num_uncvrd.end() && it->second > 0) {
                        int correct_bucket = getBucket(it->second, bucket_boundaries);
                        if (correct_bucket == highest_bucket) {
                            // This set is valid and in the correct bucket
                            s = candidate;
                            cnt = it->second;
                            break;
                        } else if (correct_bucket >= 0 && correct_bucket < num_buckets) {
                            // Move to correct bucket
                            buckets[correct_bucket].push_back(candidate);
                        }
                    }
                }
                
                if (s >= 0) {
                    break; // Found a valid set
                }
                
                // This bucket is now empty, move to next
                highest_bucket--;
            }
            
            if (s < 0 || cnt <= 0) {
                break; // No more sets available
            }
            
            // choose level for this set
            int chosen_level = floor_log_beta(cnt, beta);
            // Apply truncation if specified
            if (max_truncate_level >= 0 && chosen_level > max_truncate_level) {
                chosen_level = max_truncate_level;
            }
            if (chosen_level >= num_levels) chosen_level = num_levels - 1;
            
            set_levels[s] = chosen_level;
            
            // iterate over that set's elements; only handle ones not yet covered
            auto set_lst_it = set_lst.find(s);
            if (set_lst_it == set_lst.end()) continue;
            for (int e : set_lst_it->second) {
                if (cvrd.find(e) != cvrd.end()) continue;    // skip elements already covered by earlier chosen sets.

                // mark element covered and assign its level
                cvrd.insert(e);
                
                // Update level_to_elements if provided
                if (level_to_elements != nullptr) {
                    int old_elem_level = (e < (int)element_levels.size()) ? element_levels[e] : -1;
                    updateLevelToElements(e, old_elem_level, chosen_level, num_levels, *level_to_elements);
                }
                
                if (e >= (int)element_levels.size()) {
                    element_levels.resize(e + 1, -1);
                }
                element_levels[e] = chosen_level;

                // add to cov set of s
                addToCovSet(s, e, cov_sets, cov_sets_size, cov_sets_index);

                // for every participating set that contains e, decrement uncovered counters
                // Sets will be moved to correct buckets lazily when we encounter them
                for (int s_prime : adj[e]) {
                    if (participating_sets.find(s_prime) == participating_sets.end()) {
                        continue;
                    }

                    // decrement uncovered count for s_prime
                    auto it = num_uncvrd.find(s_prime);
                    if (it != num_uncvrd.end() && it->second > 0) {
                        it->second--;
                    }
                }

                uncovered--;
            } // end for each element in set s
            
            // Update highest_bucket if current bucket is now empty
            while (highest_bucket >= 0 && buckets[highest_bucket].empty()) {
                highest_bucket--;
            }
        } // end while uncovered > 0
    }


    // Check if invariant is violated: sum of all D entries >= eps/beta * |S|
    // |S| is the number of sets at level >= 0 (sc_size)
    bool checkDirtInvariant(int num_levels,
                         double eps,
                         double beta,
                         const std::vector<double>& D,
                         int sc_size) {
        double sum_D = 0.0;
        for (int i = 0; i < num_levels; ++i) {
            sum_D += D[i];
        }
        
        double threshold = (eps / beta) * sc_size;
        
        //cout << "Total D is " << sum_D << endl;
        //cout << "Total S is " << sc_size << endl;
        
        return (sum_D >= threshold);
    }

    // Find the highest critical level
    // A level i is critical if for all j in [0, i]: D[j] + ... + D[i] >= eps/beta * (S[j] + ... + S[i])
    // Check in O(num_levels^2) time, starting from num_levels-1 down to 0
    int findHighestCriticalLevel(int num_levels,
                              double eps,
                              double beta,
                              const std::vector<double>& D,
                              const std::vector<int>& S) {
        int critical_level = -1;
        
        // Check from highest level down to 0
        for (int i = num_levels - 1; i >= 0; --i) {
            bool is_critical = true;
            
            // Check for all j in [0, i]
            for (int j = 0; j <= i; ++j) {
                // Compute sum D[j] + ... + D[i]
                double sum_D = 0.0;
                for (int k = j; k <= i; ++k) {
                    sum_D += D[k];
                }
                
                // Compute sum S[j] + ... + S[i]
                int sum_S = 0;
                for (int k = j; k <= i; ++k) {
                    sum_S += S[k];
                }
                
                // Check if D[j]+...+D[i] >= eps/beta * (S[j]+...+S[i])
                double threshold = (eps / beta) * sum_S;
                if (sum_D < threshold) {
                    is_critical = false;
                    break;
                }
            }
            
            if (is_critical) {
                critical_level = i;
                break;
            }
        }
        
        return critical_level;
    }


    // Reset procedure: run msc on elements at level <=k and sets containing them
    // Truncates levels to k+1
    // Resets D counters for k' <= k (zeros D[0] through D[k])
    // Uses level_indexed data structures to optimize to O(n_k) instead of O(n)
    void resetAtLevel(int k,
                  int n,
                  int m,
                  double beta,
                  int num_levels,
                  const std::vector<std::vector<int>>& adj,
                  std::vector<int>& element_levels,
                  std::vector<int>& set_levels,
                  int& sc_size,
                  std::vector<std::vector<int>>& cov_sets,
                  std::vector<int>& cov_sets_size,
                  std::unordered_map<int, std::pair<int, int>>& cov_sets_index,
                  std::vector<int>& S,
                  std::vector<double>& D,
                  std::vector<std::vector<std::vector<int>>>& NjS,
                  std::vector<std::vector<int>>& NjS_size,
                  std::unordered_map<long long, int>& NjS_index,
                  std::vector<std::vector<int>>& level_to_elements,
                  std::function<void(int)>& trackSetChange) {
        double sum_D = 0.0;
        int sum_S = 0;
        for (int i = 0 ; i < num_levels ; i++) {
            sum_D += D[i];
            sum_S += S[i];
        }

        // Collect elements at level <= k using level_to_elements for O(n_k) efficiency
        std::unordered_set<int> participating_elements;
        for (int lev = 0; lev <= k && lev < num_levels; ++lev) {
            for (int e : level_to_elements[lev]) {
                if (e < (int)adj.size()) { // Ensure element still exists
                    participating_elements.insert(e);
                }
            }
        }

        // Collect sets at level <= k that contain at least one participating element. Takes time O(n_k * f).
        std::unordered_set<int> participating_sets;
        for (int lev = 0; lev <= k && lev < num_levels; ++lev) {
            for (int e : level_to_elements[lev]) {
                if (e < (int)adj.size()) {
                    for (int s : adj[e]) {
                        // Only mark sets at level <= k
                        int set_lev = (s < (int)set_levels.size()) ? set_levels[s] : -1;
                        if (set_lev <= k) {
                            participating_sets.insert(s);
                        }
                    }
                }
            }
        }
        
        // Store old element levels before msc (needed to remove from NjS)
        std::unordered_map<int, int> old_element_levels;
        for (int e : participating_elements) {
            if (e < (int)element_levels.size()) {
                old_element_levels[e] = element_levels[e];
            }
        }
        
        // Store previous set_levels for participating sets to calculate recourse in O(|participating_sets|)
        std::unordered_map<int, int> prev_set_levels_map;
        for (int s : participating_sets) {
            if (s < (int)set_levels.size()) {
                prev_set_levels_map[s] = set_levels[s];
                // Add all participating sets to tracking, regardless of whether their level changes
                trackSetChange(s);
            }
        }
        
        // Run msc without truncation, passing level_indexed structure for maintenance
        msc(n, m, beta, num_levels, adj, element_levels, set_levels, sc_size,
            cov_sets, cov_sets_size, cov_sets_index,
            -1, participating_elements, participating_sets,
            &level_to_elements);

        // Rebuild NjS data structures to match new element levels
        // Only need to update participating elements (non-participating elements' levels didn't change)
        // Ensure NjS and NjS_size are properly sized
        if (NjS.size() < m) {
            NjS.resize(m);
        }
        if (NjS_size.size() < m) {
            NjS_size.resize(m);
        }
        
        // Update NjS for participating elements
        for (int e : participating_elements) {
            if (e >= (int)adj.size()) continue;
            
            // Cache adj[e] to avoid repeated vector lookups
            const std::vector<int>& sets_containing_e = adj[e];
            
            // Get old and new levels once
            auto it_old = old_element_levels.find(e);
            int old_level = (it_old != old_element_levels.end()) ? it_old->second : -1;
            int new_level = (e < (int)element_levels.size()) ? element_levels[e] : -1;
            
            // Remove from old level
            if (old_level >= 0 && old_level < num_levels) {
                for (int s : sets_containing_e) {
                    deleteFromNjS(s, old_level, e, num_levels, n, NjS, NjS_size, NjS_index);
                }
            }
            
            // Add to new level
            if (new_level >= 0 && new_level < num_levels) {
                for (int s : sets_containing_e) {
                    addToNjS(s, new_level, e, num_levels, n, NjS, NjS_size, NjS_index);
                }
            }
        }

        // Update S counter incrementally for participating sets only
        // Subtract old contributions and add new contributions
        for (int s : participating_sets) {
            if (s < (int)set_levels.size()) {
                // Subtract old level contribution
                auto it = prev_set_levels_map.find(s);
                if (it != prev_set_levels_map.end()) {
                    int old_level = it->second;
                    if (old_level >= 0 && old_level < num_levels) {
                        S[old_level]--;
                    }
                }
                // Add new level contribution
                int new_level = set_levels[s];
                if (new_level >= 0 && new_level < num_levels) {
                    S[new_level]++;
                }
            }
        }
        
        // Recompute sc_size to match sum of all S[j] (includes sets at level > k)
        sc_size = 0;
        for (int i = 0; i < num_levels; ++i) {
            sc_size += S[i];
        }
        
        // Reset D counters for k' <= k: zero D[0] through D[k]
        for (int k_prime = 0; k_prime <= k && k_prime < num_levels; ++k_prime) {
            D[k_prime] = 0.0;
        }
    }


    // Rise function: Raise the level of set s to j+1
    // Takes all elements in set s at level less than j and creates a new cov set at level j+1
    // Updates set_levels[s] and element_levels[e] for all such elements
    // Updates NjS and cov_sets structures accordingly
    // Adds dirt when elements leave their previous sets: 1/beta^j for each element at level j
    void rise(int s, int j,
          int num_levels, int n, int m, double beta, double eps,
          std::vector<int>& set_levels,
          std::vector<int>& element_levels,
          int& sc_size,
          std::vector<int>& S,
          std::vector<double>& D,
          std::vector<std::vector<int>>& cov_sets,
          std::vector<int>& cov_sets_size,
          const std::vector<std::vector<int>>& adj,
          std::vector<std::vector<std::vector<int>>>& NjS,
          std::vector<std::vector<int>>& NjS_size,
          std::unordered_map<long long, int>& NjS_index,
          std::unordered_map<int, std::pair<int, int>>& cov_sets_index,
          std::vector<std::vector<int>>& level_to_elements,
          std::function<void(int)>& trackSetChange) {
        int new_level = j + 1;

        // Clamp to valid range
        if (new_level >= num_levels) {
            new_level = num_levels - 1;
        }
        
        // Collect all elements in set s at level less than j
        std::vector<int> elements_to_rise;
        for (int i = 0; i < j; ++i) {
            if (i < (int)NjS[s].size()) {
                // Add all elements at level i in set s
                for (int e : NjS[s][i]) {
                    elements_to_rise.push_back(e);
                }
            }
        }

        // Bounds check for set index
        if (s < 0 || s >= (int)set_levels.size()) {
            return; // Invalid set index
        }
        
        // Get old set level for S counter update
        int old_set_level = set_levels[s];
        
        // Update set level
        trackSetChange(s); // Track this set change
        set_levels[s] = new_level;
        
        // Update S counter: decrement old level, increment new level
        if (old_set_level >= 0 && old_set_level < num_levels) {
            S[old_set_level]--;
            sc_size --;
        }
        if (new_level >= 0 && new_level < num_levels) {
            S[new_level]++;
            sc_size ++;
        }
        
        // Process each element to rise
        for (int e : elements_to_rise) {
            // Bounds check for element index
            if (e < 0 || e >= (int)adj.size()) {
                continue; // Skip invalid element index
            }
            
            // Get the element's old level (should be >= 0 since e is in NjS[s][i])
            int old_level = (e < (int)element_levels.size()) ? element_levels[e] : -1;
            
            // Add dirt: element at level old_level leaves its previous set
            // Dirt = 1/beta^old_level
            if (old_level >= 0 && old_level < num_levels && old_level < (int)D.size()) {
                D[old_level] += 1.0 / std::pow(beta, old_level);
            }
            
            auto old_pos = getCovSetPosition(e, cov_sets_index);
            if (old_pos.first != -1 && old_pos.first != s) {
                // Element was covered by a different set, remove it
                deleteFromCovSet(e, cov_sets, cov_sets_size, cov_sets_index);
                
                // Check if the old covering set's cover set is now empty
                removeSetFromCoverIfEmpty(old_pos.first, S, set_levels, cov_sets_size, sc_size, trackSetChange);
            } else if (old_pos.first == s) {
                // Element was already in cov_sets[s], just remove it (we'll re-add it below)
                deleteFromCovSet(e, cov_sets, cov_sets_size, cov_sets_index);
            }
            
            // Remove element from NjS[s_prime][old_level] for all sets s_prime containing e
            for (int s_prime : adj[e]) {
                deleteFromNjS(s_prime, old_level, e, num_levels, n, NjS, NjS_size, NjS_index);
            }
            
            // Add element to cov_sets[s] at new level
            addToCovSet(s, e, cov_sets, cov_sets_size, cov_sets_index);
            
            // Update element level
            if (e >= (int)element_levels.size()) {
                element_levels.resize(e + 1, -1);
            }
            element_levels[e] = new_level;
            
            // Update level_to_elements
            updateLevelToElements(e, old_level, new_level, num_levels, level_to_elements);
            
            // Add element to NjS[s_prime][new_level] for all sets s_prime containing e
            for (int s_prime : adj[e]) {
                // Ensure NjS[s_prime] has enough levels
                if (new_level >= (int)NjS[s_prime].size()) {
                    NjS[s_prime].resize(new_level + 1);
                    NjS_size[s_prime].resize(new_level + 1, 0);
                }
                addToNjS(s_prime, new_level, e, num_levels, n, NjS, NjS_size, NjS_index);
            }
        }
        
        // Check reset invariant after rise
        if (checkDirtInvariant(num_levels, eps, beta, D, sc_size)) {
            int critical_level = findHighestCriticalLevel(num_levels, eps, beta, D, S);
            if (critical_level >= 0) {
                resetAtLevel(critical_level, n, m, beta, num_levels, adj, element_levels, set_levels,
                            sc_size, cov_sets, cov_sets_size, cov_sets_index,
                            S, D, NjS, NjS_size, NjS_index, level_to_elements, trackSetChange);
            }
        }
    }
}


// Apply update logic (renamed from updateElement)
void Algorithm3::applyUpdate(const Update& update) {
    
    // Track sets that change during this update step for recourse calculation
    std::unordered_set<int> sets_for_recourse_check;
    
    // Store previous set_levels for sets that change
    std::unordered_map<int, int> prev_set_levels_map;
    
    // Helper function to track set level changes
    // Store lambda in std::function so it can be passed by reference to other functions
    std::function<void(int)> trackSetChange = [&](int s) {
        if (s >= 0 && s < (int)set_levels.size()) {
            sets_for_recourse_check.insert(s);
            // Store initial state if not already stored
            if (prev_set_levels_map.find(s) == prev_set_levels_map.end()) {
                prev_set_levels_map[s] = set_levels[s];
            }
        }
    };
    
    if (update.operation == 0) {
        // === INSERTION ===
        if ((int)adj.size() >= n) {
            std::cout << "ERROR: Maximum number of elements reached: " << n << std::endl;
            return; // Skip insertion
        }
        
        // Add the element with its sets from the dataset
        adj.push_back(update.sets);
        int e = adj.size() - 1; // Index of the newly inserted element
        
        // Map element_id to its index in adj
        element_id_to_index[update.element_id] = e;
        
        // Resize element_levels if needed
        if (e >= (int)element_levels.size()) {
            element_levels.resize(e + 1, -1);
        }
        element_levels[e] = -1; // Initialize as uncovered

        // Find the set with the highest level among sets containing element e
        int max_level_set = -1;
        int chosen_set = -1;
        for (int s : adj[e]) {  // check all sets containing the newly inserted element
            if (set_levels[s] > max_level_set) {
                max_level_set = set_levels[s];
                chosen_set = s;
            }
        }

        if (chosen_set != -1 && max_level_set >= 0) {
            // Found a set in the cover at level max_level_set
            // Add element e to cov_sets[chosen_set]
            addToCovSet(chosen_set, e, cov_sets, cov_sets_size, cov_sets_index);
            
            // ASSIGN ELEMENT LEVEL TO THE NEW LEVEL
            int old_level = element_levels[e]; // Will be -1 for newly inserted element
            element_levels[e] = max_level_set;
        
            // UPDATE level_to_elements
            updateLevelToElements(e, old_level, max_level_set, num_levels, level_to_elements);

            // Update NjS for all sets containing e
            for (int s_prime : adj[e]) {
                // Ensure NjS[s_prime] has enough levels
                if (max_level_set >= (int)NjS[s_prime].size()) {
                    NjS[s_prime].resize(max_level_set + 1);
                    NjS_size[s_prime].resize(max_level_set + 1, 0);
                }
                // Add element e to NjS[s_prime] at level max_level_set
                addToNjS(s_prime, max_level_set, e, num_levels, n, NjS, NjS_size, NjS_index);
            }
        } else {
            
            // No set in the cover contains this element
            // Add an arbitrary set containing e to the cover
            if (!adj[e].empty()) {
                int chosen_set = adj[e][0]; // Choose the first set containing e
                int new_level = 0; // New sets start at level 0
                // cout << "Performing naive rise" << endl;
                // Update set_levels: from -1 to 0
                trackSetChange(chosen_set); // Track this set change
                set_levels[chosen_set] = new_level;
                
                // Increment sc_size since this set was not previously in the cover
                sc_size ++;

                // Increment S[0] too
                S[0] ++;
                
                // Add element e to cov_sets[chosen_set]
                addToCovSet(chosen_set, e, cov_sets, cov_sets_size, cov_sets_index);
                
                // Set element level to 0
                element_levels[e] = new_level;
                
                // Update level_to_elements
                updateLevelToElements(e, -1, new_level, num_levels, level_to_elements);
                
                //cout << "Inserting element to level 0" << endl; 

                // Update NjS for all sets containing e
                for (int s_prime : adj[e]) {
                    // Ensure NjS[s_prime] has enough levels
                    if (new_level >= (int)NjS[s_prime].size()) {
                        NjS[s_prime].resize(new_level + 1);
                        NjS_size[s_prime].resize(new_level + 1, 0);
                    }
                    // Add element e to NjS[s_prime] at level new_level
                    addToNjS(s_prime, new_level, e, num_levels, n, NjS, NjS_size, NjS_index);
                }
            }
        }
        
        // Check if any set containing e is PD after insertion
        auto pd_result = isSetPD(e, 0, num_levels, adj, NjS_size, beta);
        if (pd_result.first != -1) {
            // Found a PD set, call rise on it
            int pd_set = pd_result.first;
            int pd_level = pd_result.second;
            rise(pd_set, pd_level, num_levels, n, m, beta, eps, set_levels, element_levels, sc_size,
                 S, D, cov_sets, cov_sets_size, adj, NjS, NjS_size, NjS_index, cov_sets_index,
                 level_to_elements, trackSetChange);
        }
    } 
    else {
        // === DELETION ===
        if (adj.empty()) return; // Nothing to delete
        
        // Find the index of the element to delete using the mapping
        auto it = element_id_to_index.find(update.element_id);
        if (it == element_id_to_index.end()) {
            // Element not found, skip this deletion
            return;
        }
        
        int del_idx = it->second;
        int deleted_elem = del_idx; // Element ID to delete (element IDs are indices)
        int moved_elem = adj.size() - 1; // Element that will be moved to del_idx

        // Get information about the deleted element before removing it
        int deleted_level = (deleted_elem < (int)element_levels.size()) ? element_levels[deleted_elem] : -1;
        
        // Remove deleted element from cov_sets if it was covered
        int covering_set = -1; // Track which set covered the deleted element

        if (deleted_level >= 0) {
            // Update relevant dirt counter: deletion of element at level j adds 1/beta^j to D[j]
            if (deleted_level < (int)D.size()) {
                D[deleted_level] += 1.0 / std::pow(beta, deleted_level);
            }
            
            // Remove from level_to_elements
            removeFromLevelToElements(deleted_elem, deleted_level, num_levels, level_to_elements);
            
            // Get which set covers this element
            auto pos = getCovSetPosition(deleted_elem, cov_sets_index);
            if (pos.first != -1) {
                covering_set = pos.first;
                // Remove element from cov_sets[covering_set]
                deleteFromCovSet(deleted_elem, cov_sets, cov_sets_size, cov_sets_index);
                // Don't call removeSetFromCoverIfEmpty here - do it after resetAtLevel to avoid double counting
            }
        }
        
        // Remove deleted element from NjS for all sets containing it
        if (deleted_level >= 0 && deleted_elem < (int)adj.size()) {
            for (int s : adj[deleted_elem]) {
                // Remove element from NjS[s][deleted_level]
                deleteFromNjS(s, deleted_level, deleted_elem, num_levels, n, NjS, NjS_size, NjS_index);
            }
        }
        
        // Now handle the element that will be moved (from position moved_elem to del_idx)
        if (moved_elem != del_idx && moved_elem < (int)adj.size()) {
            int moved_level = (moved_elem < (int)element_levels.size()) ? element_levels[moved_elem] : -1;
            
            // Update cov_sets_index: change element ID from moved_elem to del_idx
            if (moved_level >= 0) {
                auto moved_pos = getCovSetPosition(moved_elem, cov_sets_index);
                if (moved_pos.first != -1) {
                    // Update the index map: change key from moved_elem to del_idx
                    cov_sets_index[del_idx] = moved_pos;
                    cov_sets_index.erase(moved_elem);
                    // Update the actual vector entry
                    cov_sets[moved_pos.first][moved_pos.second] = del_idx;
                }
            }
            
            // Update NjS: change all references from moved_elem to del_idx
            if (moved_level >= 0) {
                for (int s : adj[moved_elem]) {
                    long long old_key = encodeKey(s, moved_level, moved_elem, num_levels, n);
                    long long new_key = encodeKey(s, moved_level, del_idx, num_levels, n);
                    
                    if (NjS_index.find(old_key) != NjS_index.end()) {
                        int idx = NjS_index[old_key];
                        NjS_index[new_key] = idx;
                        NjS_index.erase(old_key);
                        // Update the actual vector entry
                        NjS[s][moved_level][idx] = del_idx;
                    }
                }
            }
            
            // Update element_levels
            if (moved_elem < (int)element_levels.size()) {
                if (del_idx >= (int)element_levels.size()) {
                    element_levels.resize(del_idx + 1, -1);
                }
                element_levels[del_idx] = element_levels[moved_elem];
            }
            
            // Update level_to_elements: remove moved_elem from its level, add del_idx to same level
            if (moved_level >= 0) {
                removeFromLevelToElements(moved_elem, moved_level, num_levels, level_to_elements);
                addToLevelToElements(del_idx, moved_level, num_levels, level_to_elements);
            }
        }
        
        // Update the mapping for the element that will be moved to del_idx
        if (del_idx != moved_elem) {
            // Find which element_id corresponds to the last element
            for (auto& pair : element_id_to_index) {
                if (pair.second == moved_elem) {
                    pair.second = del_idx;
                    break;
                }
            }
        }
        
        // Remove the mapping for the deleted element
        element_id_to_index.erase(it);
        
        // Swap and remove from adj
        adj[del_idx] = std::move(adj.back());
        adj.pop_back();
        
        // Clear element_levels entry for deleted element if it was the last one
        if (deleted_elem == (int)adj.size() && deleted_elem < (int)element_levels.size()) {
            // The deleted element was the last one, so we can just resize
            element_levels.resize(adj.size(), -1);
        } else if (deleted_elem < (int)element_levels.size()) {

        }
        
        // Check reset invariant after deletion (no falls or ND sets)
        if (checkDirtInvariant(num_levels, eps, beta, D, sc_size)) {
            int critical_level = findHighestCriticalLevel(num_levels, eps, beta, D, S);
            if (critical_level >= 0) {
                resetAtLevel(critical_level, n, m, beta, num_levels, adj, element_levels, set_levels,
                            sc_size, cov_sets, cov_sets_size, cov_sets_index,
                            S, D, NjS, NjS_size, NjS_index, level_to_elements, trackSetChange);
            }
        }
        
        // Now check if covering_set needs to be removed (after reset, if any)
        // This avoids double counting: if reset happened, it already handled set changes
        if (covering_set >= 0) {
            removeSetFromCoverIfEmpty(covering_set, S, set_levels, cov_sets_size, sc_size, trackSetChange);
        }
    }
    
    // Calculate recourse and update set_cover vector for sets that changed during this update step
    // (runs for both insertion and deletion operations)
    int update_recourse = 0;
    for (int s : sets_for_recourse_check) {
        if (s >= 0 && s < (int)set_levels.size()) {
            auto it = prev_set_levels_map.find(s);
            if (it != prev_set_levels_map.end()) {
                int prev_level = it->second;
                int curr_level = set_levels[s];
                // Count if set changed from not in cover to in cover, or vice versa
                if ((prev_level == -1 && curr_level >= 0) || (prev_level >= 0 && curr_level == -1)) {
                    update_recourse++;
                    // Update set_cover vector based on whether set was added or removed
                    if (prev_level < 0 && curr_level >= 0) {
                        // Set was added to cover
                        addSetToCover(s);
                    } else if (prev_level >= 0 && curr_level < 0) {
                        // Set was removed from cover
                        removeSetFromCover(s);
                    }
                }
            }
        }
    }
    
    // Store values for fast getter access
    last_update_recourse = update_recourse;
    last_set_cover_size = sc_size;
}