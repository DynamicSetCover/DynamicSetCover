#include "alg4.h"
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
        // Linear search from highest bucket (fast since num_buckets is small ~log(n))
        for (int j = bucket_boundaries.size() - 1; j >= 0; j--) {
            if (uncovered_count >= bucket_boundaries[j]) {
                return j;
            }
        }
        return 0; // Should not reach here, but return 0 as fallback
    }
} // end anonymous namespace

// Constructor: initialize data structures
Algorithm4::Algorithm4(int n, int m, double eps)
    : n(n), m(m), eps(eps), beta(1 + eps), sc_size(0), last_update_recourse(0), last_set_cover_size(0) {
    // Calculate number of levels
    num_levels = floor_log_beta(n, beta) + 1;
    
    // Initialize data structures
    element_levels.resize(n, -1);
    passive_levels.resize(n, num_levels);  // Initialize to max level + 1
    set_levels.resize(m, -1);
    cov_sets.resize(m);
    cov_sets_size.resize(m, 0);
    level_to_elements.resize(num_levels + 1);
    A.resize(num_levels, 0);
    P.resize(num_levels, 0);
    D.resize(num_levels, 0);
    set_cover.clear(); // Initially empty set cover
}

// Process a single update from update sequence
void Algorithm4::updateElement(const std::vector<int>& update_row) {
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
void Algorithm4::addSetToCover(int s) {
    // Check if set is already in the cover vector
    for (int set_id : set_cover) {
        if (set_id == s) {
            return; // Already in cover
        }
    }
    set_cover.push_back(s);
}

// Remove set s from set_cover vector (maintains O(1) access)
void Algorithm4::removeSetFromCover(int s) {
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
                                std::vector<int>& set_levels,
                                std::vector<int>& cov_sets_size,
                                int& sc_size,
                                std::function<void(int)>& trackSetChange) {
        if (set_levels[s] >= 0 && cov_sets_size[s] == 0) {
            trackSetChange(s); // Track this set change
            set_levels[s] = -1;
            sc_size--;
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

    // Update counters when inserting element e at level lev with plev(e) = lev
    void updateCountersOnInsert(int lev,
                            int num_levels,
                            std::vector<int>& A,
                            std::vector<int>& P) {
        // Since lev = plev, element contributes to P(k) for k >= lev
        // But does NOT contribute to A(k) since lev(e) = plev(e) = lev (not lev(e) <= k and plev(e) > k)
        for (int k = lev; k < num_levels; ++k) {
            P[k]++;
        }
    }

    // Update counters when deleting element e at level i with passive level j
    void updateCountersOnDelete(int i,
                            int j,
                            int num_levels,
                            std::vector<int>& A,
                            std::vector<int>& P,
                            std::vector<int>& D) {
        // Lower A(k) for i <= k < j
        for (int k = i; k < j && k < num_levels; ++k) {
            A[k]--;  // Remove safety check to detect bugs
            if (A[k] < 0) {
                std::cerr << "WARNING: A[" << k << "] became negative after deletion of element at level " << i << " with plev " << j << std::endl;
                A[k] = 0;  // Clamp to 0 but warn
            }
        }
        
        // Lower P(k) for k >= j
        for (int k = j; k < num_levels; ++k) {
            P[k]--;  // Remove safety check to detect bugs
            if (P[k] < 0) {
                std::cerr << "WARNING: P[" << k << "] became negative after deletion of element at level " << i << " with plev " << j << std::endl;
                P[k] = 0;  // Clamp to 0 but warn
            }
        }
        
        // Note: We only increment D[i], not all D[k] for k >= i, because the invariant
        // uses sum_{i=0}^k D(i), so incrementing at each level would double-count
        D[i]++;
    }

    // Check invariant: P(k) + sum_{i=0}^k D(i) <= 2*eps*A(k)
    // Returns the largest k such that the invariant is violated, or -1 if no violation
    int checkInvariant(int num_levels,
                   double eps,
                   const std::vector<int>& A,
                   const std::vector<int>& P,
                   const std::vector<int>& D) {
        int violating_k = -1;
        
        // Compute sum of all D(i) for highest level first (O(num_levels))
        int sum_D = 0;
        for (int i = 0; i < num_levels; ++i) {
            sum_D += D[i];
        }
        
        // Check from high to low levels, subtracting D[k+1] from sum_D each iteration
        for (int k = num_levels - 1; k >= 0; --k) {
            // Check if P(k) + sum_D > 2*eps*A(k)
            // If A(k) = 0, the invariant requires P(k) + sum_D <= 0, so both must be 0
            if (A[k] == 0) {
                if (P[k] + sum_D > 0) {
                    violating_k = k;
                    break;
                }
            } else {
                // Normal case: check if P(k) + sum_D > 2*eps*A(k)
                if (P[k] + sum_D > 2 * eps * A[k]) {
                    violating_k = k;
                    break;
                }
            }
            
            // Subtract D[k] from sum_D for next iteration (k-1)
            // Note: for k=0, we subtract D[0] but then exit, so it's fine
            if (k > 0) {
                sum_D -= D[k];
            }
        }
        
        return violating_k;
    }

    // Reset procedure: run msc on elements at level <=k and sets containing them
    // Truncates levels to k+1, sets plev(e) = max(previous_plev, k+1)
    // Resets counters for k' <= k
    // Uses level_indexed data structures to optimize to O(n_k) instead of O(n)
    void resetAtLevel(int k,
                  int n,
                  int m,
                  double beta,
                  int num_levels,
                  std::vector<std::vector<int>>& adj,
                  std::vector<int>& element_levels,
                  std::vector<int>& set_levels,
                  int& sc_size,
                  std::vector<std::vector<int>>& cov_sets,
                  std::vector<int>& cov_sets_size,
                  std::unordered_map<int, std::pair<int, int>>& cov_sets_index,
                  std::vector<int>& passive_levels,
                  std::vector<int>& A,
                  std::vector<int>& P,
                  std::vector<int>& D,
                  std::vector<std::vector<int>>& level_to_elements,
                  std::function<void(int)>& trackSetChange,
                  const std::vector<int>& prev_set_levels_for_reset) {
                    
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
        
        // Store previous passive levels for elements being reset using existing participating_elements set - O(n_k) time.
        std::vector<int> prev_plev;
        std::vector<int> participating_elements_vec;
        for (int e : participating_elements) {
            if (e < (int)adj.size()) {
                prev_plev.push_back((e < (int)passive_levels.size()) ? passive_levels[e] : num_levels + 1);
                participating_elements_vec.push_back(e);
            }
        }
        
        // Store previous set_levels for participating sets to calculate recourse in O(|participating_sets|)
        // Use prev_set_levels_for_reset (state before updateElement started) instead of current state
        // Add all participating sets to tracking, regardless of whether their level changes
        std::unordered_map<int, int> prev_set_levels_map;
        for (int s : participating_sets) {
            if (s < (int)prev_set_levels_for_reset.size()) {
                prev_set_levels_map[s] = prev_set_levels_for_reset[s];
            } else if (s < (int)set_levels.size()) {
                prev_set_levels_map[s] = -1; // Set didn't exist before, treat as -1
            }
            // Add all participating sets to tracking
            trackSetChange(s);
        }
        
        // Count how many participating sets were in the cover before msc (for efficient sc_size maintenance)
        int participating_sets_in_cover_before = 0;
        for (int s : participating_sets) {
            if (s < (int)set_levels.size() && set_levels[s] >= 0) {
                participating_sets_in_cover_before++;
            }
        }
        
        // Subtract participating sets from sc_size before msc
        sc_size -= participating_sets_in_cover_before;
        
        // Run msc with truncation to k+1, passing level_indexed structure for maintenance
        msc(n, m, beta, num_levels, adj, element_levels, set_levels, sc_size,
            cov_sets, cov_sets_size, cov_sets_index,
            k + 1, participating_elements, participating_sets,
            &level_to_elements);
        
        // Count how many participating sets are in the cover after msc and add to sc_size
        int participating_sets_in_cover_after = 0;
        for (int s : participating_sets) {
            if (s < (int)set_levels.size() && set_levels[s] >= 0) {
                participating_sets_in_cover_after++;
            }
        }
        sc_size += participating_sets_in_cover_after;
        
        // Update passive levels: plev(e) = max(previous_plev, k+1) using participating_elements_vec list - O(n_k) time.
        for (size_t i = 0; i < participating_elements_vec.size(); ++i) {
            int e = participating_elements_vec[i];
            if (e >= (int)passive_levels.size()) {
                passive_levels.resize(e + 1, num_levels + 1);
            }
            int old_plev = (i < prev_plev.size()) ? prev_plev[i] : num_levels + 1;
            passive_levels[e] = std::max(old_plev, k + 1);
        }
        
        // Reset counters for k' <= k according to definitions
        for (int k_prime = 0; k_prime <= k && k_prime < num_levels; ++k_prime) {
            P[k_prime] = 0; // P(k') = 0 for k' <= k (since all participating elements have plev >= k+1 > k')
            D[k_prime] = 0; // D(k') = 0 following the reset
            
            // A(k') = number of elements e such that lev(e) <= k' and plev(e) > k', only iterate through elements at level <= k' using level_to_elements for O(n_k) time.
            A[k_prime] = 0;
            for (int lev = 0; lev <= k_prime && lev < num_levels; ++lev) {
                for (int e : level_to_elements[lev]) {
                    if (e < (int)adj.size() && e < (int)passive_levels.size()) {
                        int plev = passive_levels[e];
                        if (plev > k_prime) {
                            A[k_prime]++;
                        }
                    }
                }
            }
        }
    }
} // end anonymous namespace


// Apply update logic (renamed from updateElement)
void Algorithm4::applyUpdate(const Update& update) {
    
    // Store initial set_levels state before any changes (for resetAtLevel to use correct previous state)
    std::vector<int> initial_set_levels = set_levels;
    
    // Track sets that change during this update step for recourse calculation
    std::unordered_set<int> sets_for_recourse_check;
    
    // Store previous set_levels for sets that change
    // Use initial_set_levels as the baseline (state before updateElement started)
    std::unordered_map<int, int> prev_set_levels_map;
    
    // Helper function to track set level changes
    // Store lambda in std::function so it can be passed by reference to other functions
    std::function<void(int)> trackSetChange = [&](int s) {
        if (s >= 0 && s < (int)set_levels.size()) {
            sets_for_recourse_check.insert(s);
            // Store initial state if not already stored (use initial_set_levels as baseline)
            if (prev_set_levels_map.find(s) == prev_set_levels_map.end()) {
                if (s < (int)initial_set_levels.size()) {
                    prev_set_levels_map[s] = initial_set_levels[s];
                } else {
                    prev_set_levels_map[s] = -1; // Set didn't exist before, treat as -1
                }
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
        int e = adj.size() - 1; // index of the newly inserted element
        
        // Map element_id to its index in adj
        element_id_to_index[update.element_id] = e;
        
        // RESIZE ELEMENT_LEVELS IF NEEDED
        if (e >= (int)element_levels.size()) {
            element_levels.resize(e + 1, -1);
        }
        element_levels[e] = -1; // Initialize as uncovered

        // FIND THE SET WITH THE HIGHEST LEVEL AMONG SETS CONTAINING ELEMENT E
        int max_level_set = -1;
        int chosen_set = -1;
        for (int s : adj[e]) {  // check all sets containing the newly inserted element
            if (set_levels[s] > max_level_set) {
                max_level_set = set_levels[s];
                chosen_set = s;
            }
        }

        int new_level; // new level for the element
        int chosen_set_final; // the set that will contain the element at the new level
        bool naive_set_added = false; // Track if we added a naive set
        
        if (chosen_set != -1 && max_level_set >= 0) {
            // Found a set in the cover at level max_level_set
            chosen_set_final = chosen_set;
            new_level = max_level_set;
        } else {
            // No set in the cover contains this element, so we choose an arbitrary set and take it to level 0
            if (!adj[e].empty()) {
                chosen_set_final = adj[e][0]; // Choose the first set containing e
                new_level = 0; // New sets start at level 0
                trackSetChange(chosen_set_final); // Track this set change
                set_levels[chosen_set_final] = new_level; // Update the set level
                sc_size++; // This set is now being added to the cover
                naive_set_added = true; // Mark that we added a naive set
            } else {
                // No sets containing e, this is an error
                std::cerr << "ERROR: Element " << e << " has no sets containing it!" << std::endl;
                return;
            }
        }
        
        // ADD ELEMENT e TO cov_sets[chosen_set_final]
        addToCovSet(chosen_set_final, e, cov_sets, cov_sets_size, cov_sets_index);
        
        // ASSIGN ELEMENT LEVEL TO THE NEW LEVEL
        int old_level = element_levels[e]; // Will be -1 for newly inserted element
        element_levels[e] = new_level;
        
        // UPDATE level_to_elements
        updateLevelToElements(e, old_level, new_level, num_levels, level_to_elements);
        
        // SET PASSIVE LEVEL TO THE NEW LEVEL (SINCE lev = plev)
        if (e >= (int)passive_levels.size()) {
            passive_levels.resize(e + 1, num_levels + 1);
        }
        passive_levels[e] = new_level;
        
        // UPDATE COUNTERS: ELEMENT AT LEVEL new_level WITH plev = new_level
        updateCountersOnInsert(new_level, num_levels, A, P);
        
        // CHECK INVARIANT AND RESET IF NEEDED (check for highest violating k)
        int violating_k = checkInvariant(num_levels, eps, A, P, D);
        bool reset_occurred = false;
        if (violating_k >= 0) {
            resetAtLevel(violating_k, n, m, beta, num_levels, adj, element_levels, set_levels,
                        sc_size, cov_sets, cov_sets_size, cov_sets_index,
                        passive_levels, A, P, D, level_to_elements, trackSetChange, initial_set_levels);
            reset_occurred = true;
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
        int moved_elem = adj.size() - 1; // Element that will be moved to del_idx to avoid gaps in the adjacency list

        // GET INFORMATION ABOUT THE DELETED ELEMENT BEFORE REMOVING IT
        int deleted_level = element_levels[deleted_elem];
        int deleted_plev = passive_levels[deleted_elem];
        
        // UPDATE COUNTERS BEFORE REMOVING THE ELEMENT
        updateCountersOnDelete(deleted_level, deleted_plev, num_levels, A, P, D);
        
        // REMOVE DELETED ELEMENT FROM level_to_elements
        removeFromLevelToElements(deleted_elem, deleted_level, num_levels, level_to_elements);
        
        // REMOVE DELETED ELEMENT FROM cov_sets
        int covering_set = -1; // Track which set covered the deleted element
        auto pos = getCovSetPosition(deleted_elem, cov_sets_index); // Get which set covers this element
        if (pos.first != -1) {
            covering_set = pos.first;
            deleteFromCovSet(deleted_elem, cov_sets, cov_sets_size, cov_sets_index); // Remove element from cov_sets[covering_set]
            // Don't call removeSetFromCoverIfEmpty here - do it after resetAtLevel to avoid double counting
        }
        
        // Now handle the element that will be moved (from position moved_elem to del_idx)
        if (moved_elem != del_idx && moved_elem < (int)adj.size()) {
            int moved_level = element_levels[moved_elem];
            
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
            
            // Update element_levels
            if (del_idx >= (int)element_levels.size()) {
                element_levels.resize(del_idx + 1, -1);
            }
            int moved_new_level = element_levels[moved_elem];
            element_levels[del_idx] = moved_new_level;
            
            // Update level_to_elements: remove moved_elem from its level, add del_idx to same level
            removeFromLevelToElements(moved_elem, moved_level, num_levels, level_to_elements);
            addToLevelToElements(del_idx, moved_new_level, num_levels, level_to_elements);
            
            // Update passive_levels
            if (del_idx >= (int)passive_levels.size()) {
                passive_levels.resize(del_idx + 1, num_levels + 1);
            }
            passive_levels[del_idx] = passive_levels[moved_elem];
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
        if (del_idx < (int)adj.size()) {
            adj[del_idx] = std::move(adj.back());
        }
        adj.pop_back();
        
        // Clear element_levels entry for deleted element if it was the last one
        if (deleted_elem == (int)adj.size() && deleted_elem < (int)element_levels.size()) {
            // The deleted element was the last one, so we can just resize
            element_levels.resize(adj.size(), -1);
        } else if (deleted_elem < (int)element_levels.size()) {
            // Mark as deleted (though this entry might now be used by moved element)
            // Actually, if moved_elem != del_idx, we already updated it above
        }
        
        // Check invariant and reset if needed (check for highest violating k, after all deletion bookkeeping is complete)
        int violating_k = checkInvariant(num_levels, eps, A, P, D);
        if (violating_k >= 0) {
            resetAtLevel(violating_k, n, m, beta, num_levels, adj, element_levels, set_levels,
                        sc_size, cov_sets, cov_sets_size, cov_sets_index,
                        passive_levels, A, P, D, level_to_elements, trackSetChange, initial_set_levels);
        }
        
        // Now check if covering_set needs to be removed (after reset, if any)
        // This avoids double counting: if reset happened, it already handled set changes
        if (covering_set >= 0) {
            removeSetFromCoverIfEmpty(covering_set, set_levels, cov_sets_size, sc_size, trackSetChange);
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