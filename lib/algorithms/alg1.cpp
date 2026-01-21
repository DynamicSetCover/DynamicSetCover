#include "alg1.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>

// Anonymous namespace for helper functions (to avoid multiple definition errors)
namespace {
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

    // Function to perform MSC with bucket-based relaxation on participating sets only
    // Only considers sets in participating_sets, reducing complexity from O(m) to O(n*f)
    std::pair<int, std::vector<bool>> msc(int n,
                                          int m,
                                          const std::vector<std::vector<int>>& element_lst,
                                          double beta,
                                          const std::unordered_set<int>& participating_sets) {
        
        int actual_n = element_lst.size(); // actual number of active elements
        std::unordered_map<int, int> num_uncvrd; // number of uncovered elements in each participating set
        std::unordered_map<int, std::vector<int>> set_lst; // vector of elements in each participating set
        std::vector<bool> cvrd(actual_n, false); // true if element is covered
        std::unordered_map<int, bool> in_sc; // true if participating set is in the SC

        // Build set_lst and num_uncvrd for participating sets only (O(n*f))
        for (int e = 0; e < actual_n; e++) {
            for (int s : element_lst[e]) {
                if (participating_sets.find(s) != participating_sets.end()) {
                    set_lst[s].push_back(e);
                    num_uncvrd[s]++; // initially all elements are uncovered
                }
            }
        }

        // Find maximum uncovered count to determine number of buckets needed
        int max_uncvrd = 0;
        for (const auto& [s, count] : num_uncvrd) {
            if (count > max_uncvrd) {
                max_uncvrd = count;
            }
        }

        // Calculate number of buckets: log_{beta}(max_uncvrd) = log(max_uncvrd) / log(beta)
        int num_buckets = (max_uncvrd > 0) ? (int)std::ceil(std::log(max_uncvrd) / std::log(beta)) + 1 : 1;
        
        // Precompute bucket boundaries: bucket_boundaries[j] = beta^j with iterative multiplication to avoid expensive pow() calls
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
        
        std::vector<std::vector<int>> buckets(num_buckets); // buckets[j] contains set indices in bucket j

        // Initialize buckets for participating sets only
        for (const auto& [s, count] : num_uncvrd) {
            if (count > 0) {
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

        int uncovered = actual_n;

        while(uncovered > 0) {
            // Find highest non-empty bucket with a valid set (lazy update)
            int s = -1;
            while (highest_bucket >= 0) {
                // Clean up invalid sets from this bucket (sets that moved to different buckets)
                while (!buckets[highest_bucket].empty()) {
                    int candidate = buckets[highest_bucket].back();
                    buckets[highest_bucket].pop_back();
                    
                    // Check if this set is still valid (not in SC, has uncovered elements, in correct bucket)
                    auto in_sc_it = in_sc.find(candidate);
                    bool is_in_sc = (in_sc_it != in_sc.end() && in_sc_it->second);
                    auto num_uncvrd_it = num_uncvrd.find(candidate);
                    int uncvrd_count = (num_uncvrd_it != num_uncvrd.end()) ? num_uncvrd_it->second : 0;
                    
                    if (!is_in_sc && uncvrd_count > 0) {
                        int correct_bucket = getBucket(uncvrd_count, bucket_boundaries);
                        if (correct_bucket == highest_bucket) {
                            // This set is valid and in the correct bucket
                            s = candidate;
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
            
            if (s < 0) {
                break; // No more sets available
            }

            in_sc[s] = true;

            // Update uncovered counts (but don't move sets to new buckets yet - lazy update)
            auto set_lst_it = set_lst.find(s);
            if (set_lst_it != set_lst.end()) {
                for (int e : set_lst_it->second) {       
                    if (!cvrd[e]) {   
                        cvrd[e] = true;
                        uncovered--;
                        
                        // For each other participating set containing this element, decrease uncovered count
                        // Sets will be moved to correct buckets lazily when we encounter them
                        for (int s_prime : element_lst[e]) {
                            if (participating_sets.find(s_prime) != participating_sets.end()) {
                                auto it = num_uncvrd.find(s_prime);
                                if (it != num_uncvrd.end() && it->second > 0) {
                                    it->second--;
                                }
                            }
                        }
                    }
                }
            }
        }  
        
        // Build result vector<bool> for all m sets (only participating sets will be true)
        std::vector<bool> result_in_sc(m, false);
        int sc_size = 0;
        for (const auto& [s, is_covered] : in_sc) {
            if (is_covered && s >= 0 && s < m) {
                result_in_sc[s] = true;
                sc_size++;
            }
        }
        
        return {sc_size, result_in_sc};
    }
} // end anonymous namespace

// Structure to represent an update operation (internal use)
struct Update {
    int operation;  // 0 = insert, 1 = delete
    int element_id; // element identifier
    std::vector<int> sets; // sets containing the element (only for insertions)
};

// Constructor: initialize data structures
Algorithm1::Algorithm1(int n, int m, double eps)
    : n(n), m(m), eps(eps), beta(1 + eps), sc_size(0), updates_since_reset(0), batch_size(0), last_update_recourse(0), last_set_cover_size(0) {
        
    // Initialize data structures, initially set cover is empty
    in_sc.resize(m, false);
    set_cover.clear(); // Initially empty set cover
    prev_set_cover.clear();
}

// Process a single update from update sequence
void Algorithm1::updateElement(const std::vector<int>& update_row) {
    if (update_row.empty()) {
        return; // Skip empty updates
    }
    
    // Convert update_row to Update struct
    Update update;
    update.operation = update_row[0];
    update.element_id = update_row[1];
    
    // For insertions, extract sets (convert from 1-indexed to 0-indexed since datasets are 1-indexed)
    if (update.operation == 0 && update_row.size() > 2) {
        for (size_t i = 2; i < update_row.size(); ++i) {
            update.sets.push_back(update_row[i] - 1); // Convert to 0-indexed
        }
    }
    
    // Call applyUpdate with the Update struct
    applyUpdate(update);
}

// Add set s to set_cover vector
void Algorithm1::addSetToCover(int s) {
    // Check if set is already in the cover vector
    for (int set_id : set_cover) {
        if (set_id == s) {
            return; // Already in cover
        }
    }
    set_cover.push_back(s);
}

// Remove set s from set_cover vector
void Algorithm1::removeSetFromCover(int s) {
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


// Function to determine the next batch size after MSC
int Algorithm1::getBatchSize(int sc_size, double eps) {
    return std::max(1, static_cast<int>(eps * sc_size));
}

// Apply update logic (renamed from applyUpdate)
void Algorithm1::applyUpdate(const Update& update) {
    // Save prev_set_cover at the start, before any changes are made
    prev_set_cover = set_cover;
    
    // Check if reset is needed (batch size reached or first update)
    // Check BEFORE incrementing, so we know if reset will occur after this update
    bool reset_will_occur = (batch_size == 0 || updates_since_reset + 1 >= batch_size);
    
    // Track if a naive set was added (for recourse calculation when no reset)
    bool naive_set_added = false;
    
    if (update.operation == 0) {
        // === INSERTION ===
        if ((int)adj.size() >= n) {
            std::cout << "ERROR: Maximum number of elements reached: " << n << std::endl;
            // Skip insertion but still count as an update
        } else {
            // Add the element with its sets
            adj.push_back(update.sets);
            int element_index = adj.size() - 1;
            
            // Map element_id to its index in adj
            element_id_to_index[update.element_id] = element_index;
            
            // Ensure coverage only if reset won't occur (reset will compute proper cover anyway)
            if (!reset_will_occur) {
                bool covered = false;
                for (int s : adj.back()) {  // check the sets of the most recently inserted element
                    if (s >= 0 && s < m && in_sc[s]) {
                        covered = true;
                        break;
                    }
                }

                if (!covered && !adj.back().empty()) {
                    int chosen_set = adj.back()[0]; // arbitrarily choose the first set
                    if (chosen_set >= 0 && chosen_set < m) {
                        in_sc[chosen_set] = true;
                        sc_size++;
                        addSetToCover(chosen_set); // Add to set_cover vector
                        naive_set_added = true; // Track that we added a naive set
                    }
                }
            }
        }
    } 
    else {
        // === DELETION ===
        if (!adj.empty()) {
            // Find the index of the element to delete using the mapping
            auto it = element_id_to_index.find(update.element_id);
            if (it != element_id_to_index.end()) {
                int del_idx = it->second;
                int last_idx = adj.size() - 1;
                
                // Update the mapping for the element that will be moved to del_idx
                if (del_idx != last_idx) {
                    // Find which element_id corresponds to the last element
                    for (auto& pair : element_id_to_index) {
                        if (pair.second == last_idx) {
                            pair.second = del_idx;
                            break;
                        }
                    }
                }
                
                // Remove the mapping for the deleted element
                element_id_to_index.erase(it);
                
                // Swap the element-to-delete with the last one and remove it
                adj[del_idx] = std::move(adj.back());
                adj.pop_back();
            }
            // If element not found or adj is empty, just skip (but still count as update)
        }
    }
    
    // Increment updates since last reset (regardless of insertion or deletion outcome)
    updates_since_reset++;
    
    if (reset_will_occur) {
        performReset();
        
        // Calculate recourse by comparing set_cover and prev_set_cover
        // Convert both to sets for efficient comparison
        std::unordered_set<int> prev_set(prev_set_cover.begin(), prev_set_cover.end());
        std::unordered_set<int> curr_set(set_cover.begin(), set_cover.end());
        
        int update_recourse = 0;
        
        // Count sets in prev_set_cover but not in set_cover (removed)
        for (int s : prev_set_cover) {
            if (curr_set.find(s) == curr_set.end()) {
                update_recourse++;
            }
        }
        
        // Count sets in set_cover but not in prev_set_cover (added)
        for (int s : set_cover) {
            if (prev_set.find(s) == prev_set.end()) {
                update_recourse++;
            }
        }
        
        last_update_recourse = update_recourse;
    } else {
        // No reset: recourse is 1 if naive set was added, 0 otherwise
        last_update_recourse = naive_set_added ? 1 : 0;
    }
    
    // Store set cover size for fast getter access
    last_set_cover_size = sc_size;
}

// Helper function to perform MSC reset
void Algorithm1::performReset() {
    // Collect participating sets: sets that contain at least one currently active element
    // This takes O(n*f) time where f is average number of sets per element
    std::unordered_set<int> participating_sets;
    for (const auto& element_sets : adj) {
        for (int s : element_sets) {
            if (s >= 0 && s < m) {
                participating_sets.insert(s);
            }
        }
    }
    
    // Run MSC only on participating sets (O(n*f) instead of O(m))
    std::tie(sc_size, in_sc) = msc(n, m, adj, beta, participating_sets);
    
    // Update set_cover vector after reset (only check participating sets, O(n*f) instead of O(m))
    set_cover.clear();
    for (int s : participating_sets) {
        if (s >= 0 && s < m && in_sc[s]) {
            set_cover.push_back(s);
        }
    }
    
    // Set batch size for next batch (fixed until next MSC)
    batch_size = getBatchSize(sc_size, eps);
    updates_since_reset = 0;
}
