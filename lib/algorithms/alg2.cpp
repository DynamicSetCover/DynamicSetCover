#include "alg2.h"
#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <random>
#include <chrono>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <unordered_map>
#include <utility>
#include <sstream>
#include <string>
#include <functional>
#include <unordered_set>

// Structure to represent an update operation (internal use)
struct Update {
    int operation;  // 0 = insert, 1 = delete
    int element_id; // element identifier
    std::vector<int> sets; // sets containing the element (only for insertions)
};

// Helper: floor log base beta for n >= 1
static inline int floor_log_beta(int n, double beta) {
    if (n <= 1) return 0;
    return static_cast<int>(std::log(n) / std::log(beta)); // safe for n>=1
}

// Constructor: initialize data structures
Algorithm2::Algorithm2(int n, int m, double eps) 
    : n(n), m(m), beta(1 + eps), sc_size(0), last_update_recourse(0), last_set_cover_size(0) {
    // Calculate number of levels
    num_levels = floor_log_beta(n, beta) + 1;
    
    // Initialize data structures
    element_levels.resize(n, -1);
    set_levels.resize(m, -1);
    cov_sets.resize(m);
    cov_sets_size.resize(m, 0);
    NjS.resize(m, std::vector<std::vector<int>>(num_levels));
    NjS_size.resize(m, std::vector<int>(num_levels, 0));
    set_cover.clear(); // Initially empty set cover
}

// Add set s to set_cover vector (maintains O(1) access)
void Algorithm2::addSetToCover(int s) {
    // Check if set is already in the cover vector
    for (int set_id : set_cover) {
        if (set_id == s) {
            return; // Already in cover
        }
    }
    set_cover.push_back(s);
}

// Remove set s from set_cover vector (maintains O(1) access)
void Algorithm2::removeSetFromCover(int s) {
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

// Process a single update from update sequence
void Algorithm2::updateElement(const std::vector<int>& update_row) {
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

// Helper: Encode (s, j, e) into a single integer key for the hash map
// Assumes: s < m, j < num_levels, e < n
static inline long long encodeKey(int s, int j, int e, int num_levels, int n) {
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
    long long key = encodeKey(s, j, e, num_levels, n);
    
    // Check if element exists in the index map
    if (NjS_index.find(key) == NjS_index.end()) {
        return false; // Element not found
    }
    
    int idx = NjS_index[key]; // Get the index of element e
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
    NjS_size[s][j]--;
    
    return true;
}

// Helper: Add element e to NjS[s][j] and maintain index map
void addToNjS(int s, int j, int e,
               int num_levels, int n,
               std::vector<std::vector<std::vector<int>>>& NjS,
               std::vector<std::vector<int>>& NjS_size,
               std::unordered_map<long long, int>& NjS_index) {
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
// Note: This is a standalone function, but we need to call Algorithm2::removeSetFromCover
// We'll handle this by checking in applyUpdate after calling this function
void removeSetFromCoverIfEmpty(int s,
                                std::vector<int>& set_levels,
                                std::vector<int>& cov_sets_size,
                                int& sc_size,
                                std::function<void(int)>& trackSetChange) {
    if (set_levels[s] >= 0 && cov_sets_size[s] == 0) {
        trackSetChange(s); // Track this set change
        set_levels[s] = -1;
        sc_size--;
        //cout << "Removed set from cover." << endl;
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

// Returns number of elements in s at level less than j
int numElementsLessThan(int s, int j, const std::vector<std::vector<int>>& NjS_size) {
    int sum = 0;
    for (int i = 0 ; i < j ; i ++) {
        if (i < (int)NjS_size[s].size()) {
            sum += NjS_size[s][i];
        }
    }
    return (sum);
}

// Check if a set s is ND (Negative Dirty)
// A set s is ND if it is in the cover and cov_sets[s].size() < beta^(set_levels[s] - 1)
// Returns: true if set s is ND, false otherwise (or if s is not in the cover)
bool isSetND(int s,
              const std::vector<int>& set_levels,
              const std::vector<int>& cov_sets_size,
              double beta) {
    // Check if set s is in the cover
    if (set_levels[s] < 0) {
        return false; // Set is not in the cover
    }
    
    int level = set_levels[s];
    int cov_size = cov_sets_size[s];
    
    // Calculate beta^(level - 1)
    double threshold = std::pow(beta, level - 1);
    
    // Check if set is ND: cov_size < beta^(level - 1)
    return (cov_size < threshold);
}

// Forward declarations for rise and fall
std::vector<int> rise(int s, int j,
          int num_levels, int n, double beta,
          std::vector<int>& set_levels,
          std::vector<int>& element_levels,
          int& sc_size,
          std::vector<std::vector<int>>& cov_sets,
          std::vector<int>& cov_sets_size,
          const std::vector<std::vector<int>>& adj,
          std::vector<std::vector<std::vector<int>>>& NjS,
          std::vector<std::vector<int>>& NjS_size,
          std::unordered_map<long long, int>& NjS_index,
          std::unordered_map<int, std::pair<int, int>>& cov_sets_index,
          std::function<void(int)>& trackSetChange);

std::vector<int> fall(int s,
          int num_levels, int n, double beta,
          std::vector<int>& set_levels,
          std::vector<int>& element_levels,
          int& sc_size,
          std::vector<std::vector<int>>& cov_sets,
          std::vector<int>& cov_sets_size,
          const std::vector<std::vector<int>>& adj,
          std::vector<std::vector<std::vector<int>>>& NjS,
          std::vector<std::vector<int>>& NjS_size,
          std::unordered_map<long long, int>& NjS_index,
          std::unordered_map<int, std::pair<int, int>>& cov_sets_index,
          std::function<void(int)>& trackSetChange);

// Forward declaration for fallingPhase
void fallingPhase(std::vector<int>& participating_sets,
                  int num_levels, int n, double beta,
                  std::vector<int>& set_levels,
                  std::vector<int>& element_levels,
                  int& sc_size,
                  std::vector<std::vector<int>>& cov_sets,
                  std::vector<int>& cov_sets_size,
                  const std::vector<std::vector<int>>& adj,
                  std::vector<std::vector<std::vector<int>>>& NjS,
                  std::vector<std::vector<int>>& NjS_size,
                  std::unordered_map<long long, int>& NjS_index,
                  std::unordered_map<int, std::pair<int, int>>& cov_sets_index,
                  std::function<void(int)>& trackSetChange);

void risingPhase(std::vector<int>& participating_sets,
                 int num_levels, int n, double beta,
                 std::vector<int>& set_levels,
                 std::vector<int>& element_levels,
                 int& sc_size,
                 std::vector<std::vector<int>>& cov_sets,
                 std::vector<int>& cov_sets_size,
                 const std::vector<std::vector<int>>& adj,
                 std::vector<std::vector<std::vector<int>>>& NjS,
                 std::vector<std::vector<int>>& NjS_size,
                 std::unordered_map<long long, int>& NjS_index,
                 std::unordered_map<int, std::pair<int, int>>& cov_sets_index,
                 std::function<void(int)>& trackSetChange) {
                    
    std::vector<int> rising_elements; // Collect all elements whose level changed due to rises

    // Check levels from high to low
    for (int j = num_levels; j > 0; --j) {
        // For each set s in the participating sets
        for (size_t idx = 0 ; idx < participating_sets.size() ; ++idx) {
            int s = participating_sets[idx];
            // Check if set s is PD with respect to level j
            if (numElementsLessThan(s, j, NjS_size) >= std::pow(beta, j + 1)) {
                // Call rise and collect the elements that participated
                std::vector<int> elements_from_rise = rise(s, j, num_levels, n, beta, set_levels, element_levels, sc_size,
                 cov_sets, cov_sets_size, adj, NjS, NjS_size, NjS_index, cov_sets_index, trackSetChange);
                // Add to rising_elements
                rising_elements.insert(rising_elements.end(), elements_from_rise.begin(), elements_from_rise.end());
            }
        }
    }
    
    // Create list of all sets that in their cov set contain an element from rising_elements
    std::unordered_map<int, bool> sets_to_fall_map; // Use map to avoid duplicates
    for (int e : rising_elements) {
        auto pos = getCovSetPosition(e, cov_sets_index);
        if (pos.first != -1) {
            sets_to_fall_map[pos.first] = true;
        }
    }
    
    // Convert to vector
    std::vector<int> sets_for_falling_phase;
    for (const auto& pair : sets_to_fall_map) {
        sets_for_falling_phase.push_back(pair.first);
    }
    
    // Call fallingPhase with this list of sets if non-empty
    if (!sets_for_falling_phase.empty()) {
        fallingPhase(sets_for_falling_phase, num_levels, n, beta, set_levels, element_levels, sc_size,
                     cov_sets, cov_sets_size, adj, NjS, NjS_size, NjS_index, cov_sets_index, trackSetChange);
    }
}

void fallingPhase(std::vector<int>& participating_sets,
                  int num_levels, int n, double beta,
                  std::vector<int>& set_levels,
                  std::vector<int>& element_levels,
                  int& sc_size,
                  std::vector<std::vector<int>>& cov_sets,
                  std::vector<int>& cov_sets_size,
                  const std::vector<std::vector<int>>& adj,
                  std::vector<std::vector<std::vector<int>>>& NjS,
                  std::vector<std::vector<int>>& NjS_size,
                  std::unordered_map<long long, int>& NjS_index,
                  std::unordered_map<int, std::pair<int, int>>& cov_sets_index,
                  std::function<void(int)>& trackSetChange) {
                    
    std::vector<int> falling_elements; // Collect all elements whose level changed due to falls
    
    for (size_t idx = 0 ; idx < participating_sets.size() ; ++idx) {
        int s = participating_sets[idx];
        // Check if set s is ND
        if (isSetND(s, set_levels, cov_sets_size, beta)) {
            // Call fall and collect the elements that participated
            std::vector<int> elements_from_fall = fall(s, num_levels, n, beta, set_levels, element_levels, sc_size, 
             cov_sets, cov_sets_size, adj, NjS, NjS_size, NjS_index, cov_sets_index, trackSetChange);
            // Add to falling_elements
            falling_elements.insert(falling_elements.end(), elements_from_fall.begin(), elements_from_fall.end());
        }
    }
    
    // Create list of all sets that contain an element from falling_elements
    // (not just cov sets, all sets containing these elements)
    std::unordered_map<int, bool> sets_to_rise_map; // Use map to avoid duplicates
    for (int e : falling_elements) {
        if (e < (int)adj.size()) {
            for (int s : adj[e]) {
                sets_to_rise_map[s] = true;
            }
        }
    }
    
    // Convert to vector
    std::vector<int> sets_for_rising_phase;
    for (const auto& pair : sets_to_rise_map) {
        sets_for_rising_phase.push_back(pair.first);
    }
    
    // Call risingPhase with this list of sets if non-empty
    if (!sets_for_rising_phase.empty()) {
        risingPhase(sets_for_rising_phase, num_levels, n, beta, set_levels, element_levels, sc_size,
                   cov_sets, cov_sets_size, adj, NjS, NjS_size, NjS_index, cov_sets_index, trackSetChange);
    }
}


// Apply update logic (renamed from updateElement)
void Algorithm2::applyUpdate(const Update& update) {
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
            
            // Assign element level to the level of the covering set
            element_levels[e] = max_level_set;
            
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
    
                // Update set_levels: from -1 to 0
                trackSetChange(chosen_set); // Track this set change
                set_levels[chosen_set] = new_level;
                
                // Increment sc_size since this set was not previously in the cover
                sc_size++;
                
                // Add element e to cov_sets[chosen_set]
                addToCovSet(chosen_set, e, cov_sets, cov_sets_size, cov_sets_index);
                
                // Set element level to 0
                element_levels[e] = new_level;
                
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
        // Call the rising phase with all sets containing e
        std::vector<int> initial_sets;
        for (int s : adj[e]) {
            initial_sets.push_back(s);
        }
        if (!initial_sets.empty()) {
            risingPhase(initial_sets, num_levels, n, beta, set_levels, element_levels, sc_size,
                       cov_sets, cov_sets_size, adj, NjS, NjS_size, NjS_index, cov_sets_index, trackSetChange);
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
        bool needs_fall = false; // Track if covering_set needs to fall after all updates
        if (deleted_level >= 0) {
            // Get which set covers this element
            auto pos = getCovSetPosition(deleted_elem, cov_sets_index);
            if (pos.first != -1) {
                covering_set = pos.first;
                // Remove element from cov_sets[covering_set]
                deleteFromCovSet(deleted_elem, cov_sets, cov_sets_size, cov_sets_index);
                
                // Check if the set's cover set is now empty
                removeSetFromCoverIfEmpty(covering_set, set_levels, cov_sets_size, sc_size, trackSetChange);
                
                // Check if the set is ND (will check again at end after all updates)
                if (set_levels[covering_set] >= 0 && isSetND(covering_set, set_levels, cov_sets_size, beta)) {
                    needs_fall = true;
                }
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
        
        // Check if covering_set needs to fall after all data structures are updated
        if (needs_fall && covering_set >= 0 && set_levels[covering_set] >= 0 && isSetND(covering_set, set_levels, cov_sets_size, beta)) {
            // The set covering_set is now ND (Negative Dirty)
            // Call falling phase with this set
            std::vector<int> initial_sets;
            initial_sets.push_back(covering_set);
            fallingPhase(initial_sets, num_levels, n, beta, set_levels, element_levels, sc_size,
                        cov_sets, cov_sets_size, adj, NjS, NjS_size, NjS_index, cov_sets_index, trackSetChange);
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
                
                // Check if set changed from not in cover to in cover, or vice versa
                if ((prev_level == -1 && curr_level >= 0) || (prev_level >= 0 && curr_level == -1)) {
                    update_recourse++;
                }
                
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
    
    // Store values for fast getter access
    last_update_recourse = update_recourse;
    last_set_cover_size = sc_size;
}


// Fall function: Lower the level of set s to appropriate level j
// Finds level j such that beta^j <= cov_sets_size[s] < beta^{j+1}
// Updates set_levels[s] and element_levels[e] for all e in cov_sets[s]
// Updates NjS structures accordingly
// Returns: vector of elements whose level changed due to this fall
std::vector<int> fall(int s,
          int num_levels, int n, double beta,
          std::vector<int>& set_levels,
          std::vector<int>& element_levels,
          int& sc_size,
          std::vector<std::vector<int>>& cov_sets,
          std::vector<int>& cov_sets_size,
          const std::vector<std::vector<int>>& adj,
          std::vector<std::vector<std::vector<int>>>& NjS,
          std::vector<std::vector<int>>& NjS_size,
          std::unordered_map<long long, int>& NjS_index,
          std::unordered_map<int, std::pair<int, int>>& cov_sets_index,
          std::function<void(int)>& trackSetChange) {
    std::vector<int> affected_elements; // Elements whose level changed
    
    // Get current level of set s
    int old_level = set_levels[s];
    if (old_level < 0) {
        return affected_elements; // Set is not in the cover, nothing to do
    }
    
    // Calculate new level j such that beta^j <= cov_sets_size[s] < beta^{j+1}
    int cov_size = cov_sets_size[s];
    int new_level = floor_log_beta(cov_size, beta);
    // cout << "Performing fall from level " << old_level << " to level " << new_level << endl;
    
    // Clamp to valid range
    if (new_level >= num_levels) {
        new_level = num_levels - 1;
    }
    if (new_level < 0) {
        new_level = 0; // Ensure non-negative
    }
    
    // If level hasn't changed, nothing to do
    if (new_level == old_level) {
        return affected_elements;
    }
    
    // Update set level
    trackSetChange(s); // Track this set change
    set_levels[s] = new_level;
    
    // Update all elements in cov_sets[s]
    for (int e : cov_sets[s]) {
        // Get the element's previous level (should equal old_level since all elements in cov_sets[s] have the set's level)
        int element_old_level = element_levels[e];
        
        // Remove element e from NjS[s_prime][element_old_level] for all sets s_prime containing e
        for (int s_prime : adj[e]) {
            deleteFromNjS(s_prime, element_old_level, e, num_levels, n, NjS, NjS_size, NjS_index);
        }
        
        // Add element e to NjS[s_prime][new_level] for all sets s_prime containing e
        for (int s_prime : adj[e]) {
            // Ensure NjS[s_prime] has enough levels
            if (new_level >= (int)NjS[s_prime].size()) {
                NjS[s_prime].resize(new_level + 1);
                NjS_size[s_prime].resize(new_level + 1, 0);
            }
            addToNjS(s_prime, new_level, e, num_levels, n, NjS, NjS_size, NjS_index);
        }
        
        // Update element level
        if (e < (int)element_levels.size()) {
            element_levels[e] = new_level;
        }
        
        // Add to affected elements (all elements in cov_sets[s] changed level)
        affected_elements.push_back(e);
    }
    
    // Check if set s's cover set is now empty after fall
    removeSetFromCoverIfEmpty(s, set_levels, cov_sets_size, sc_size, trackSetChange);
    
    return affected_elements;
}


// Rise function: Raise the level of set s to j+1
// Takes all elements in set s at level less than j and creates a new cov set at level j+1
// Updates set_levels[s] and element_levels[e] for all such elements
// Updates NjS and cov_sets structures accordingly
// Returns: vector of elements whose level changed due to this rise
std::vector<int> rise(int s, int j,
          int num_levels, int n, double beta,
          std::vector<int>& set_levels,
          std::vector<int>& element_levels,
          int& sc_size,
          std::vector<std::vector<int>>& cov_sets,
          std::vector<int>& cov_sets_size,
          const std::vector<std::vector<int>>& adj,
          std::vector<std::vector<std::vector<int>>>& NjS,
          std::vector<std::vector<int>>& NjS_size,
          std::unordered_map<long long, int>& NjS_index,
          std::unordered_map<int, std::pair<int, int>>& cov_sets_index,
          std::function<void(int)>& trackSetChange) {
    std::vector<int> affected_elements; // Elements whose level changed
    
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

    // Get old set level for S counter update
    int old_set_level = set_levels[s];
    
    // Update set level
    trackSetChange(s); // Track this set change
    set_levels[s] = new_level;

    // Update sc_size
    if (old_set_level == -1) {
        sc_size ++;
    }
  
    // Process each element to rise
    for (int e : elements_to_rise) {
        // Get the element's old level (should be >= 0 since e is in NjS[s][i])
        int old_level = element_levels[e];
        auto old_pos = getCovSetPosition(e, cov_sets_index);
        if (old_pos.first != -1 && old_pos.first != s) {
            // Element was covered by a different set, remove it
            deleteFromCovSet(e, cov_sets, cov_sets_size, cov_sets_index);
            
            // Check if the old covering set's cover set is now empty
            removeSetFromCoverIfEmpty(old_pos.first, set_levels, cov_sets_size, sc_size, trackSetChange);
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
        
        // Add element to NjS[s_prime][new_level] for all sets s_prime containing e
        for (int s_prime : adj[e]) {
            // Ensure NjS[s_prime] has enough levels
            if (new_level >= (int)NjS[s_prime].size()) {
                NjS[s_prime].resize(new_level + 1);
                NjS_size[s_prime].resize(new_level + 1, 0);
            }
            addToNjS(s_prime, new_level, e, num_levels, n, NjS, NjS_size, NjS_index);
        }
        
        // Add to affected elements (all elements that participated in the rise)
        affected_elements.push_back(e);
    }
    
    // Check if set s's cover set is now empty after rise
    // (This could happen if all elements were removed and none were added)
    removeSetFromCoverIfEmpty(s, set_levels, cov_sets_size, sc_size, trackSetChange);
    
    return affected_elements;
}