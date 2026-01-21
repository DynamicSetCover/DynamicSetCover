#ifndef ALG2_H
#define ALG2_H

#include <vector>
#include <unordered_map>
#include <functional>

// Forward declarations
struct Update;

class Algorithm2 {
private:
    // Parameters
    int n;
    int m;
    int num_levels;
    double beta;
    
    // Data structures
    std::vector<std::vector<int>> adj;
    std::vector<int> element_levels;
    std::vector<int> set_levels;
    int sc_size;
    std::vector<std::vector<int>> cov_sets;
    std::vector<int> cov_sets_size;
    std::vector<std::vector<std::vector<int>>> NjS;
    std::vector<std::vector<int>> NjS_size;
    std::unordered_map<long long, int> NjS_index;
    std::unordered_map<int, std::pair<int, int>> cov_sets_index;
    std::unordered_map<int, int> element_id_to_index;
    
    // Per-update tracking (stored after each update)
    int last_update_recourse;
    int last_set_cover_size;
    
    // Set cover vector (maintained for O(1) access)
    // Stores 0-indexed set IDs that are currently in the cover
    std::vector<int> set_cover;
    
    // Helper functions (private)
    void applyUpdate(const Update& update);
    void addSetToCover(int s);  // Add set s to set_cover vector
    void removeSetFromCover(int s);  // Remove set s from set_cover vector
    
public:
    // Constructor: initialize data structures
    Algorithm2(int n, int m, double eps);
    
    // Process a single update from update sequence
    // update_row: [operation, element_id, ...sets] (sets only for insertions)
    void updateElement(const std::vector<int>& update_row);
    
    // Get recourse for the last update (O(1))
    int getRecourse() const { return last_update_recourse; }
    
    // Get set cover size after the last update (O(1))
    int getSetCoverSize() const { return last_set_cover_size; }
    
    // Get set cover as const reference (0-indexed, O(1), no copy)
    // Returns a const reference to the internal set_cover vector (0-indexed)
    const std::vector<int>& getSetCoverRef() const { return set_cover; }
};

#endif // ALG2_H