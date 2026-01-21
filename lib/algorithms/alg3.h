#ifndef ALG3_H
#define ALG3_H

#include <vector>
#include <unordered_map>
#include <functional>
#include <unordered_set>

// Forward declarations
struct Update;

class Algorithm3 {
private:
    // Parameters
    int n;
    int m;
    int num_levels;
    double beta;
    double eps;
    
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
    std::vector<std::vector<int>> level_to_elements;
    
    // Counter-based algorithm structures
    std::vector<int> S;  // Set counter: S(k) = number of sets covering at level k
    std::vector<double> D;  // Dirt counter: D(k) = amount of dirt at level k
    
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
    
    // Helper functions for the algorithm (declared but implemented as regular functions in .cpp)
    
public:
    // Constructor: initialize data structures
    Algorithm3(int n, int m, double eps);
    
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

#endif // ALG3_H