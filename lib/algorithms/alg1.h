#ifndef ALG1_H
#define ALG1_H

#include <vector>
#include <unordered_map>

// Forward declarations
struct Update;

class Algorithm1 {
private:
    // Parameters
    int n;
    int m;
    double beta;
    double eps;
    
    // Data structures
    std::vector<std::vector<int>> adj;
    std::vector<bool> in_sc;  // true if set is in the set cover
    int sc_size;
    std::unordered_map<int, int> element_id_to_index;
    
    // Batch management
    int updates_since_reset;
    int batch_size;
    
    // Per-update tracking (stored after each update)
    int last_update_recourse;
    int last_set_cover_size;
    
    // Set cover vector (maintained for O(1) access)
    // Stores 0-indexed set IDs that are currently in the cover
    std::vector<int> set_cover;
    std::vector<int> prev_set_cover;  // Previous set cover (for recourse calculation after reset)
    
    // Helper functions (private)
    void applyUpdate(const Update& update);
    void performReset();
    void addSetToCover(int s);
    void removeSetFromCover(int s);
    int getBatchSize(int sc_size, double eps);
    
public:
    // Constructor: initialize data structures
    Algorithm1(int n, int m, double eps);
    
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

#endif // ALG1_H