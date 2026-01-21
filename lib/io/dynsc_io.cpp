#include "dynsc_io.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>

dynsc_io::dynsc_io() {}

    // Read update sequence from file
    // Returns tuple (k, m, n) where k is number of updates, m is number of sets, n is max number of elements
std::tuple<int, int, int> dynsc_io::readUpdateSequence(const std::string& filename, std::vector<std::vector<int>>& update_sequence) {
        // Resize update_sequence to 0
        update_sequence.resize(0);
        
        std::ifstream file(filename);
        
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return {-1, -1, -1};
        }
        
        std::string line;
        bool header_read = false;
        int k = -1, m = -1, n = -1;
        
        while (std::getline(file, line)) {
            // Trim leading and trailing whitespace
            size_t start = line.find_first_not_of(" \t\r\n");
            if (start == std::string::npos) continue; // Empty line
            size_t end = line.find_last_not_of(" \t\r\n");
            line = line.substr(start, end - start + 1);
            
            // Check if this is the header line (starts with #)
            if (!line.empty() && line[0] == '#') {
                if (header_read) {
                    std::cerr << "Warning: Multiple header lines found, ignoring: " << line << std::endl;
                    continue;
                }
                std::istringstream iss(line);
                char hash;
                int f; // We don't need f, but need to read it
                if (iss >> hash >> k >> n >> m >> f) {
                    header_read = true;
                    continue;
                } else {
                    std::cerr << "Error: Could not parse header line: " << line << std::endl;
                    file.close();
                    return {-1, -1, -1};
                }
            }
            
            // Skip if header not read yet
            if (!header_read) {
                std::cerr << "Error: Header line not found before data line: " << line << std::endl;
                file.close();
                return {-1, -1, -1};
            }
            
            std::istringstream iss(line);
            std::vector<int> update_line;
            
            // Read operation (0 = insert, 1 = delete)
            int operation;
            if (!(iss >> operation)) {
                std::cerr << "Error: Could not read operation from line: " << line << std::endl;
                continue;
            }
            update_line.push_back(operation);
            
            // Read element identifier
            int element_id;
            if (!(iss >> element_id)) {
                std::cerr << "Error: Could not read element_id from line: " << line << std::endl;
                continue;
            }
            update_line.push_back(element_id);
            
            // For insertions, read the sets
            if (operation == 0) {
                int set_id;
                while (iss >> set_id) {
                    update_line.push_back(set_id);
                }
            }
            
            update_sequence.push_back(std::move(update_line));
        }
        
        file.close();
        
        if (!header_read) {
            std::cerr << "Error: Header line not found in dataset file " << filename << std::endl;
            return {-1, -1, -1};
        }
        
        return {k, m, n};
    }