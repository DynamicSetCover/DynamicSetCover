#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <chrono>
#include <string>
#include "../lib/DataTypes.h"
#include "../extern/argtable3.h"   
#include "../lib/io/dynsc_io.h"  
#include "../lib/algorithms/dynsc_ctrl.h" 
#include "../lib/algorithms/alg1.h" 
#include "../lib/algorithms/alg2.h"
#include "../lib/algorithms/alg3.h"
#include "../lib/algorithms/alg4.h"


struct Parameters {
    std::string dataset_filename;
    int algorithm; 
    double eps;
    int num_runs;
    int set_cover_output_frequency;
    int k, m, n; 
};


bool parse_arguments(int argc, char** argv, Parameters& config) {
    struct arg_int *alg = arg_int1(NULL, "algorithm", "<1-4>", "Algorithm ID to run (1 to 4)");
    struct arg_dbl *eps = arg_dbl0(NULL, "epsilon", "<n>", "Trade-off parameter epsilon (default: 0.5)");
    struct arg_int *runs = arg_int0(NULL, "runs", "<n>", "Number of benchmark runs (default: 1)");
    struct arg_int *dump_freq = arg_int0(NULL, "dump_freq", "<n>", "Set cover output frequency (steps) (default: 0 = never)");
    struct arg_file *dataset = arg_file1(NULL, "dataset", "<file>", "Path to the input dataset file.");
    struct arg_lit *help = arg_lit0(NULL, "help", "Display this help message.");
    struct arg_end *end = arg_end(20);

    void* argtable[] = {alg, eps, runs, dump_freq, dataset, help, end};
    int nerrors;
    
    // Set default values
    eps->dval[0] = 0.5;
    runs->ival[0] = 1; 
    dump_freq->ival[0] = 0; 

    nerrors = arg_parse(argc, argv, argtable);

    if (help->count > 0 || nerrors > 0) {
        if (nerrors > 0) {
            std::cerr << "Syntax Errors detected:" << std::endl;
            arg_print_errors(stdout, end, argv[0]);
        }
        std::cout << "Usage: " << argv[0] << std::endl;
        arg_print_glossary(stdout, argtable, "  %-25s %s\n");
        arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
        return false;
    }

    // 3. MAP Argtable3 VALUES TO Parameters STRUCT
    config.algorithm = alg->ival[0];
    config.eps = eps->dval[0];
    config.num_runs = runs->ival[0];
    config.set_cover_output_frequency = dump_freq->ival[0];
    config.dataset_filename = dataset->filename[0];

    // Input validation
    if (config.algorithm < 1 || config.algorithm > 4) {
        std::cerr << "Error: Algorithm must be 1, 2, 3, or 4. Got: " << config.algorithm << std::endl;
        arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
        return false;
    }

    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
    return true;
}


int main(int argc, char* argv[]) {
    // 1. Initialize and Parse using Argtable3
    Parameters params;
    if (!parse_arguments(argc, argv, params)) {
        return 1; // Exit on error or help request
    }

    // Initialize I/O handler
    auto updateIO = dynsc_io();
    auto start = std::chrono::steady_clock::now();
    
    // Initialize the update sequence
    std::vector<std::vector<int>> update_sequence;
    
    // Read update sequence from file
    auto [k, m, n] = updateIO.readUpdateSequence(params.dataset_filename, update_sequence);
    auto stop = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time taken to read update sequence: " << duration.count() << " nanoseconds" << std::endl;

    // Check if reading was successful
    if (k == -1 || m == -1 || n == -1) {
        std::cerr << "Error: Failed to read update sequence from " << params.dataset_filename << std::endl;
        return 1;
    }
    
    std::cout << "Loaded dataset: k=" << k << ", m=" << m << ", n=" << n << std::endl;
    std::cout << "Running algorithm " << params.algorithm << " on dataset " << params.dataset_filename << " with eps=" << params.eps << " for " << params.num_runs << " runs" << std::endl;

    // Run the algorithm multiple times
    for (int run = 0; run < params.num_runs; ++run) {
        // Run the dynamic algorithm and get results
        //std::vector<UpdateResult> results = run_dynamic_algorithm(update_sequence, params.algorithm, params.eps, k, m, n);
        std::vector<UpdateResult> results;

        switch (params.algorithm) {
            case 1:
                results = run_algorithm_loop<Algorithm1>(update_sequence, k, n, m, params.eps, params.set_cover_output_frequency > 0);
                break;
            case 2:
                results = run_algorithm_loop<Algorithm2>(update_sequence, k, n, m, params.eps, params.set_cover_output_frequency > 0);
                break;
            case 3:
                results = run_algorithm_loop<Algorithm3>(update_sequence, k, n, m, params.eps, params.set_cover_output_frequency > 0);
                break;
            case 4:
                results = run_algorithm_loop<Algorithm4>(update_sequence, k, n, m, params.eps, params.set_cover_output_frequency > 0);
                break;
            // Default case should be unreachable due to Argtable3 validation
            default: 
                std::cerr << "Internal Error: Algorithm ID not in range [1, 4]." << std::endl;
                return 1;
        }

        // Validate results
        std::cout << "Finished run number " << run + 1 << std::endl;
        auto start = std::chrono::steady_clock::now();

        // Create results directory if it doesn't exist
        std::filesystem::create_directories("results");
        
        // Generate output filename based on dataset and algorithm
        std::string dataset_name = params.dataset_filename;
        // Extract filename without path
        size_t last_slash = dataset_name.find_last_of("/\\");
        if (last_slash != std::string::npos) {
            dataset_name = dataset_name.substr(last_slash + 1);
        }
        // Remove extension if present
        size_t last_dot = dataset_name.find_last_of(".");
        if (last_dot != std::string::npos) {
            dataset_name = dataset_name.substr(0, last_dot);
        }
        
        std::ostringstream filename;
        filename << "results/" << dataset_name << "_alg" << params.algorithm << "_eps" << params.eps;
        if (params.num_runs > 1) {
            filename << "_run" << (run + 1);
        }
        filename << ".txt";
        
        // Write results to file
        std::ofstream outfile(filename.str());
        if (!outfile.is_open()) {
            std::cerr << "Error: Could not open output file " << filename.str() << std::endl;
            return 1;
        }
        
        // Calculate statistics
        if (results.empty()) {
            std::cerr << "Error: No results to write" << std::endl;
            return 1;
        }
        
        long long max_time = 0;
        long long total_time = 0;
        int max_recourse = 0;
        long long total_recourse = 0;
        int max_set_cover_size = 0;
        long long total_set_cover_size = 0;
        
        for (const auto& result : results) {
            // Time statistics
            if (result.time_nanoseconds > max_time) {
                max_time = result.time_nanoseconds;
            }
            total_time += result.time_nanoseconds;
            
            // Recourse statistics
            if (result.recourse > max_recourse) {
                max_recourse = result.recourse;
            }
            total_recourse += result.recourse;
            
            // Set cover size statistics
            if (result.set_cover_size > max_set_cover_size) {
                max_set_cover_size = result.set_cover_size;
            }
            total_set_cover_size += result.set_cover_size;
        }
        
        double avg_time = static_cast<double>(total_time) / results.size();
        double avg_recourse = static_cast<double>(total_recourse) / results.size();
        double avg_set_cover_size = static_cast<double>(total_set_cover_size) / results.size();
        
        // Write summary header
        outfile << "# Summary: max_time(ns) avg_time(ns) max_recourse avg_recourse max_set_cover_size avg_set_cover_size\n";
        outfile << "# " << max_time << " " << avg_time << " " 
                << max_recourse << " " << avg_recourse << " "
                << max_set_cover_size << " " << avg_set_cover_size << "\n";
        
        // Write data header
        outfile << "# Update results: time(nanoseconds) recourse set_cover_size set_cover\n";
        if (params.set_cover_output_frequency == 0) {
            outfile << "# Set cover output frequency: 0 (set cover never included)\n";
        } else {
            outfile << "# Set cover output frequency: " << params.set_cover_output_frequency << " (set cover included every " << params.set_cover_output_frequency << " updates)\n";
        }
        
        // Write results for each update
        for (size_t i = 0; i < results.size(); ++i) {
            const auto& result = results[i];
            outfile << result.time_nanoseconds << " " 
                    << result.recourse << " " 
                    << result.set_cover_size << " ";
            
            // Write set cover only every X updates (if X=0, never write; if X>0, write every X updates)
            if (params.set_cover_output_frequency > 0 && i % params.set_cover_output_frequency == 0) {
                // Write set cover (space-separated, convert from 0-indexed to 1-indexed)
                for (size_t j = 0; j < result.set_cover.size(); ++j) {
                    if (j > 0) outfile << " ";
                    outfile << (result.set_cover[j] + 1); // Convert from 0-indexed to 1-indexed
                }
            }
            // If set cover is not included, the line ends after set_cover_size
            outfile << "\n";
        }
        
        outfile.close();

        auto stop = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
        std::cout << "Time taken for post-processing: " << duration.count() << " nanoseconds" << std::endl;

        std::cout << "Results written to " << filename.str() << std::endl;
    }
    
    return 0;
}