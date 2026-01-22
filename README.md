# Dynamic Set Cover

This repository contains the implementation for dynamic set cover algorithms from the paper "Engineering Algorithms for Dynamic Greedy Set Cover".

## Table of Contents

1. [System Requirements](#system-requirements)
2. [Downloading and Extracting](#downloading-and-extracting)
3. [Dataset Files](#dataset-files)
4. [Dependencies](#dependencies)
5. [Compilation](#compilation)
6. [Input Data Format](#input-data-format)
7. [Running Experiments](#running-experiments)
8. [Output Format](#output-format)
9. [License](#license)

## System Requirements

- **Operating System**: Linux (tested on Ubuntu 22.04) or WSL (Windows Subsystem for Linux)
- **Compiler**: GCC 7.0 or later with C++17 support
- **Build System**: CMake 3.10 or later
- **Memory**: Sufficient RAM for dataset sizes
- **Disk Space**: ~500MB for code and build artifacts, plus space for results

## Downloading and Extracting

1. Clone the repository from GitHub:
   ```bash
   git clone https://github.com/DynamicSetCover/DynamicSetCover.git
   ```
2. Navigate to the project directory:
   ```bash
   cd DynamicSetCover
   ```

## Dataset Files

The repository includes a subset of the benchmark datasets for testing and validation. The complete dataset collection is available separately.

### Included Datasets

This repository contains the first 10 datasets (dataset001.hgr through dataset010.hgr) in the `data/` directory. These are sufficient for initial testing and verification of the implementation.

### Complete Dataset Collection

The full set of 120 datasets (datasets 001-120) is available for download from Zenodo:

**Zenodo Repository:** https://zenodo.org/records/18268181

After downloading the complete dataset archive, extract it to the `data/` directory to replace or supplement the included datasets.

### Dataset Metadata Files

The `data/` directory also contains two metadata files:

- **`dataset_mapping.txt`**: Maps the standardized dataset filenames (dataset001.hgr, dataset002.hgr, etc.) to their original source file names. This is useful for referencing the original datasets and understanding their provenance.

- **`dataset_summary.txt`**: Contains statistical information about each dataset.

These files provide useful metadata for understanding the characteristics of each dataset in the benchmark suite.

## Dependencies

The project uses only standard C++ libraries and includes the Argtable3 library (for command-line parsing) as a bundled dependency. No external C++ libraries need to be installed.

## Compilation

### Step 1: Create Build Directory

```bash
mkdir build
cd build
```

### Step 2: Configure with CMake

```bash
cmake ..
```

This will configure the build system. By default, it builds in Release mode with optimizations.

### Step 3: Compile

```bash
make
```

The compilation will create the executable at `build/app/dynsc_benchmark`.

### Step 4: Verify Build

Check that the executable was created:

```bash
ls -lh build/app/dynsc_benchmark
```

You should see the executable file. You can also test it:

```bash
./build/app/dynsc_benchmark --help
```

This should display the usage information.

## Input Data Format

The input dataset files are in `.hgr` format and contain a sequence of dynamic updates (insertions and deletions of elements).

### File Structure

Each dataset file consists of:
1. **Header line**: Contains metadata about the dataset
2. **Update lines**: Each subsequent line represents a single update operation

### Header Line Format

The first line of each file is a header that starts with `#` and contains four space-separated integers:

```
# k n m f
```

Where:
- **k**: Total number of updates in the sequence
- **n**: Maximum number of elements at any point during the sequence
- **m**: Number of sets 
- **f**: Frequency (maximum number of sets containing an element)

**Example:**
```
# 5082 254 31022 969
```
This indicates 5082 updates, up to 254 elements, 31022 sets, and frequency 969.

### Update Line Format

Each line after the header represents a single update operation. The format is:

```
<operation> <element_id> [set_ids...]
```

Where:
- **operation**: `0` = insert element, `1` = delete element
- **element_id**: The ID of the element being inserted or deleted
- **set_ids**: (Only for insertions) Space-separated list of set IDs that contain this element. This field is omitted for deletions.

**Examples:**

Insertion:
```
0 0 1 2 3 4 5 6 7 8 9 10 11 12 13
```
This inserts element 0, which is contained in sets 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, and 13.

Deletion:
```
1 2
```
This deletes element 2.

### Complete Example

```
# 5082 254 31022 969
0 0 1 2 3 4 5 6 7 8 9 10 11 12 13
0 1 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41
0 2 42 43 44
1 2
0 3 45 46 47 48 49 50 51
```

This example shows:
- Header: 5082 updates, 254 maximum elements, 31022 sets, 969 frequency
- Update 1: Insert element 0 (contained in sets 1-13)
- Update 2: Insert element 1 (contained in sets 14-41)
- Update 3: Insert element 2 (contained in sets 42-44)
- Update 4: Delete element 2
- Update 5: Insert element 3 (contained in sets 45-51)

## Running Experiments

### Command-Line Usage

To run a single experiment:

```bash
./build/app/dynsc_benchmark --algorithm <1-4> --epsilon <value> --dataset <path_to_dataset.hgr> --runs <number_of_runs> --dump_freq <frequency>
```

**Parameters:**
- `--algorithm`: Algorithm ID (1, 2, 3, or 4)
  - Algorithm 1: Robust Algorithm
  - Algorithm 2: Local Algorithm
  - Algorithm 3: Partial Algorithm
  - Algorithm 4: Global Algorithm
- `--epsilon`: Trade-off parameter (epsilon = beta - 1, so if needed beta = 1.5 use epsilon = 0.5)
- `--dataset`: Path to input dataset file (e.g., `data/dataset001.hgr`)
- `--runs`: Number of benchmark runs (default: 1). When multiple runs are specified, separate output files are created for each run.
- `--dump_freq`: Set cover output frequency in steps (default: 0 = never). If set to a positive value N, the sets in the set cover are written every N updates. Use 1 for outputting the set cover every update (useful for verification) but this creates very large output files.

**Example:**
```bash
./build/app/dynsc_benchmark --algorithm 1 --epsilon 0.5 --dataset data/dataset001.hgr --runs 5 --dump_freq 10
```

This will run the Robust Algorithm with epsilon 0.5 (beta = 1.5) on dataset001.hgr for 5 runs, creating 5 separate output files, containing the maintained set cover every 10 updates.

## Output Format

### Output File Naming

Each experiment run creates a result file with the naming pattern:
```
dataset<NNN>_alg<algorithm>_eps<epsilon>_run<run_number>.txt
```

Example: `dataset001_alg1_eps0.5_run1.txt`

The output files are created in the `results/` directory (created in the current working directory where the command is executed). If you run the command from the `build/` directory, results will be in `build/results/`.

### File Format

The output file begins with header lines (prefixed with `#`) that provide summary statistics and metadata, followed by the per-update results.

#### Header Lines

The file starts with three types of header lines:

1. **Summary statistics** (2 lines):
   - First line: Description of summary fields
   - Second line: Actual summary values (max_time, avg_time, max_recourse, avg_recourse, max_set_cover_size, avg_set_cover_size)

2. **Data column description**: Describes what each column in the data lines represents

3. **Set cover frequency**: Indicates the `--dump_freq` setting used (0 means set cover is never included in output)

#### Data Lines

Following the headers, each line in the output file contains (space-separated):
1. **Update time**: Time taken for this update in nanoseconds
2. **Recourse**: Number of sets added/removed in the set cover during this update
3. **Set cover size**: Current size of the set cover
4. **(Optional) Set cover**: If `dump_freq > 0` and this update is a multiple of `dump_freq`, the sets in the set cover are included (space-separated)

**Note**: The update number is implicit (line number minus header lines). The first data line corresponds to the first update, the second to the second update, and so on.

**Example output:**
```
# Summary: max_time(ns) avg_time(ns) max_recourse avg_recourse max_set_cover_size avg_set_cover_size
# 54321 34567 5 2 15 12
# Update results: time(nanoseconds) recourse set_cover_size set_cover
# Set cover output frequency: 3 (set cover included every 3 updates)
12345 2 10
23456 1 11
34567 0 11 1 3 5 7 9 11 13 15 17 19 21
```

In this example:
- The maximum (worst-case) update time was 54321(ns), the average (amortized) update time was 34567(ns)
- The maximum (worst-case) recourse was 5, the average (amortized) recourse was 2
- The maximum (worst-case) set cover size was 15, the average (amortized) set cover size was 12
- First update: took 12345 nanoseconds, recourse=2, set cover size=10
- Second update: took 23456 nanoseconds, recourse=1, set cover size=11
- Third update: took 34567 nanoseconds, recourse=0, set cover size=11, and the set cover is the collection of sets {1,3,5,7,9,11,13,15,17,19,21}

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


This project is currently under double-blind review. If you publish results using these algorithms, or if you have any other questions, please contact our anonymous email at: dynamicsetcover@gmail.com
