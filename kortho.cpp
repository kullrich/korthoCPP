#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
#include <zlib.h>
#include "include/ketopt.h"
#include "include/kseq.h"
#include "include/quickpool.hpp"

//KSEQ_INIT(gzFile, gzread)
KSEQ_INIT(int, read)

static void seq_to_map(
    std::vector<std::map<std::string, int>>& q_kmers_map,
    std::vector<std::vector<std::string>>& q_kmers,
    std::vector<std::vector<int>>& q_kmers_counts,
    const int k,
    std::string &seq,
    const int seq_idx) {
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    std::map<std::string, int> seq_kmer_map;
    std::map<std::string, int> seq_kmer_counts_map;
    std::vector<std::string> seq_kmers;
    std::vector<int> seq_kmers_counts_sorted;
    if (seq.length() >= k) {
        for (size_t i = 0; i < seq.length() - k + 1; ++i) {
            auto kmer = seq.substr(i, k);
            ++seq_kmer_counts_map[kmer];
            seq_kmer_map[kmer] = seq_idx;
        }
        std::map<std::string, int, std::less<std::string>> seq_kmer_map_sorted(seq_kmer_map.begin(), seq_kmer_map.end());
        std::map<std::string, int, std::less<std::string>> seq_kmer_counts_map_sorted(seq_kmer_counts_map.begin(), seq_kmer_counts_map.end());
        seq_kmers.reserve(seq_kmer_map_sorted.size());
        for (const auto& pair : seq_kmer_map_sorted) {
            seq_kmers.push_back(pair.first);
        }
        seq_kmers_counts_sorted.reserve(seq_kmer_counts_map_sorted.size());
        for (const auto& pair : seq_kmer_counts_map_sorted) {
            seq_kmers_counts_sorted.push_back(pair.second);
        }
        q_kmers_map.push_back(seq_kmer_map_sorted);
        q_kmers.push_back(seq_kmers);
        q_kmers_counts.push_back(seq_kmers_counts_sorted);
    } else {
        std::map<std::string, int> seq_kmer_map_sorted;
        q_kmers_map.push_back(seq_kmer_map_sorted);
        q_kmers.push_back(seq_kmers);
        q_kmers_counts.push_back(seq_kmers_counts_sorted);
    }
}

static void count_file(
        std::vector<std::map<std::string, int>>& q_kmers_map,
        std::vector<std::vector<std::string>>& q_kmers,
        std::vector<std::vector<int>>& q_kmers_counts,
        std::vector<std::string>& seqnames,
        const char *fn,
        int k) {
    //gzFile fp = gzopen(fn, "r");
    FILE* fp;
    fp = fopen(fn, "r");
    kseq_t *seq = kseq_init(fileno(fp));
    int counter = -1;
    while(kseq_read(seq) >= 0)
    {
        counter += 1;
        std::string cur_seq;
        cur_seq = (seq->seq.s);
        seq_to_map(q_kmers_map, q_kmers, q_kmers_counts, k, cur_seq, counter);
        std::string cur_name;
        cur_name = (seq->name.s);
        seqnames.push_back(cur_name);
    }
    kseq_destroy(seq);
}

std::map<std::string, std::vector<int>> combineMapsParallel(
    std::vector<std::map<std::string, int>>& maps,
    quickpool::ThreadPool& pool) {
    std::map<std::string, std::vector<int>> combinedMap;
    std::mutex mutex; // Mutex for thread-safe access to combinedMap
    // Define a lambda function to process each map
    auto processMap = [&](const std::map<std::string, int>& map) {
        for (const auto& pair : map) {
            std::lock_guard<std::mutex> lock(mutex); // Lock the mutex for thread-safe access
            combinedMap[pair.first].push_back(pair.second);
        }
    };
    pool.parallel_for_each(maps, [&] (std::map<std::string, int>& map) {
        processMap(map);
    });
    pool.wait();
    return combinedMap;
}

std::map<std::string, std::vector<int>> getSubsetMap(
    const std::map<std::string, std::vector<int>>& originalMap,
    int n) {
    std::map<std::string, std::vector<int>> subsetMap;
    int count = 0;
    for (const auto& pair : originalMap) {
        if (count >= n) {
            break;
        }
        subsetMap.insert(pair);
        count++;
    }
    return subsetMap;
}

static std::vector<std::vector<int>> findPairs_sparse(
        const std::map<std::string, std::vector<int>>& kmerMap_short,
        const std::map<std::string, std::vector<int>>& kmerMap_long,
        const int n) {
    auto it1 = kmerMap_short.begin();
    auto it2 = kmerMap_long.begin();
    std::vector<std::vector<int>> kmerMap_i_j_hits(n);
    while (it1 != kmerMap_short.end() && it2 != kmerMap_long.end()) {
        if (it1->first < it2->first) {
            ++it1;
        } else if (it2->first < it1->first) {
            ++it2;
        } else {
            for (const int& i : it1->second) {
                for (const int& j : it2->second) {
                    kmerMap_i_j_hits[i].push_back(j);
                }
            }
            ++it1;
            ++it2;
        }
    }
    std::vector<std::vector<int>> result(n);
    for (int i = 0; i < n; ++i) {
        if (!kmerMap_i_j_hits[i].empty()) {
            std::set<int> kmerMap_i_j_hits_set(kmerMap_i_j_hits[i].begin(), kmerMap_i_j_hits[i].end());
            std::vector<int> result_i(kmerMap_i_j_hits_set.size());
            result[i] = std::vector<int>(kmerMap_i_j_hits_set.begin(), kmerMap_i_j_hits_set.end());
        }
    }
    return result;
}

std::vector<int> findPairs(
        const std::map<std::string, int>& kmerMap_i,
        const std::map<std::string, std::vector<int>>& kmerMap) {
    auto it1 = kmerMap_i.begin();
    auto it2 = kmerMap.begin();
    std::vector<int> kmerMap_i_hits;
    while (it1 != kmerMap_i.end() && it2 != kmerMap.end()) {
        if (it1->first < it2->first) {
            ++it1;
        } else if (it2->first < it1->first) {
            ++it2;
        } else {
            for (const int& hit : it2->second) {
                kmerMap_i_hits.push_back(hit);
            }
            ++it1;
            ++it2;
        }
    }
    std::vector<int> result;
    if (!kmerMap_i_hits.empty()) {
        std::set<int> kmerMap_i_hits_set(kmerMap_i_hits.begin(), kmerMap_i_hits.end());
        //std::vector<int> result(kmerMap_i_hits_set.size());
        result = std::vector<int>(kmerMap_i_hits_set.begin(), kmerMap_i_hits_set.end());
        return result;
    }
    return result;
}

static std::tuple<std::vector<int>, std::vector<int>> flatten2tuple(
        const std::vector<std::vector<int>>& vec) {
    std::vector<int> result_i;
    std::vector<int> result_j;
    for (int i = 0; i < vec.size(); ++i) {
        for (int j = 0; j < vec[i].size(); ++j) {
            result_i.push_back(i);
            result_j.push_back(vec[i][j]);
        }
    }
    return std::make_tuple(result_i, result_j);
}

std::vector<std::vector<int>> transposeTQ(
        const std::vector<std::vector<int>>& matrix,
        const int n) {
    std::vector<std::vector<int>> transpose_matrix(n);
    size_t numRows = matrix.size();
    for (int row = 0; row < numRows; ++row) {
        if (!matrix[row].empty()) {
            for (size_t col = 0; col < matrix[row].size(); ++col) {
                transpose_matrix[matrix[row][col]].push_back(row);
            }
        }
    }
    std::vector<std::vector<int>> result(n);
    for (int i = 0; i < n; ++i) {
        if (!transpose_matrix[i].empty()) {
            std::set<int> result_set(transpose_matrix[i].begin(), transpose_matrix[i].end());
            std::vector<int> result_i(result_set.size());
            result[i] = std::vector<int>(result_set.begin(), result_set.end());
        }
    }
    return result;
}

double get_jaccard_value(
        const int common,
        const int uncommon1,
        const int uncommon2) {
    double jaccard_value;
    jaccard_value = ( static_cast<double>(common) ) / ( common + uncommon1 + uncommon2 );
    return jaccard_value;
}

double get_mash_value(
        const double jaccard_value,
        const int k) {
    double mash_value;
    mash_value = -( 1/static_cast<double>(k) ) * std::log( ( 2 * jaccard_value ) / ( 1 + jaccard_value ) );
    return mash_value;
}

double get_ani_value(
        const double mash_value) {
    double ani_value;
    ani_value = 1 - mash_value;
    return ani_value;
}

double get_sumdist_value(
        const int common1_sum,
        const int common2_sum,
        const int uncommon1_sum,
        const int uncommon2_sum) {
    double sumdist_value;
    sumdist_value = 1-( static_cast<double>(common1_sum) + common2_sum ) / ( common1_sum + common2_sum + uncommon1_sum + uncommon2_sum );
    return sumdist_value;
}

std::vector<double> getJaccardByIntegerVector(
        const std::vector<std::string>& seq_q_i_kmers_sorted,
        const std::vector<std::string>& seq_t_i_kmers_sorted,
        const std::vector<int>& seq_q_i_kmers_counts_sorted,
        const std::vector<int>& seq_t_i_kmers_counts_sorted,
        const int seq_q_i_idx,
        const int seq_t_i_idx,
        const int k) {
    std::vector<double> out(6);
    out[0] = static_cast<double>(seq_q_i_idx);
    out[1] = static_cast<double>(seq_t_i_idx);
    int i = 0;
    int j = 0;
    int common = 0;
    int uncommon1 = 0;
    int uncommon2 = 0;
    int common1_sum = 0;
    int common2_sum = 0;
    int uncommon1_sum = 0;
    int uncommon2_sum = 0;
    while (i < seq_q_i_kmers_sorted.size() && j < seq_t_i_kmers_sorted.size()) {
        if (seq_q_i_kmers_sorted[i] < seq_t_i_kmers_sorted[j]) {
            uncommon1 += 1;
            uncommon1_sum += seq_q_i_kmers_counts_sorted[i];
            ++i;
        } else if (seq_q_i_kmers_sorted[i] > seq_t_i_kmers_sorted[j]) {
            uncommon2 += 1;
            uncommon2_sum += seq_t_i_kmers_counts_sorted[j];
            ++j;
        } else {
            common += 1;
            common1_sum += seq_q_i_kmers_counts_sorted[i];
            common2_sum += seq_t_i_kmers_counts_sorted[j];
            ++i;
            ++j;
        }
    }
    for (; i < seq_q_i_kmers_sorted.size(); ++i) {
        uncommon1 += 1;
        uncommon1_sum += seq_q_i_kmers_counts_sorted[i];
    }
    for (; j < seq_t_i_kmers_sorted.size(); ++j) {
        uncommon2 += 1;
        uncommon2_sum += seq_t_i_kmers_counts_sorted[j];
    }
    //calculate jaccard, mash, ani, sumdist
    double jaccard_value;
    double mash_value;
    double ani_value;
    double sumdist_value;
    jaccard_value = get_jaccard_value(common, uncommon1, uncommon2);
    mash_value = get_mash_value(jaccard_value, k);
    ani_value = get_ani_value(mash_value);
    sumdist_value = get_sumdist_value(common1_sum, common2_sum, uncommon1_sum, uncommon2_sum);
    out[2] = jaccard_value;
    out[3] = mash_value;
    out[4] = ani_value;
    out[5] = sumdist_value;
    return out;
}

int main(int argc, char *argv[]) {
    int k = 6; // kmer size
    double m = 0.01; // min jaccard distance to report
    double s = 0.1; // sparse threshold to switch search strategy
    int n = 20; // number of kmer subset to get sparse threshold
    int p = 1; // number of parallel threads
    bool debug = false;
    const char *filename_q = nullptr; // Query fasta file
    const char *filename_t = nullptr; // Target fasta file
    const char *filename_o = "output.txt"; // Default output filename
    ketopt_t opt = KETOPT_INIT;
    int c;
    if (argc - opt.ind < 2) {
        fprintf(stderr, "Usage: korthoCPP [options] -q <query.fasta> -t <target.fasta>\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  -q FILE    query peptide fasta file\n");
        fprintf(stderr, "  -t FILE    target peptide fasta file\n");
        fprintf(stderr, "  -o FILE    output file to write jaccard distances\n");
        fprintf(stderr, "  -k INT     kmer length [default: 6]\n");
        fprintf(stderr, "  -m DOUBLE  min jaccard distance to report pair [default: 0.01]\n");
        fprintf(stderr, "  -s DOUBLE  sparse threshold to switch search strategy [defualt: 0.1]\n");
        fprintf(stderr, "  -n INT     number of kmers to check for sparse [default: 20]\n");
        fprintf(stderr, "  -p INT     number of threads [default: 1]\n");
        fprintf(stderr, "  -d         debug\n");
        return 1;
    }
    while ((c = ketopt(&opt, argc, argv, 1, "q:t:o:k:m:s:n:p:d", nullptr)) >= 0) {
        switch (c) {
            case 'd':
                debug = true;
                break;
            case 'k':
                k = std::stoi(opt.arg);
                break;
            case 'm':
                m = std::stod(opt.arg);
                break;
            case 'n':
                n = std::stoi(opt.arg);
                break;
            case 's':
                s = std::stod(opt.arg);
                break;
            case 'p':
                p = std::stoi(opt.arg);
                break;
            case 'q':
                filename_q = opt.arg;
                break;
            case 't':
                filename_t = opt.arg;
                break;
            case 'o':
                filename_o = opt.arg;
                break;
            default:
                std::cerr << "Unknown option: " << char(c) << std::endl;
                return 1;
        }
    }
    // Check if both filenames are provided
    if (filename_q == nullptr || filename_t == nullptr) {
        std::cerr << "Error: Both filenames ('-q' and '-t') must be provided." << std::endl;
        return 1;
    }
    std::map<std::string, std::vector<int>> QkmerMap;
    std::map<std::string, std::vector<int>> TkmerMap;
    std::vector<std::map<std::string, int>> QkmerMap_n;
    std::vector<std::map<std::string, int>> TkmerMap_m;
    std::vector<std::vector<std::string>> Qkmers_n;
    std::vector<std::vector<int>> Qkmers_n_counts;
    std::vector<std::vector<std::string>> Tkmers_m;
    std::vector<std::vector<int>> Tkmers_m_counts;
    std::vector<std::string> Qnames;
    std::vector<std::string> Tnames;
    int seq_q_n;
    int seq_t_m;
    bool use_sparse=false;
    // create Thread pool
    quickpool::ThreadPool pool(p);
    // create kmers vector
    auto start_kmers_vec = std::chrono::steady_clock::now();
    auto workerQ = [&QkmerMap_n, &Qkmers_n, &Qkmers_n_counts, &Qnames] (
        const char *fn,
        int k) {
        count_file(QkmerMap_n, Qkmers_n, Qkmers_n_counts, Qnames, fn, k);
    };
    auto workerT = [&TkmerMap_m, &Tkmers_m, &Tkmers_m_counts, &Tnames] (
            const char *fn,
            int k) {
        count_file(TkmerMap_m, Tkmers_m, Tkmers_m_counts, Tnames, fn, k);
    };
    pool.push(workerQ, filename_q, k);
    pool.push(workerT, filename_t, k);
    pool.wait();
    seq_q_n = static_cast<int>(Qkmers_n.size());
    seq_t_m = static_cast<int>(Tkmers_m.size());
    auto end_kmers_vec = std::chrono::steady_clock::now();
    auto duration_kmers_vec = std::chrono::duration_cast<std::chrono::milliseconds>(end_kmers_vec - start_kmers_vec);
    if (debug) {
        std::cout << "Time taken: kmer extraction " << duration_kmers_vec.count() << " milliseconds" << std::endl;
        std::cout << "number of sequences Q: " << seq_q_n << std::endl;
        std::cout << "number of sequences T: " << seq_t_m << std::endl;
    }
    // create combined kmerMaps
    auto start_QkmerMap = std::chrono::steady_clock::now();
    //QkmerMap = combineMaps(QkmerMap_n);
    QkmerMap = combineMapsParallel(QkmerMap_n, pool);
    //for (const auto& pair : QkmerMap) {
        //std::cout << "QkmerMap kmer sorted: " << pair.first << std::endl;
    //}
    auto end_QkmerMap = std::chrono::steady_clock::now();
    auto duration_QkmerMap = std::chrono::duration_cast<std::chrono::milliseconds>(end_QkmerMap - start_QkmerMap);
    if (debug) {
        std::cout << "Time taken: QkmerMap creation " << duration_QkmerMap.count() << " milliseconds" << std::endl;
        std::cout << "number of kmers QkmerMap: " << QkmerMap.size() << std::endl;
    }
    auto start_TkmerMap = std::chrono::steady_clock::now();
    //TkmerMap = combineMaps(TkmerMap_m);
    TkmerMap = combineMapsParallel(TkmerMap_m, pool);
    //for (const auto& pair : TkmerMap) {
        //std::cout << "TkmerMap kmer sorted: " << pair.first << std::endl;
    //}
    auto end_TkmerMap = std::chrono::steady_clock::now();
    auto duration_TkmerMap = std::chrono::duration_cast<std::chrono::milliseconds>(end_TkmerMap - start_TkmerMap);
    if (debug) {
        std::cout << "Time taken: TkmerMap creation " << duration_TkmerMap.count() << " milliseconds" << std::endl;
        std::cout << "number of kmers TkmerMap: " << TkmerMap.size() << std::endl;
    }
    // check sparse
    auto start_sparseCheck = std::chrono::steady_clock::now();
    std::map<std::string, std::vector<int>> subsetQkmerMap;
    subsetQkmerMap = getSubsetMap(QkmerMap, n);
    std::map<std::string, std::vector<int>> subsetTkmerMap;
    subsetTkmerMap = getSubsetMap(QkmerMap, n);
    std::vector<std::vector<int>> sparse_pairs(seq_q_n);
    sparse_pairs = findPairs_sparse(subsetQkmerMap, subsetTkmerMap, seq_q_n);
    std::tuple<std::vector<int>, std::vector<int>> sparsePairs;
    sparsePairs = flatten2tuple(sparse_pairs);
    std::vector<int> sparsePairs_i = std::get<0>(sparsePairs);
    std::vector<int> sparsePairs_j = std::get<1>(sparsePairs);
    double sparseValue = (static_cast<double>(sparsePairs_i.size())/(seq_q_n*seq_t_m));
    if (sparseValue<s) {
        use_sparse = true;
    }
    if (debug) {
        std::cout << "number of sparse candidate pairs: " << sparsePairs_i.size();
        std::cout << std::endl;
        if (use_sparse) {
            std::cout << "sparse value: " << sparseValue << " < sparse threshold: " << s << " >>> search strategy one vs one";
            std::cout << std::endl;
        } else {
            std::cout << "sparse value: " << sparseValue << " > sparse threshold: " << s << " >>> search strategy one vs many";
            std::cout << std::endl;
        }
    }
    auto end_sparseCheck = std::chrono::steady_clock::now();
    auto duration_sparseCheck = std::chrono::duration_cast<std::chrono::milliseconds>(end_sparseCheck - start_sparseCheck);
    if (debug) {
        std::cout << "Time taken: check sparse threshold " << duration_sparseCheck.count() << " milliseconds" << std::endl;
    }
    //get candidate pairs
    auto start_getCandidates = std::chrono::steady_clock::now();
    std::vector<std::vector<int>> Q_T_pairs(seq_q_n);
    if (QkmerMap.size() < TkmerMap.size()) {
        if (use_sparse) {
            Q_T_pairs = findPairs_sparse(QkmerMap, TkmerMap, seq_q_n);
        } else {
            //for (int i = 0; i < seq_q_n; ++i) {
            pool.parallel_for(0, seq_q_n, [&] (int i) {
                Q_T_pairs[i] = findPairs(QkmerMap_n[i], TkmerMap);
            });
            pool.wait();
            //}
        }
    } else {
        std::vector<std::vector<int>> T_Q_pairs(seq_t_m);
        if (use_sparse) {
            T_Q_pairs = findPairs_sparse(TkmerMap, QkmerMap, seq_t_m);
        } else {
            //for (int j = 0; j < seq_t_m; ++j) {
            pool.parallel_for(0, seq_t_m, [&] (int j) {
                T_Q_pairs[j] = findPairs(TkmerMap_m[j], QkmerMap);
            });
            pool.wait();
            //}
        }
        //transpose T_Q_hits to get Q_T_hits
        Q_T_pairs = transposeTQ(T_Q_pairs, seq_q_n);
    }
    std::tuple<std::vector<int>, std::vector<int>> candidatePairs;
    candidatePairs = flatten2tuple(Q_T_pairs);
    std::vector<int> candidatePairs_ni = std::get<0>(candidatePairs);
    std::vector<int> candidatePairs_mj = std::get<1>(candidatePairs);
    int candidatePairs_ni_mj_size;
    candidatePairs_ni_mj_size = static_cast<int>(candidatePairs_ni.size());
    pool.done();
    auto end_getCandidates = std::chrono::steady_clock::now();
    auto duration_getCandidates = std::chrono::duration_cast<std::chrono::milliseconds>(end_getCandidates - start_getCandidates);
    if (debug) {
        std::cout << "Time taken: candidate pairs creation " << duration_getCandidates.count() << " milliseconds" << std::endl;
        std::cout << "number of candidate pairs: " << candidatePairs_ni_mj_size;
        std::cout << std::endl;
    }
    //calculate distances for each candidate pair
    auto start_calcDist = std::chrono::steady_clock::now();
    std::vector<std::vector<double>> Qni_Tmj_distances(
        candidatePairs_ni_mj_size,
        std::vector<double>(6, 0.0));
    pool.parallel_for(0, candidatePairs_ni_mj_size, [&] (int ni_mj) {
        int Qni = candidatePairs_ni[ni_mj];
        int Tmj = candidatePairs_mj[ni_mj];
        Qni_Tmj_distances[ni_mj]=getJaccardByIntegerVector(
            Qkmers_n[Qni],
            Tkmers_m[Tmj],
            Qkmers_n_counts[Qni],
            Tkmers_m_counts[Tmj],
            Qni,
            Tmj,
            k);
    });
    pool.wait();
    auto end_calcDist = std::chrono::steady_clock::now();
    auto duration_calcDist = std::chrono::duration_cast<std::chrono::milliseconds>(end_calcDist - start_calcDist);
    if (debug) {
        std::cout << "Time taken: distance calculation " << duration_calcDist.count() << " milliseconds" << std::endl;
    }
    // min jaccard filter
    std::vector<std::string> out_qname;
    std::vector<std::string> out_tname;
    std::vector<double> out_jaccard;
    std::vector<double> out_mash;
    std::vector<double> out_ani;
    std::vector<double> out_sumdist;
    std::ofstream outputFile(filename_o);
    if (!outputFile) {
        std::cerr << "Error: Failed to open output file." << std::endl;
        return 0;
    }
    outputFile << "qname\ttname\tjaccard\tmash\tani\tsumdist\n";
    for (const auto& distances : Qni_Tmj_distances) {
        if ( distances[2] > m ) {
            outputFile << Qnames[static_cast<size_t>(distances[0])] << '\t'
                       << Tnames[static_cast<size_t>(distances[1])] << '\t'
                       << distances[2] << '\t'
                       << distances[3] << '\t'
                       << distances[4] << '\t'
                       << distances[5] << '\n';
        }
    }
    outputFile.close();
    return 0;
}
