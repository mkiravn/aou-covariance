/*
 * grm_class_tool -- bin phenotype cross-products from a row-chunked sharded
 * GRM, keeping designated pair classes in their OWN bin sets.
 *
 * Parent-offspring and full-sib pairs share the same expected additive
 * relatedness (a_ij ~= 0.5) so the GRM cannot separate them, and they land in
 * the same bin. Given an external classification (e.g. from IBS0), this routes
 * each listed pair into a per-class copy of the bin array instead of pooling
 * them. Unlisted pairs go to class "other".
 *
 * That keeps every pair -- nothing is discarded -- so one pass yields the mean
 * cross-product for each class separately, and the pooled result is recoverable
 * by summing classes.
 *
 * Standalone: does not modify GRM-pairs/grm_bin_sharded/grm_shard_tool, whose
 * accumulate/merge remain the reference implementation for the unclassed case.
 * The row-range recovery and accumulator maths below are ported from it
 * deliberately unchanged, so classed and unclassed runs stay comparable.
 *
 * Subcommands:
 *   accumulate  -- read one shard, accumulate per-(class, bin, jackknife-block)
 *   merge       -- sum accumulator files, emit per-class bin means/SDs/jackknife SEs
 *   ranges      -- print the row range a given shard should cover
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

static const char* OTHER_CLASS = "other";

struct Bin {
    double left;
    double right;
};

struct Acc {
    double sum = 0.0;
    double sum_sq = 0.0;
    std::uint64_t n = 0;

    void add(double x) {
        sum += x;
        sum_sq += x * x;
        n += 1;
    }
    void add(const Acc& other) {
        sum += other.sum;
        sum_sq += other.sum_sq;
        n += other.n;
    }
};

// ---------------------------------------------------------------------------
// Readers
// ---------------------------------------------------------------------------

static bool is_missing_token(const std::string& s) {
    return s == "NA" || s == "NaN" || s == "nan" || s == ".";
}

static double parse_double_or_nan(const std::string& s) {
    if (is_missing_token(s)) return std::numeric_limits<double>::quiet_NaN();
    try {
        std::size_t idx = 0;
        double x = std::stod(s, &idx);
        return idx == s.size() ? x : std::numeric_limits<double>::quiet_NaN();
    } catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

struct IdTable {
    std::vector<std::string> fid, iid;
    std::unordered_map<std::string, std::size_t> index_of_iid;
};

// IID alone is the key: plink writes FID = "0" for every row of a .grm.id built
// without pedigree, so (FID, IID) carries no more information than IID and
// forces every downstream file to reproduce plink's convention exactly.
// Uniqueness is checked rather than assumed.
static IdTable read_ids(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Could not open grm.id file: " + path);

    IdTable t;
    std::string fid, iid;
    while (in >> fid >> iid) {
        if (t.index_of_iid.count(iid)) {
            throw std::runtime_error("Duplicate IID in " + path + ": " + iid +
                                     " -- IID must be unique to key pair classes on it");
        }
        t.index_of_iid[iid] = t.fid.size();
        t.fid.push_back(fid);
        t.iid.push_back(iid);
    }
    if (t.fid.empty()) throw std::runtime_error("No IDs found in grm.id file: " + path);
    return t;
}

static std::vector<double> read_pheno_aligned(const std::string& path, const IdTable& ids) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Could not open phenotype file: " + path);

    std::unordered_map<std::string, double> pheno;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        std::string fid, iid, ystr;
        if (!(iss >> fid >> iid >> ystr)) continue;
        if (fid == "FID" && iid == "IID") continue;
        pheno[iid] = parse_double_or_nan(ystr);
    }

    std::vector<double> y(ids.iid.size(), std::numeric_limits<double>::quiet_NaN());
    for (std::size_t i = 0; i < ids.iid.size(); ++i) {
        auto it = pheno.find(ids.iid[i]);
        if (it != pheno.end()) y[i] = it->second;
    }
    return y;
}

static std::vector<Bin> read_bins(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Could not open bins file: " + path);

    std::vector<Bin> bins;
    double a, b;
    while (in >> a >> b) {
        if (!(a < b)) throw std::runtime_error("Invalid bin with left >= right");
        bins.push_back({a, b});
    }
    if (bins.empty()) throw std::runtime_error("No bins found in file: " + path);

    std::sort(bins.begin(), bins.end(), [](const Bin& x, const Bin& y) { return x.left < y.left; });
    return bins;
}

static int get_bin(double g, const std::vector<Bin>& bins) {
    int lo = 0, hi = static_cast<int>(bins.size()) - 1;
    while (lo <= hi) {
        int mid = lo + (hi - lo) / 2;
        const Bin& b = bins[mid];
        bool in_bin = (mid == static_cast<int>(bins.size()) - 1)
                          ? (g >= b.left && g <= b.right)
                          : (g >= b.left && g < b.right);
        if (in_bin) return mid;
        if (g < b.left) hi = mid - 1; else lo = mid + 1;
    }
    return -1;
}

static std::vector<int> make_random_blocks(std::size_t n, int nblocks, unsigned int seed) {
    if (static_cast<std::size_t>(nblocks) > n) {
        throw std::runtime_error("Number of blocks exceeds number of individuals");
    }
    std::vector<std::size_t> perm(n);
    for (std::size_t i = 0; i < n; ++i) perm[i] = i;

    std::mt19937 rng(seed);
    std::shuffle(perm.begin(), perm.end(), rng);

    std::vector<int> block_of(n, -1);
    for (std::size_t rank = 0; rank < n; ++rank) {
        int b = static_cast<int>((static_cast<long double>(rank) * nblocks) / n);
        if (b >= nblocks) b = nblocks - 1;
        block_of[perm[rank]] = b;
    }
    return block_of;
}

// ---------------------------------------------------------------------------
// Pair class table
// ---------------------------------------------------------------------------

struct ClassedPartner {
    std::uint32_t j;
    std::uint16_t cls;
};

struct PairClasses {
    std::vector<std::string> names;                       // names[0] is always OTHER_CLASS
    std::vector<std::vector<ClassedPartner>> by_row;      // by_row[i] sorted ascending by j
    std::uint64_t n_read = 0, n_mapped = 0, n_unmapped = 0, n_duplicate = 0;
};

static std::vector<std::string> split_tabs_or_ws(const std::string& line) {
    std::vector<std::string> out;
    std::istringstream iss(line);
    std::string tok;
    while (iss >> tok) out.push_back(tok);
    return out;
}

// TSV with a header naming id_col1, id_col2 and class_col. Any pair not listed
// falls through to OTHER_CLASS. Pairs are stored under row max(i,j) since the
// accumulate loop only ever visits j <= i.
static PairClasses read_pair_classes(const std::string& path, const IdTable& ids,
                                     const std::string& id_col1, const std::string& id_col2,
                                     const std::string& class_col) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Could not open pair-classes file: " + path);

    PairClasses pc;
    pc.names.push_back(OTHER_CLASS);
    pc.by_row.assign(ids.iid.size(), {});

    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("Empty pair-classes file: " + path);
    std::vector<std::string> header = split_tabs_or_ws(line);

    auto col_index = [&](const std::string& want) -> int {
        for (std::size_t c = 0; c < header.size(); ++c) {
            std::string h = header[c];
            if (!h.empty() && h[0] == '#') h = h.substr(1);
            if (h == want) return static_cast<int>(c);
        }
        throw std::runtime_error("pair-classes file " + path + " has no column '" + want + "'");
    };
    const int c1 = col_index(id_col1);
    const int c2 = col_index(id_col2);
    const int cc = col_index(class_col);
    const int cmax = std::max({c1, c2, cc});

    std::map<std::string, std::uint16_t> class_id;
    std::vector<std::pair<std::uint64_t, std::uint64_t>> seen_guard;

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> f = split_tabs_or_ws(line);
        if (static_cast<int>(f.size()) <= cmax) continue;
        pc.n_read++;

        auto it1 = ids.index_of_iid.find(f[c1]);
        auto it2 = ids.index_of_iid.find(f[c2]);
        if (it1 == ids.index_of_iid.end() || it2 == ids.index_of_iid.end()) {
            pc.n_unmapped++;
            continue;
        }
        std::size_t a = it1->second, b = it2->second;
        if (a == b) { pc.n_unmapped++; continue; }
        std::size_t hi = std::max(a, b), lo = std::min(a, b);

        const std::string& name = f[cc];
        auto ins = class_id.emplace(name, static_cast<std::uint16_t>(pc.names.size()));
        if (ins.second) pc.names.push_back(name);
        std::uint16_t cid = ins.first->second;

        bool dup = false;
        for (const auto& p : pc.by_row[hi]) {
            if (p.j == static_cast<std::uint32_t>(lo)) { dup = true; break; }
        }
        if (dup) { pc.n_duplicate++; continue; }

        pc.by_row[hi].push_back({static_cast<std::uint32_t>(lo), cid});
        pc.n_mapped++;
    }

    for (auto& v : pc.by_row) {
        std::sort(v.begin(), v.end(),
                  [](const ClassedPartner& x, const ClassedPartner& y) { return x.j < y.j; });
    }
    return pc;
}

// ---------------------------------------------------------------------------
// Row-range recovery for one shard of `--parallel k n` (ported unchanged)
// ---------------------------------------------------------------------------

static uint64_t cumulative_entries(uint64_t i) { return i * (i + 1) / 2; }

struct ShardRange {
    uint64_t row_start;
    uint64_t row_end;  // exclusive
};

static uint64_t triangle_divide_off_diag(uint64_t target) {
    if (target == 0) return 1;
    double approx = (std::sqrt(static_cast<double>(target)) + 1.0);
    int64_t v = static_cast<int64_t>(approx) + 1;
    if (v < 1) v = 1;
    while (v > 1 && static_cast<uint64_t>(v - 1) * static_cast<uint64_t>(v - 2) >= target) v--;
    while (static_cast<uint64_t>(v) * static_cast<uint64_t>(v - 1) < target) v++;
    return static_cast<uint64_t>(v);
}

static uint64_t plink_parallel_row_start(uint64_t n_ids, int parallel_idx, int n_shards) {
    uint64_t ct_tot = n_ids * (n_ids - 1);
    uint64_t target = (ct_tot * static_cast<uint64_t>(parallel_idx)) / static_cast<uint64_t>(n_shards);
    uint64_t v = triangle_divide_off_diag(target);
    return v == 1 ? 0 : v;
}

static bool try_row_start(uint64_t candidate_start, uint64_t n_ids, uint64_t shard_floats,
                          uint64_t& row_end_out) {
    uint64_t consumed = 0;
    uint64_t i = candidate_start;
    while (consumed < shard_floats && i < n_ids) {
        consumed += (i + 1);
        i++;
    }
    if (consumed == shard_floats) {
        row_end_out = i;
        return true;
    }
    return false;
}

static ShardRange resolve_shard_range(uint64_t n_ids, int k, int n_shards, uint64_t shard_floats) {
    uint64_t guess = plink_parallel_row_start(n_ids, k - 1, n_shards);
    static const int offsets[] = {0, -1, 1, -2, 2, -3, 3};
    for (int off : offsets) {
        int64_t candidate = static_cast<int64_t>(guess) + off;
        if (candidate < 0 || static_cast<uint64_t>(candidate) >= n_ids) continue;
        uint64_t row_end;
        if (try_row_start(static_cast<uint64_t>(candidate), n_ids, shard_floats, row_end)) {
            return {static_cast<uint64_t>(candidate), row_end};
        }
    }
    std::ostringstream oss;
    oss << "Could not determine this shard's row range from its file size.\n"
        << "  n_ids=" << n_ids << " parallel=" << k << "/" << n_shards
        << " shard_floats=" << shard_floats << " (best guess row_start=" << guess << ")";
    throw std::runtime_error(oss.str());
}

// ---------------------------------------------------------------------------
// accumulate
// ---------------------------------------------------------------------------

struct AccumulateArgs {
    std::string grm_id, shard, pheno, bins, out, pair_classes;
    std::string id_col1 = "IID1", id_col2 = "IID2", class_col = "cls";
    int parallel_k = 0, parallel_n = 0, nblocks = 0;
    unsigned int seed = 1;
};

static AccumulateArgs parse_accumulate_args(int argc, char** argv) {
    AccumulateArgs a;
    for (int i = 2; i < argc; ++i) {
        std::string key = argv[i];
        auto need = [&](const std::string& flag) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("Missing value for " + flag);
            return argv[++i];
        };
        if (key == "--grm-id") a.grm_id = need(key);
        else if (key == "--shard") a.shard = need(key);
        else if (key == "--pheno") a.pheno = need(key);
        else if (key == "--bins") a.bins = need(key);
        else if (key == "--out") a.out = need(key);
        else if (key == "--pair-classes") a.pair_classes = need(key);
        else if (key == "--id-col1") a.id_col1 = need(key);
        else if (key == "--id-col2") a.id_col2 = need(key);
        else if (key == "--class-col") a.class_col = need(key);
        else if (key == "--nblocks") a.nblocks = std::stoi(need(key));
        else if (key == "--seed") a.seed = static_cast<unsigned int>(std::stoul(need(key)));
        else if (key == "--parallel") {
            a.parallel_k = std::stoi(need(key));
            a.parallel_n = std::stoi(need(key));
        } else throw std::runtime_error("Unknown argument: " + key);
    }
    if (a.grm_id.empty() || a.shard.empty() || a.pheno.empty() || a.bins.empty() ||
        a.out.empty() || a.parallel_k == 0 || a.parallel_n == 0 || a.nblocks <= 0) {
        throw std::runtime_error(
            "Usage: grm_class_tool accumulate --grm-id FILE --shard FILE --parallel K N "
            "--pheno FILE --bins FILE --nblocks INT [--seed INT] "
            "[--pair-classes FILE [--id-col1 IID1] [--id-col2 IID2] [--class-col cls]] "
            "--out FILE");
    }
    return a;
}

static int run_accumulate(int argc, char** argv) {
    AccumulateArgs args = parse_accumulate_args(argc, argv);

    IdTable ids = read_ids(args.grm_id);
    const uint64_t n_ids = ids.iid.size();
    std::vector<double> y = read_pheno_aligned(args.pheno, ids);
    std::vector<Bin> bins = read_bins(args.bins);
    const int nbins = static_cast<int>(bins.size());
    std::vector<int> block_of = make_random_blocks(n_ids, args.nblocks, args.seed);

    PairClasses pc;
    if (!args.pair_classes.empty()) {
        pc = read_pair_classes(args.pair_classes, ids, args.id_col1, args.id_col2, args.class_col);
        std::cerr << "[INFO] pair classes: " << pc.n_read << " read, " << pc.n_mapped
                  << " mapped, " << pc.n_unmapped << " unmapped, " << pc.n_duplicate
                  << " duplicate\n[INFO] classes:";
        for (const auto& nm : pc.names) std::cerr << " " << nm;
        std::cerr << "\n";
        if (pc.n_mapped == 0) {
            std::cerr << "[WARN] no pairs mapped -- check that " << args.id_col1 << "/"
                      << args.id_col2 << " hold IIDs matching the .grm.id\n";
        }
    } else {
        pc.names.push_back(OTHER_CLASS);
        pc.by_row.assign(n_ids, {});
    }
    const int ncls = static_cast<int>(pc.names.size());

    std::ifstream shard(args.shard, std::ios::binary | std::ios::ate);
    if (!shard) throw std::runtime_error("Could not open shard file: " + args.shard);
    std::streamsize byte_size = shard.tellg();
    if (byte_size % 4 != 0) {
        throw std::runtime_error("Shard file size is not a multiple of 4 bytes: " + args.shard);
    }
    uint64_t shard_floats = static_cast<uint64_t>(byte_size) / 4;
    shard.seekg(0);

    ShardRange range = resolve_shard_range(n_ids, args.parallel_k, args.parallel_n, shard_floats);
    std::cerr << "[INFO] Shard " << args.parallel_k << "/" << args.parallel_n
              << " covers rows [" << range.row_start << ", " << range.row_end << ")\n";

    // full[c][k], drop[c][b][k]
    std::vector<std::vector<Acc>> full(ncls, std::vector<Acc>(nbins));
    std::vector<std::vector<std::vector<Acc>>> drop(
        ncls, std::vector<std::vector<Acc>>(args.nblocks, std::vector<Acc>(nbins)));
    std::vector<std::uint64_t> class_hits(ncls, 0);

    float g_f = 0.0f;
    for (uint64_t i = range.row_start; i < range.row_end; ++i) {
        const std::vector<ClassedPartner>& ex = pc.by_row[i];
        std::size_t p = 0;
        for (uint64_t j = 0; j <= i; ++j) {
            shard.read(reinterpret_cast<char*>(&g_f), sizeof(float));
            if (!shard) throw std::runtime_error("Unexpected end of shard while reading row " +
                                                  std::to_string(i));
            if (i == j) continue;

            // one integer compare on rows with no classed partners (the vast majority)
            int cls = 0;
            if (p < ex.size() && ex[p].j == static_cast<std::uint32_t>(j)) {
                cls = ex[p].cls;
                ++p;
            }

            const double yi = y[i], yj = y[j];
            if (std::isnan(yi) || std::isnan(yj)) continue;

            const int k = get_bin(static_cast<double>(g_f), bins);
            if (k < 0) continue;

            const double prod = yi * yj;
            full[cls][k].add(prod);
            class_hits[cls]++;

            const int bi = block_of[i], bj = block_of[j];
            drop[cls][bi][k].add(prod);
            if (bj != bi) drop[cls][bj][k].add(prod);
        }
    }

    std::ofstream out(args.out);
    if (!out) throw std::runtime_error("Could not open output file: " + args.out);
    out << "class\tscope\tblock\tbin_index\tsum\tsum_sq\tn\n";
    out.precision(17);
    for (int c = 0; c < ncls; ++c) {
        for (int k = 0; k < nbins; ++k) {
            if (full[c][k].n == 0) continue;
            out << pc.names[c] << "\tfull\t-1\t" << k << '\t' << full[c][k].sum << '\t'
                << full[c][k].sum_sq << '\t' << full[c][k].n << '\n';
        }
        for (int b = 0; b < args.nblocks; ++b) {
            for (int k = 0; k < nbins; ++k) {
                if (drop[c][b][k].n == 0) continue;
                out << pc.names[c] << "\tdrop\t" << b << '\t' << k << '\t' << drop[c][b][k].sum
                    << '\t' << drop[c][b][k].sum_sq << '\t' << drop[c][b][k].n << '\n';
            }
        }
    }

    std::cerr << "[INFO] binned pairs by class:";
    for (int c = 0; c < ncls; ++c) std::cerr << " " << pc.names[c] << "=" << class_hits[c];
    std::cerr << "\n[INFO] Wrote accumulator to " << args.out << "\n";
    return 0;
}

// ---------------------------------------------------------------------------
// merge
// ---------------------------------------------------------------------------

static double bin_midpoint(const Bin& b) { return 0.5 * (b.left + b.right); }

static double safe_mean(double sum, std::uint64_t n) {
    return n == 0 ? std::numeric_limits<double>::quiet_NaN() : sum / static_cast<double>(n);
}

static double safe_sample_variance(double sum, double sum_sq, std::uint64_t n) {
    if (n < 2) return std::numeric_limits<double>::quiet_NaN();
    const double nd = static_cast<double>(n);
    double var = (sum_sq - (sum * sum) / nd) / static_cast<double>(n - 1);
    if (var < 0.0 && std::abs(var) < 1e-12) var = 0.0;
    return var;
}

static double safe_sample_sd(double sum, double sum_sq, std::uint64_t n) {
    double var = safe_sample_variance(sum, sum_sq, n);
    return std::isnan(var) ? var : std::sqrt(var);
}

static double safe_se_of_mean(double sum, double sum_sq, std::uint64_t n) {
    double var = safe_sample_variance(sum, sum_sq, n);
    return std::isnan(var) ? var : std::sqrt(var / static_cast<double>(n));
}

struct MergeArgs {
    std::string acc_list, bins, out_prefix;
    int nblocks = 0;
};

static MergeArgs parse_merge_args(int argc, char** argv) {
    MergeArgs a;
    for (int i = 2; i < argc; ++i) {
        std::string key = argv[i];
        auto need = [&](const std::string& flag) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("Missing value for " + flag);
            return argv[++i];
        };
        if (key == "--acc-list") a.acc_list = need(key);
        else if (key == "--bins") a.bins = need(key);
        else if (key == "--out-prefix") a.out_prefix = need(key);
        else if (key == "--nblocks") a.nblocks = std::stoi(need(key));
        else throw std::runtime_error("Unknown argument: " + key);
    }
    if (a.acc_list.empty() || a.bins.empty() || a.out_prefix.empty() || a.nblocks <= 0) {
        throw std::runtime_error(
            "Usage: grm_class_tool merge --acc-list FILE --bins FILE --nblocks INT "
            "--out-prefix PREFIX");
    }
    return a;
}

static void jackknife_summary(int k, const std::vector<Acc>& full,
                              const std::vector<std::vector<Acc>>& drop, int nblocks,
                              double& jk_mean, double& jk_var, double& jk_se) {
    std::vector<double> means;
    means.reserve(nblocks);
    for (int b = 0; b < nblocks; ++b) {
        const double jk_sum = full[k].sum - drop[b][k].sum;
        const uint64_t jk_n = full[k].n - drop[b][k].n;
        if (jk_n == 0) {
            jk_mean = jk_var = jk_se = std::numeric_limits<double>::quiet_NaN();
            return;
        }
        means.push_back(jk_sum / static_cast<double>(jk_n));
    }
    double mean_theta = 0.0;
    for (double x : means) mean_theta += x;
    mean_theta /= static_cast<double>(nblocks);

    double ssd = 0.0;
    for (double x : means) {
        const double d = x - mean_theta;
        ssd += d * d;
    }
    jk_mean = mean_theta;
    jk_var = ((static_cast<double>(nblocks) - 1.0) / static_cast<double>(nblocks)) * ssd;
    jk_se = std::sqrt(jk_var);
}

static int run_merge(int argc, char** argv) {
    MergeArgs args = parse_merge_args(argc, argv);
    std::vector<Bin> bins = read_bins(args.bins);
    const int nbins = static_cast<int>(bins.size());

    std::map<std::string, std::vector<Acc>> full;
    std::map<std::string, std::vector<std::vector<Acc>>> drop;

    auto ensure_class = [&](const std::string& c) {
        if (!full.count(c)) {
            full[c] = std::vector<Acc>(nbins);
            drop[c] = std::vector<std::vector<Acc>>(args.nblocks, std::vector<Acc>(nbins));
        }
    };

    std::ifstream list(args.acc_list);
    if (!list) throw std::runtime_error("Could not open shard list file: " + args.acc_list);

    std::string acc_path;
    int n_read = 0;
    while (std::getline(list, acc_path)) {
        if (acc_path.empty()) continue;
        std::ifstream in(acc_path);
        if (!in) throw std::runtime_error("Could not open accumulator file: " + acc_path);

        std::string line;
        std::getline(in, line);  // header
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            std::istringstream iss(line);
            std::string cls, scope;
            int block, bin_index;
            double sum, sum_sq;
            uint64_t n;
            if (!(iss >> cls >> scope >> block >> bin_index >> sum >> sum_sq >> n)) continue;

            if (bin_index < 0 || bin_index >= nbins) {
                throw std::runtime_error("Bin index out of range in " + acc_path +
                                         " -- bins file mismatch between accumulate and merge?");
            }
            ensure_class(cls);
            Acc piece{sum, sum_sq, n};
            if (scope == "full") {
                full[cls][bin_index].add(piece);
            } else if (scope == "drop") {
                if (block < 0 || block >= args.nblocks) {
                    throw std::runtime_error("Block index out of range in " + acc_path +
                                             " -- nblocks mismatch?");
                }
                drop[cls][block][bin_index].add(piece);
            }
        }
        n_read++;
    }
    std::cerr << "[INFO] Merged " << n_read << " accumulator files, " << full.size()
              << " classes\n";

    // "pooled" reproduces what an unclassed run would have produced
    ensure_class("pooled");
    for (const auto& kv : full) {
        if (kv.first == "pooled") continue;
        for (int k = 0; k < nbins; ++k) full["pooled"][k].add(kv.second[k]);
    }
    for (const auto& kv : drop) {
        if (kv.first == "pooled") continue;
        for (int b = 0; b < args.nblocks; ++b)
            for (int k = 0; k < nbins; ++k) drop["pooled"][b][k].add(kv.second[b][k]);
    }

    std::ofstream out_full(args.out_prefix + ".full.tsv");
    std::ofstream out_jk(args.out_prefix + ".jk.tsv");
    out_full << "class\tbin_index\tbin_left\tbin_right\tbin_midpoint\tfull_sum\tfull_sum_sq\t"
                "full_n\tfull_mean\tfull_sd\tfull_se\tjk_mean\tjk_var\tjk_se\n";
    out_jk << "class\tblock\tbin_index\tbin_left\tbin_right\tbin_midpoint\tjk_sum\tjk_sum_sq\t"
              "jk_n\tjk_mean\tjk_sd\tjk_se\n";
    out_full.precision(17);
    out_jk.precision(17);

    for (const auto& kv : full) {
        const std::string& cls = kv.first;
        for (int k = 0; k < nbins; ++k) {
            if (full[cls][k].n == 0) continue;
            double jk_mean = std::numeric_limits<double>::quiet_NaN();
            double jk_var = jk_mean, jk_se = jk_mean;
            jackknife_summary(k, full[cls], drop[cls], args.nblocks, jk_mean, jk_var, jk_se);
            out_full << cls << '\t' << (k + 1) << '\t' << bins[k].left << '\t' << bins[k].right
                     << '\t' << bin_midpoint(bins[k]) << '\t' << full[cls][k].sum << '\t'
                     << full[cls][k].sum_sq << '\t' << full[cls][k].n << '\t'
                     << safe_mean(full[cls][k].sum, full[cls][k].n) << '\t'
                     << safe_sample_sd(full[cls][k].sum, full[cls][k].sum_sq, full[cls][k].n) << '\t'
                     << safe_se_of_mean(full[cls][k].sum, full[cls][k].sum_sq, full[cls][k].n)
                     << '\t' << jk_mean << '\t' << jk_var << '\t' << jk_se << '\n';
        }
        for (int b = 0; b < args.nblocks; ++b) {
            for (int k = 0; k < nbins; ++k) {
                if (full[cls][k].n == 0) continue;
                const double jk_sum = full[cls][k].sum - drop[cls][b][k].sum;
                const double jk_sum_sq = full[cls][k].sum_sq - drop[cls][b][k].sum_sq;
                const uint64_t jk_n = full[cls][k].n - drop[cls][b][k].n;
                out_jk << cls << '\t' << b << '\t' << (k + 1) << '\t' << bins[k].left << '\t'
                       << bins[k].right << '\t' << bin_midpoint(bins[k]) << '\t' << jk_sum << '\t'
                       << jk_sum_sq << '\t' << jk_n << '\t' << safe_mean(jk_sum, jk_n) << '\t'
                       << safe_sample_sd(jk_sum, jk_sum_sq, jk_n) << '\t'
                       << safe_se_of_mean(jk_sum, jk_sum_sq, jk_n) << '\n';
            }
        }
    }

    std::cerr << "[INFO] Wrote " << args.out_prefix << ".full.tsv and " << args.out_prefix
              << ".jk.tsv\n";
    return 0;
}

// ---------------------------------------------------------------------------

static int run_ranges(int argc, char** argv) {
    uint64_t n_ids = 0;
    int k = 0, n_shards = 0;
    for (int i = 2; i < argc; ++i) {
        std::string key = argv[i];
        auto need = [&](const std::string& flag) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("Missing value for " + flag);
            return argv[++i];
        };
        if (key == "--n-ids") n_ids = std::stoull(need(key));
        else if (key == "--parallel") { k = std::stoi(need(key)); n_shards = std::stoi(need(key)); }
        else throw std::runtime_error("Unknown argument: " + key);
    }
    if (n_ids == 0 || k == 0 || n_shards == 0) {
        throw std::runtime_error("Usage: grm_class_tool ranges --n-ids INT --parallel K N");
    }
    uint64_t row_start = plink_parallel_row_start(n_ids, k - 1, n_shards);
    uint64_t row_end = (k == n_shards) ? n_ids : plink_parallel_row_start(n_ids, k, n_shards);
    uint64_t n_floats = cumulative_entries(row_end) - cumulative_entries(row_start);
    std::cout << "shard " << k << "/" << n_shards << ": rows [" << row_start << ", " << row_end
              << "), " << n_floats << " floats (" << n_floats * 4 << " bytes)\n";
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: grm_class_tool <accumulate|merge|ranges> [options]\n";
        return 1;
    }
    std::string cmd = argv[1];
    try {
        if (cmd == "accumulate") return run_accumulate(argc, argv);
        if (cmd == "merge") return run_merge(argc, argv);
        if (cmd == "ranges") return run_ranges(argc, argv);
        std::cerr << "Unknown subcommand: " << cmd << "\n";
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}
