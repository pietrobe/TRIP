#include "domain_decomposition_3D.hpp"

static_assert(__cplusplus >= 201703L, "domain_decomposition_3D requires C++17 or later");

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <stdexcept>

namespace domain_decomposition_3D {

  // ---------------------------------------------------------------------------
  // prime_factors
  // ---------------------------------------------------------------------------

  std::vector<int>
  prime_factors(int n) {
    std::vector<int> factors;

    for (int i = 2; static_cast<long long>(i) * i <= n; ++i) {
      while (n % i == 0) {
        factors.push_back(i);
        n /= i;
      }
    }

    if (n > 1) {
      factors.push_back(n);
    }

    return factors;  // NRVO / move
  }

  // ---------------------------------------------------------------------------
  // DomainDecomposition class definition
  // ---------------------------------------------------------------------------
  void
  get_1d_decomposition(int n_global, int p_dims, int p_coord, int &n_local, int &start) noexcept {
    const int base      = n_global / p_dims;
    const int remainder = n_global % p_dims;

    if (p_coord < remainder) {
      n_local = base + 1;
      start   = p_coord * (base + 1);
    } else {
      n_local = base;
      start   = remainder * (base + 1) + (p_coord - remainder) * base;
    }
  }

  // ---------------------------------------------------------------------------
  // Helpers: read an environment variable as a positive int, or fall back to a
  // default value.
  // ---------------------------------------------------------------------------
  static int
  read_env_int_or_default(const char *name, int default_value) noexcept {
    const char *value = std::getenv(name);
    if (value == nullptr || *value == '\0') {
      return default_value;
    }

    char *end    = nullptr;
    long  parsed = std::strtol(value, &end, 10);
    if (end == value || *end != '\0' || parsed <= 0 || parsed > INT_MAX) {
      return default_value;
    }

    return static_cast<int>(parsed);
  }

  // ---------------------------------------------------------------------------
  // DomainInfo member functions
  // ---------------------------------------------------------------------------
  int
  DomainInfo::rank_of(int px, int py, int pz) const noexcept {
    return px * (Py * Pz) + py * Pz + pz;
  }

  int
  DomainInfo::neighbor_rank(int dix, int diy, int diz) const noexcept {
    const int nx = this->ix + dix;
    const int ny = this->iy + diy;
    const int nz = this->iz + diz;

    if (nx < 0 || nx >= Px)
      return MPI_PROC_NULL;
    if (ny < 0 || ny >= Py)
      return MPI_PROC_NULL;
    if (nz < 0 || nz >= Pz)
      return MPI_PROC_NULL;

    return this->rank_of(nx, ny, nz);
  }

  // ---------------------------------------------------------------------------
  // DomainDecomposition — private static helpers
  // ---------------------------------------------------------------------------
  void
  DomainDecomposition::factorize_procs(int size, int Nx, int Ny, int Nz, int &Px, int &Py, int &Pz) {
    // Enumerate divisors of size once (O(√size)) so the double loop visits
    // only valid candidates instead of scanning 1..size × 1..size.
    std::vector<int> divs;
    divs.reserve(64);
    for (int i = 1; static_cast<long long>(i) * i <= size; ++i) {
      if (size % i == 0) {
        divs.push_back(i);
        if (i != size / i)
          divs.push_back(size / i);
      }
    }
    std::sort(divs.begin(), divs.end());

    int    best_px = 1, best_py = 1, best_pz = size;
    double best_score = 1e18;

    for (const int px : divs) {
      if (px > Nx)
        break;  // divs is sorted; no subsequent px can be valid

      const int rem = size / px;

      for (const int py : divs) {
        if (py > Ny || py > rem)
          break;  // divs is sorted; no subsequent py can be valid
        if (rem % py != 0)
          continue;

        const int pz = rem / py;
        if (pz > Nz)
          continue;

        // Local subdomain extents for this factorization
        const double lx = static_cast<double>(Nx) / px;
        const double ly = static_cast<double>(Ny) / py;
        const double lz = static_cast<double>(Nz) / pz;

        // Score: sum of pairwise aspect ratios (minimum = 6.0 for a perfect cube)
        const double score = (lx / ly + ly / lx) + (ly / lz + lz / ly) + (lx / lz + lz / lx);

        if (score <= best_score) {
          best_score = score;
          best_px    = px;
          best_py    = py;
          best_pz    = pz;

          DomainInfo info{};
          info.Px     = px;
          info.Py     = py;
          info.Pz     = pz;
          info.score_ = score;
          info.N_x    = Nx;
          info.N_y    = Ny;
          info.N_z    = Nz;
          this->factorization_cache_.emplace_back(std::move(info));
        }

        if (this->factorization_cache_.size() > 15)
          this->factorization_cache_.pop_front();
      }
    }

    Px = best_px;
    Py = best_py;
    Pz = best_pz;
  }

  DomainInfo
  DomainDecomposition::info_for_rank(int rank) const {
    if (rank < 0 || rank >= mpi_size_) {
      throw std::out_of_range("DomainDecomposition::info_for_rank rank out of range");
    }

    DomainInfo result{};

    result.set_global_decomposition(info_.N_x, info_.N_y, info_.N_z, info_.Px, info_.Py, info_.Pz);
    // result.Px = info_.Px;
    // result.Py = info_.Py;
    // result.Pz = info_.Pz;

    // result.N_x = info_.N_x;
    // result.N_y = info_.N_y;
    // result.N_z = info_.N_z;

    // Map linear rank -> 3D processor coordinate (row-major)
    // rank = ix*(Py*Pz) + iy*Pz + iz
    result.set_rank(rank);
    // result.ix = rank / (result.Py * result.Pz);
    // result.iy = (rank % (result.Py * result.Pz)) / result.Pz;
    // result.iz = rank % result.Pz;

    // get_1d_decomposition(result.N_x, result.Px, result.ix, result.N_local_x, result.local_start_x);
    // get_1d_decomposition(result.N_y, result.Py, result.iy, result.N_local_y, result.local_start_y);
    // get_1d_decomposition(result.N_z, result.Pz, result.iz, result.N_local_z, result.local_start_z);

    return result;
  }

  std::deque<DomainInfo>
  DomainDecomposition::best_infos(int rank) const {

    std::deque<DomainInfo> up_infos;

    for (auto info_itr = this->factorization_cache_.begin(); info_itr != this->factorization_cache_.end(); ++info_itr) {

      DomainInfo new_info;
      new_info = *info_itr;  // copy the cached info

      // Map linear rank -> 3D processor coordinate (row-major)
      // rank = ix*(Py*Pz) + iy*Pz + iz
      new_info.set_rank(rank);

      get_1d_decomposition(info_itr->N_x, info_itr->Px, new_info.ix, new_info.N_local_x, new_info.local_start_x);
      get_1d_decomposition(info_itr->N_y, info_itr->Py, new_info.iy, new_info.N_local_y, new_info.local_start_y);
      get_1d_decomposition(info_itr->N_z, info_itr->Pz, new_info.iz, new_info.N_local_z, new_info.local_start_z);

      up_infos.push_back(std::move(new_info));
    }
    return up_infos;
  }

  // ---------------------------------------------------------------------------
  // DomainDecomposition — constructor
  // ---------------------------------------------------------------------------

  DomainDecomposition::DomainDecomposition(int mpi_rank, int mpi_size, int N_x, int N_y, int N_z)
      : mpi_rank_(mpi_rank), mpi_size_(mpi_size) {
    // 1. Initialise global extents
    info_.N_x = N_x;
    info_.N_y = N_y;
    info_.N_z = N_z;

    // 2. Find the best (Px, Py, Pz) factorisation of mpi_size
    factorize_procs(mpi_size, N_x, N_y, N_z, info_.Px, info_.Py, info_.Pz);

    if (mpi_rank == 0) {
      std::printf("[DomainDecomposition] grid: %d x %d x %d procs over %d x %d x %d domain\n", info_.Px, info_.Py,
                  info_.Pz, N_x, N_y, N_z);
    }

    // 3. Compute local info for this rank via the generic helper.
    info_ = info_for_rank(mpi_rank);
  }

  // ---------------------------------------------------------------------------
  // DomainInfo — global load imbalance
  // ---------------------------------------------------------------------------

  std::tuple<double, double, double, double>
  DomainInfo::global_load_imbalance() const {
    const int total_ranks = Px * Py * Pz;
    int       max_cells   = 0;
    int       min_cells   = INT_MAX;

    for (int rank = 0; rank < total_ranks; ++rank) {
      DomainInfo info_rank{};
      info_rank.set_global_decomposition(this->N_x, this->N_y, this->N_z, this->Px, this->Py, this->Pz);
      info_rank.set_rank(rank);
      const int num_cells = info_rank.number_of_cells();
      max_cells           = std::max(max_cells, num_cells);
      min_cells           = std::min(min_cells, num_cells);
    }

    const double mean             = static_cast<double>(N_x * N_y * N_z) / static_cast<double>(total_ranks);
    const double max_over_mean    = static_cast<double>(max_cells) / mean;
    const double max_over_min     = static_cast<double>(max_cells) / static_cast<double>(min_cells);
    const double mean_over_min    = mean / static_cast<double>(min_cells);
    const double spread_over_mean = (mean + std::max(std::abs(static_cast<double>(min_cells) - mean),     //
                                                     std::abs(static_cast<double>(max_cells) - mean))) /  //
                                    mean;                                                                 //
    return {max_over_mean, max_over_min, mean_over_min, spread_over_mean};
  }

  // ---------------------------------------------------------------------------
  // DomainDecomposition — analysis and demo
  // ---------------------------------------------------------------------------
  bool
  analyze_decomposition(const DomainInfo &info, int mpi_size, bool verbose) {
    if (mpi_size <= 0) {
      std::fprintf(stderr, "Invalid mpi_size: %d\n", mpi_size);
      return false;
    }

    if (info.Px * info.Py * info.Pz != mpi_size) {
      std::fprintf(stderr, "Invalid decomposition: Px*Py*Pz=%d does not match mpi_size=%d\n", info.Px * info.Py * info.Pz,
                   mpi_size);
      return false;
    }

    int num_cells_max = 0;
    int num_cells_min = INT_MAX;
    int average_cells = 0;
    int total_cells   = 0;

    std::vector<DomainInfo> all_info;
    all_info.reserve(mpi_size);

    for (int rank_i = 0; rank_i < mpi_size; ++rank_i) {
      DomainInfo info_i{};
      info_i.Px  = info.Px;
      info_i.Py  = info.Py;
      info_i.Pz  = info.Pz;
      info_i.N_x = info.N_x;
      info_i.N_y = info.N_y;
      info_i.N_z = info.N_z;

      info_i.set_rank(rank_i);

      // info_i.ix = rank_i / (info_i.Py * info_i.Pz);
      // info_i.iy = (rank_i % (info_i.Py * info_i.Pz)) / info_i.Pz;
      // info_i.iz = rank_i % info_i.Pz;

      const auto fill_1d = [](int n_global, int p_dims, int p_coord, int &n_local, int &start) {
        const int base      = n_global / p_dims;
        const int remainder = n_global % p_dims;

        if (p_coord < remainder) {
          n_local = base + 1;
          start   = p_coord * (base + 1);
        } else {
          n_local = base;
          start   = remainder * (base + 1) + (p_coord - remainder) * base;
        }
      };

      fill_1d(info_i.N_x, info_i.Px, info_i.ix, info_i.N_local_x, info_i.local_start_x);
      fill_1d(info_i.N_y, info_i.Py, info_i.iy, info_i.N_local_y, info_i.local_start_y);
      fill_1d(info_i.N_z, info_i.Pz, info_i.iz, info_i.N_local_z, info_i.local_start_z);

      const int num_cells = info_i.number_of_cells();
      num_cells_max       = std::max(num_cells_max, num_cells);
      num_cells_min       = std::min(num_cells_min, num_cells);
      average_cells += num_cells;
      total_cells += num_cells;
      all_info.push_back(info_i);
    }

    average_cells /= mpi_size;

    const auto [imb_max_mean, imb_max_min, imb_mean_min, imb_spread] = info.global_load_imbalance();

    const int  total_local_domains = static_cast<int>(all_info.size());
    const bool domain_count_ok     = (total_local_domains == mpi_size);

    const long long expected_total_cells = 1LL * info.N_x * info.N_y * info.N_z;
    const bool      cell_sum_ok          = (static_cast<long long>(total_cells) == expected_total_cells);

    int overlap_count = 0;
    for (int i = 0; i < mpi_size; ++i) {
      const DomainInfo &a = all_info[i];
      for (int j = i + 1; j < mpi_size; ++j) {
        const DomainInfo &b = all_info[j];

        const bool ox =
            (a.local_start_x < b.local_start_x + b.N_local_x) && (b.local_start_x < a.local_start_x + a.N_local_x);
        const bool oy =
            (a.local_start_y < b.local_start_y + b.N_local_y) && (b.local_start_y < a.local_start_y + a.N_local_y);
        const bool oz =
            (a.local_start_z < b.local_start_z + b.N_local_z) && (b.local_start_z < a.local_start_z + a.N_local_z);

        if (ox && oy && oz) {
          if (verbose)
            std::printf(
                "  OVERLAP: rank %d X[%d+%d] Y[%d+%d] Z[%d+%d]"
                " <-> rank %d X[%d+%d] Y[%d+%d] Z[%d+%d]\n",
                i, a.local_start_x, a.N_local_x, a.local_start_y, a.N_local_y, a.local_start_z, a.N_local_z, j,
                b.local_start_x, b.N_local_x, b.local_start_y, b.N_local_y, b.local_start_z, b.N_local_z);
          ++overlap_count;
        }
      }
    }

    const bool ok = domain_count_ok && cell_sum_ok && (overlap_count == 0);

    if (verbose) {
      std::printf("  Cells: min = %d, max = %d, average = %d, total = %d\n", num_cells_min, num_cells_max, average_cells,
                  total_cells);
      std::printf("  Load imbalance (max/mean):              %.3f\n", imb_max_mean);
      std::printf("  Load imbalance (max/min):               %.3f\n", imb_max_min);
      std::printf("  Load imbalance (mean/min):              %.3f\n", imb_mean_min);
      std::printf("  Load imbalance (spread/mean):           %.3f\n", imb_spread);

      if (domain_count_ok)
        std::printf("  Local-domain count check: OK - %d local domains for %d MPI ranks.\n", total_local_domains,
                    mpi_size);
      else
        std::printf("  Local-domain count check: FAILED - %d local domains for %d MPI ranks.\n", total_local_domains,
                    mpi_size);

      if (cell_sum_ok)
        std::printf("  Cell-sum check: OK - sum(local cells) = %d matches global cells = %lld.\n", total_cells,
                    expected_total_cells);
      else
        std::printf("  Cell-sum check: FAILED - sum(local cells) = %d, expected global cells = %lld.\n", total_cells,
                    expected_total_cells);

      if (overlap_count == 0)
        std::printf("  Overlap check: OK - no overlapping subdomains.\n");
      else
        std::printf("  Overlap check: FAILED - %d overlapping pair(s) found.\n", overlap_count);

      FILE *csv_file = std::fopen("decomposition_info.csv", "w");
      if (csv_file != nullptr) {
        std::fprintf(csv_file,
                     "rank,ix,iy,iz,N_local_x,N_local_y,N_local_z,local_start_x,local_start_y,local_start_z,"
                     "num_cells,N_x,N_y,N_z,Px,Py,Pz\n");
        for (int i = 0; i < static_cast<int>(all_info.size()); ++i) {
          const DomainInfo &d = all_info[i];
          std::fprintf(csv_file, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", i, d.ix, d.iy, d.iz, d.N_local_x,
                       d.N_local_y, d.N_local_z, d.local_start_x, d.local_start_y, d.local_start_z, d.number_of_cells(),
                       d.N_x, d.N_y, d.N_z, d.Px, d.Py, d.Pz);
        }
        std::fclose(csv_file);
        std::printf("  CSV file 'decomposition_info.csv' written.\n");
      } else {
        std::fprintf(stderr, "  Failed to open 'decomposition_info.csv' for writing.\n");
      }
    }

    return ok;
  }

  // ---------------------------------------------------------------------------
  // select_best_decomposition
  // ---------------------------------------------------------------------------

  // Load-imbalance model used for ranking decompositions.
  // Pick one index into the tuple returned by global_load_imbalance():
  //   0 → max / mean   (default)
  //   1 → max / min
  //   2 → mean / min
  //   3 → (mean + max(|min-mean|, |max-mean|)) / mean
  DecompositionSelection
  select_best_decomposition(const std::deque<DomainInfo> &candidates, bool verbose, const double li_rel_tol,
                            LI_RankingModel li_model) {
    DecompositionSelection result{};

    const auto pick_li = [li_model](const std::tuple<double, double, double, double> &t) {
      switch (li_model) {
        case LI_RankingModel::MAX_OVER_MEAN:
          return std::get<0>(t);
        case LI_RankingModel::MAX_OVER_MIN:
          return std::get<1>(t);
        case LI_RankingModel::MEAN_OVER_MIN:
          return std::get<2>(t);
        case LI_RankingModel::SPREAD_OVER_MEAN:
          return std::get<3>(t);
      }
      return std::get<0>(t);
    };

    for (const auto &info : candidates) {
      const double aspect_ratio   = info.aspect_ratio_score();
      const int    p_max          = std::max(info.Px, std::max(info.Py, info.Pz));
      const int    p_min          = std::min(info.Px, std::min(info.Py, info.Pz));
      const double grid_ratio     = static_cast<double>(p_max) / p_min;
      const double load_imbalance = pick_li(info.global_load_imbalance());

      const double result_load_imbalance = result.load_imbalance + result.load_imbalance * li_rel_tol;
      const bool   is_better =
          !result.found || std::tie(load_imbalance, aspect_ratio, info.score_, grid_ratio) <
                               std::tie(result_load_imbalance, result.aspect_ratio, result.score, result.grid_ratio);

      if (verbose) {
        std::printf(
            "C: aspect ratio: %.3f || score: %.3f || g-ratio max/min: %.3f || imbalance: %.3f || grid: %d x %d x %d%s\n",
            aspect_ratio, info.score_, grid_ratio, load_imbalance, info.Px, info.Py, info.Pz, is_better ? "  *" : "");
      }

      if (is_better) {
        result.best_info      = info;
        result.aspect_ratio   = aspect_ratio;
        result.score          = info.score_;
        result.grid_ratio     = grid_ratio;
        result.load_imbalance = load_imbalance;
        result.found          = true;
      }
    }

    return result;
  }

  // ---------------------------------------------------------------------------
  // demo_decompose_domain_3d
  // ---------------------------------------------------------------------------
  void
  demo_decompose_domain_3d() {
    const int default_mpi_rank = 0;
    const int default_mpi_size = 4096 * 6;

    const int mpi_rank = read_env_int_or_default("DD_MPI_RANK", default_mpi_rank);
    const int mpi_size = read_env_int_or_default("DD_MPI_SIZE", default_mpi_size);

    const int N_x = read_env_int_or_default("DD_HDF_N_X", 128);
    const int N_y = read_env_int_or_default("DD_HDF_N_Y", 128);
    const int N_z = read_env_int_or_default("DD_HDF_N_Z", 128);

    std::printf("Environment overrides (defaults in use when unset/invalid):\n");
    std::printf("export DD_MPI_RANK=%d\n", mpi_rank);
    std::printf("export DD_MPI_SIZE=%d\n\n\n", mpi_size);
    std::printf("export DD_HDF_N_X=%d\n", N_x);
    std::printf("export DD_HDF_N_Y=%d\n", N_y);
    std::printf("export DD_HDF_N_Z=%d\n\n\n", N_z);

    auto primes    = prime_factors(mpi_size);
    int  max_prime = primes.empty() ? 1 : std::max(0, primes.back());
    std::printf("MPI size: %d  prime factors:", mpi_size);
    for (int p : primes) {
      std::printf(" %d", p);
    }
    std::printf("  max prime factor: %d\n", max_prime);

    if (max_prime > std::min({N_x, N_y, N_z})) {
      std::fprintf(stderr, "Invalid domain size: N_x=%d, N_y=%d, N_z=%d, for MPI size=%d with max prime factor %d\n", N_x,
                   N_y, N_z, mpi_size, max_prime);
      return;
    }

    std::printf("Decomposing domain of size %d x %d x %d = %d across %d MPI ranks...\n", N_x, N_y, N_z, N_x * N_y * N_z,
                mpi_size);

    DomainDecomposition dd(mpi_rank, mpi_size, N_x, N_y, N_z);
    const DomainInfo   &d = dd.info();

    std::printf("Rank %d/%d: local domain X[%d+%d] Y[%d+%d] Z[%d+%d]  grid: %d x %d x %d\n", mpi_rank, mpi_size,
                d.local_start_x, d.N_local_x, d.local_start_y, d.N_local_y, d.local_start_z, d.N_local_z, d.Px, d.Py,
                d.Pz);

    analyze_decomposition(d, mpi_size);

    auto best_infos = dd.best_infos(mpi_rank);
    auto selection  = select_best_decomposition(best_infos, true, 0.01, LI_RankingModel::MAX_OVER_MIN);

    analyze_decomposition(selection.best_info, mpi_size);

    if (selection.found) {
      std::printf("\n===> Selected best decomposition for rank %d: =============== \n", mpi_rank);
      std::printf(">   grid: %d x %d x %d\n", selection.best_info.Px, selection.best_info.Py, selection.best_info.Pz);
      std::printf(
          ">   aspect ratio: %.2f || score: %.2f || grid ratio max/min: %.2f || load imbalance: "
          "%.2f\n\n",
          selection.aspect_ratio, selection.score, selection.grid_ratio, selection.load_imbalance);
    }
  }

  // ---------------------------------------------------------------------------
  // demo_compare_ranking_models
  // ---------------------------------------------------------------------------
  void
  demo_compare_ranking_models() {
    const int mpi_rank = read_env_int_or_default("DD_MPI_RANK", 0);
    const int mpi_size = read_env_int_or_default("DD_MPI_SIZE", 4096 * 6);
    const int N_x      = read_env_int_or_default("DD_HDF_N_X", 128);
    const int N_y      = read_env_int_or_default("DD_HDF_N_Y", 128);
    const int N_z      = read_env_int_or_default("DD_HDF_N_Z", 128);

    std::printf("=== Ranking Model Comparison ===\n");
    std::printf("Domain: %d x %d x %d | MPI size: %d | Rank: %d\n\n", N_x, N_y, N_z, mpi_size, mpi_rank);

    const auto primes    = prime_factors(mpi_size);
    const int  max_prime = primes.empty() ? 1 : primes.back();

    if (max_prime > std::min({N_x, N_y, N_z})) {
      std::fprintf(stderr, "Invalid domain: N_x=%d, N_y=%d, N_z=%d for MPI size=%d (max prime=%d)\n", N_x, N_y, N_z,
                   mpi_size, max_prime);
      return;
    }

    DomainDecomposition dd(mpi_rank, mpi_size, N_x, N_y, N_z);
    const auto          candidates = dd.best_infos(mpi_rank);

    if (candidates.empty()) {
      std::printf("No candidates found.\n");
      return;
    }

    struct ModelEntry {
      const char     *name;
      LI_RankingModel model;
    };

    constexpr ModelEntry models[] = {
        {"MAX_OVER_MEAN", LI_RankingModel::MAX_OVER_MEAN},
        {"MAX_OVER_MIN", LI_RankingModel::MAX_OVER_MIN},
        {"MEAN_OVER_MIN", LI_RankingModel::MEAN_OVER_MIN},
        {"SPREAD_OVER_MEAN", LI_RankingModel::SPREAD_OVER_MEAN},
    };

    std::printf("\n\n  %-20s  %-13s  %-10s  %-10s  %-10s  %-12s  %-12s  %-8s  %-8s\n", "Model", "Grid", "max/mean",
                "max/min", "mean/min", "spread/mean", "aspect_ratio", "score", "g-ratio");
    for (int i = 0; i < 110; ++i) std::putchar('-');
    std::putchar('\n');

    for (const auto &entry : models) {
      const auto sel = select_best_decomposition(candidates, false, 0.01, entry.model);
      if (!sel.found) {
        std::printf("  %-20s  (no selection)\n", entry.name);
        continue;
      }

      char grid_str[24];
      std::snprintf(grid_str, sizeof(grid_str), "%dx%dx%d", sel.best_info.Px, sel.best_info.Py, sel.best_info.Pz);

      const auto [m0, m1, m2, m3] = sel.best_info.global_load_imbalance();
      const double metrics[4]     = {m0, m1, m2, m3};
      const int    widths[4]      = {10, 10, 10, 12};
      const int    active_idx     = static_cast<int>(entry.model);

      std::printf("  %-20s  %-13s  ", entry.name, grid_str);
      for (int k = 0; k < 4; ++k) {
        if (k == active_idx)
          std::fputs("\033[1;33m", stdout);
        std::printf("%-*.4f", widths[k], metrics[k]);
        if (k == active_idx)
          std::fputs("\033[0m", stdout);
        std::fputs("  ", stdout);
      }
      std::printf("%-12.4f  %-8.3f  %-8.3f\n", sel.aspect_ratio, sel.score, sel.grid_ratio);
    }
    std::putchar('\n');
  }

}  // namespace domain_decomposition_3D
