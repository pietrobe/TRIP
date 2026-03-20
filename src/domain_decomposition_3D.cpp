#include "domain_decomposition_3D.hpp"

static_assert(__cplusplus >= 201703L, "domain_decomposition_3D requires C++17 or later");

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <stdexcept>

namespace domain_decomposition_3D
{

	// ---------------------------------------------------------------------------
	// prime_factors
	// ---------------------------------------------------------------------------

	std::vector<int>
	prime_factors(int n)
	{
		std::vector<int> factors;

		for (int i = 2; static_cast<long long>(i) * i <= n; ++i)
		{
			while (n % i == 0)
			{
				factors.push_back(i);
				n /= i;
			}
		}

		if (n > 1)
		{
			factors.push_back(n);
		}

		return factors; // NRVO / move
	}

	// ---------------------------------------------------------------------------
	// Helpers: read an environment variable as a positive int, or fall back to a
	// default value.
	// ---------------------------------------------------------------------------

	static int
	read_env_int_or_default(const char *name, int default_value) noexcept
	{
		const char *value = std::getenv(name);
		if (value == nullptr || *value == '\0')
		{
			return default_value;
		}

		char *end	 = nullptr;
		long  parsed = std::strtol(value, &end, 10);
		if (end == value || *end != '\0' || parsed <= 0 || parsed > INT_MAX)
		{
			return default_value;
		}

		return static_cast<int>(parsed);
	}

	// ---------------------------------------------------------------------------
	// DomainInfo member functions
	// ---------------------------------------------------------------------------

	int
	DomainInfo::rank_of(int px, int py, int pz) const noexcept
	{
		return px * (Py * Pz) + py * Pz + pz;
	}

	int
	DomainInfo::neighbor_rank(int dix, int diy, int diz) const noexcept
	{
		const int nx = ix + dix;
		const int ny = iy + diy;
		const int nz = iz + diz;

		if (nx < 0 || nx >= Px) return MPI_PROC_NULL;
		if (ny < 0 || ny >= Py) return MPI_PROC_NULL;
		if (nz < 0 || nz >= Pz) return MPI_PROC_NULL;

		return rank_of(nx, ny, nz);
	}

	// ---------------------------------------------------------------------------
	// DomainDecomposition — private static helpers
	// ---------------------------------------------------------------------------

	void
	DomainDecomposition::factorize_procs(int size, int Nx, int Ny, int Nz, int &Px, int &Py, int &Pz)
	{
		int	   best_px = 1, best_py = 1, best_pz = size;
		double best_score = 1e18;

		for (int px = 1; px <= size; ++px)
		{
			if (size % px != 0) continue;

			for (int py = 1; py <= size / px; ++py)
			{
				if ((size / px) % py != 0) continue;

				const int pz = size / (px * py);

				// Local subdomain extents for this factorization
				const double lx = static_cast<double>(Nx) / px;
				const double ly = static_cast<double>(Ny) / py;
				const double lz = static_cast<double>(Nz) / pz;

				// Score: sum of pairwise aspect ratios (minimum = 6.0 for a perfect cube)
				const double score = (lx / ly + ly / lx) + (ly / lz + lz / ly) + (lx / lz + lz / lx);

				// std::printf("  (%d x %d x %d)  local=(%.1f x %.1f x %.1f)  score=%.3f%s\n",
				//             px, py, pz, lx, ly, lz, score,
				//             score < best_score ? "  <-- new best" : "");

				if (score <= best_score)
				{
					best_score = score;
					best_px	   = px;
					best_py	   = py;
					best_pz	   = pz;

					DomainInfo info{};
					info.Px		= px;
					info.Py		= py;
					info.Pz		= pz;
					info.score_ = score;
					info.N_x	= Nx;
					info.N_y	= Ny;
					info.N_z	= Nz;
					this->factorization_cache_.emplace_back(std::move(info));
				}

				if (this->factorization_cache_.size() > 10)
				{
					// Assuming factorization_cache_ is sorted by score_ ascending,
					// the worst score is at the back.
					this->factorization_cache_.pop_front();
				}
			}
		}

		Px = best_px;
		Py = best_py;
		Pz = best_pz;
	}

	void
	DomainDecomposition::get_1d_decomposition(int n_global, int p_dims, int p_coord, int &n_local, int &start) noexcept
	{
		const int base		= n_global / p_dims;
		const int remainder = n_global % p_dims;

		if (p_coord < remainder)
		{
			n_local = base + 1;
			start	= p_coord * (base + 1);
		}
		else
		{
			n_local = base;
			start	= remainder * (base + 1) + (p_coord - remainder) * base;
		}
	}

	DomainInfo
	DomainDecomposition::info_for_rank(int rank) const
	{
		if (rank < 0 || rank >= mpi_size_)
		{
			throw std::out_of_range("DomainDecomposition::info_for_rank rank out of range");
		}

		DomainInfo result{};
		result.Px = info_.Px;
		result.Py = info_.Py;
		result.Pz = info_.Pz;

		result.N_x = info_.N_x;
		result.N_y = info_.N_y;
		result.N_z = info_.N_z;

		// Map linear rank -> 3D processor coordinate (row-major)
		// rank = ix*(Py*Pz) + iy*Pz + iz
		result.ix = rank / (result.Py * result.Pz);
		result.iy = (rank % (result.Py * result.Pz)) / result.Pz;
		result.iz = rank % result.Pz;

		get_1d_decomposition(result.N_x, result.Px, result.ix, result.N_local_x, result.local_start_x);
		get_1d_decomposition(result.N_y, result.Py, result.iy, result.N_local_y, result.local_start_y);
		get_1d_decomposition(result.N_z, result.Pz, result.iz, result.N_local_z, result.local_start_z);

		return result;
	}

	std::deque<DomainInfo>
	DomainDecomposition::best_infos(int rank) const
	{
		std::deque<DomainInfo> up_infos;

		for (auto info_itr = this->factorization_cache_.begin(); info_itr != this->factorization_cache_.end(); ++info_itr)
		{
			DomainInfo new_info;
			new_info = *info_itr; // copy the cached info

			// Map linear rank -> 3D processor coordinate (row-major)
			// rank = ix*(Py*Pz) + iy*Pz + iz
			new_info.ix = rank / (info_itr->Py * info_itr->Pz);
			new_info.iy = (rank % (info_itr->Py * info_itr->Pz)) / info_itr->Pz;
			new_info.iz = rank % info_itr->Pz;

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
		: mpi_rank_(mpi_rank), mpi_size_(mpi_size)
	{
		// 1. Initialise global extents
		info_.N_x = N_x;
		info_.N_y = N_y;
		info_.N_z = N_z;

		// 2. Find the best (Px, Py, Pz) factorisation of mpi_size
		factorize_procs(mpi_size, N_x, N_y, N_z, info_.Px, info_.Py, info_.Pz);

		if (mpi_rank == 0)
		{
			std::printf("[DomainDecomposition] grid: %d x %d x %d procs over %d x %d x %d domain\n", info_.Px, info_.Py,
						info_.Pz, N_x, N_y, N_z);
		}

		// 3. Compute local info for this rank via the generic helper.
		info_ = info_for_rank(mpi_rank);
	}

	void
	analyze_decomposition(const DomainInfo &info, int mpi_size)
	{
		if (mpi_size <= 0)
		{
			std::fprintf(stderr, "Invalid mpi_size: %d\n", mpi_size);
			return;
		}

		if (info.Px * info.Py * info.Pz != mpi_size)
		{
			std::fprintf(stderr, "Invalid decomposition: Px*Py*Pz=%d does not match mpi_size=%d\n",
						 info.Px * info.Py * info.Pz, mpi_size);
			return;
		}

		int num_cells_max = 0;
		int num_cells_min = INT_MAX;
		int average_cells = 0;
		int total_cells	  = 0;

		std::vector<DomainInfo> all_info;
		all_info.reserve(mpi_size);

		for (int rank_i = 0; rank_i < mpi_size; ++rank_i)
		{
			DomainInfo info_i{};
			info_i.Px  = info.Px;
			info_i.Py  = info.Py;
			info_i.Pz  = info.Pz;
			info_i.N_x = info.N_x;
			info_i.N_y = info.N_y;
			info_i.N_z = info.N_z;

			info_i.ix = rank_i / (info_i.Py * info_i.Pz);
			info_i.iy = (rank_i % (info_i.Py * info_i.Pz)) / info_i.Pz;
			info_i.iz = rank_i % info_i.Pz;

			const auto fill_1d = [](int n_global, int p_dims, int p_coord, int &n_local, int &start)
			{
				const int base		= n_global / p_dims;
				const int remainder = n_global % p_dims;

				if (p_coord < remainder)
				{
					n_local = base + 1;
					start	= p_coord * (base + 1);
				}
				else
				{
					n_local = base;
					start	= remainder * (base + 1) + (p_coord - remainder) * base;
				}
			};

			fill_1d(info_i.N_x, info_i.Px, info_i.ix, info_i.N_local_x, info_i.local_start_x);
			fill_1d(info_i.N_y, info_i.Py, info_i.iy, info_i.N_local_y, info_i.local_start_y);
			fill_1d(info_i.N_z, info_i.Pz, info_i.iz, info_i.N_local_z, info_i.local_start_z);

			const int num_cells = info_i.number_of_cells();
			num_cells_max		= std::max(num_cells_max, num_cells);
			num_cells_min		= std::min(num_cells_min, num_cells);
			average_cells += num_cells;
			total_cells += num_cells;
			all_info.push_back(info_i);
		}

		average_cells /= mpi_size;
		std::printf("  Cells: min = %d, max = %d, average = %d, total = %d\n", num_cells_min, num_cells_max,
					average_cells, total_cells);

		const double imbalance = static_cast<double>(num_cells_max) / average_cells;
		std::printf("  Load imbalance (max/average): %.3f\n", imbalance);

		int overlap_count = 0;
		for (int i = 0; i < mpi_size; ++i)
		{
			const DomainInfo &a = all_info[i];
			for (int j = i + 1; j < mpi_size; ++j)
			{
				const DomainInfo &b = all_info[j];

				const bool ox = (a.local_start_x < b.local_start_x + b.N_local_x)
								&& (b.local_start_x < a.local_start_x + a.N_local_x);
				const bool oy = (a.local_start_y < b.local_start_y + b.N_local_y)
								&& (b.local_start_y < a.local_start_y + a.N_local_y);
				const bool oz = (a.local_start_z < b.local_start_z + b.N_local_z)
								&& (b.local_start_z < a.local_start_z + a.N_local_z);

				if (ox && oy && oz)
				{
					std::printf(
						"  OVERLAP: rank %d X[%d+%d] Y[%d+%d] Z[%d+%d]"
						" <-> rank %d X[%d+%d] Y[%d+%d] Z[%d+%d]\n",
						i, a.local_start_x, a.N_local_x, a.local_start_y, a.N_local_y, a.local_start_z, a.N_local_z, j,
						b.local_start_x, b.N_local_x, b.local_start_y, b.N_local_y, b.local_start_z, b.N_local_z);
					++overlap_count;
				}
			}
		}

		if (overlap_count == 0)
			std::printf("  Overlap check: OK - no overlapping subdomains.\n");
		else
			std::printf("  Overlap check: FAILED - %d overlapping pair(s) found.\n", overlap_count);
	}

	// ---------------------------------------------------------------------------
	// demo_decompose_domain_3d
	// ---------------------------------------------------------------------------

	void
	demo_decompose_domain_3d()
	{
		constexpr int mpi_rank = 0;
		constexpr int mpi_size = 4096 * 6;

		auto primes	   = prime_factors(mpi_size);
		int	 max_prime = primes.empty() ? 1 : std::max(0, primes.back());
		std::printf("MPI size: %d  prime factors:", mpi_size);
		for (int p : primes)
		{
			std::printf(" %d", p);
		}
		std::printf("  max prime factor: %d\n", max_prime);

		const int N_x = read_env_int_or_default("HDF_N_X", 128);
		const int N_y = read_env_int_or_default("HDF_N_Y", 128);
		const int N_z = read_env_int_or_default("HDF_N_Z", 128);

		if (max_prime > std::min({N_x, N_y, N_z}))
		{
			std::fprintf(stderr,
						 "Invalid domain size: N_x=%d, N_y=%d, N_z=%d, for MPI size=%d with max prime factor %d\n", N_x,
						 N_y, N_z, mpi_size, max_prime);
			return;
		}

		std::printf("Decomposing domain of size %d x %d x %d = %d across %d MPI ranks...\n", N_x, N_y, N_z,
					N_x * N_y * N_z, mpi_size);

		DomainDecomposition dd(mpi_rank, mpi_size, N_x, N_y, N_z);
		const DomainInfo   &d = dd.info();

		std::printf("Rank %d/%d: local domain X[%d+%d] Y[%d+%d] Z[%d+%d]  grid: %d x %d x %d\n", mpi_rank, mpi_size,
					d.local_start_x, d.N_local_x, d.local_start_y, d.N_local_y, d.local_start_z, d.N_local_z, d.Px, d.Py,
					d.Pz);

		analyze_decomposition(d, mpi_size);

		auto best_infos = dd.best_infos(mpi_rank);

		DomainInfo best_info{};
		bool	   has_best_info	   = false;
		double	   best_aspect_ratio   = 0.0;
		double	   best_score		   = 0.0;
		double	   best_grid_ratio	   = 0.0;
		double	   best_load_imbalance = 0.0;

		for (const auto &info : best_infos)
		{
			const double aspect_ratio	= info.aspect_ratio_score();
			const int	 p_max			= std::max(info.Px, std::max(info.Py, info.Pz));
			const int	 p_min			= std::min(info.Px, std::min(info.Py, info.Pz));
			const double grid_ratio		= static_cast<double>(p_max) / p_min;
			const double load_imbalance = info.load_imbalance();

			const bool is_better = !has_best_info
								   || std::tie(load_imbalance, aspect_ratio, info.score_, grid_ratio)
										  < std::tie(best_load_imbalance, best_aspect_ratio, best_score, best_grid_ratio);

			std::printf(
				"C: aspect ratio: %.3f || score: %.3f || g-ratio max/min: %.3f || imbalance: %.3f || grid: %d x %d x "
				"%d%s\n",
				aspect_ratio, info.score_, grid_ratio, load_imbalance, info.Px, info.Py, info.Pz, is_better ? "  *" : "");

			if (is_better)
			{
				best_info			= info;
				best_aspect_ratio	= aspect_ratio;
				best_score			= info.score_;
				best_grid_ratio		= grid_ratio;
				best_load_imbalance = load_imbalance;
				has_best_info		= true;
			}
		}

		analyze_decomposition(best_info, mpi_size);

		if (has_best_info)
		{
			std::printf(
				"\n\n  Selected best info: aspect ratio = %.3f || score = %.3f || grid ratio max/min = %.3f || load "
				"imbalance = %.3f || grid: %d x %d x %d\n",
				best_aspect_ratio, best_score, best_grid_ratio, best_load_imbalance, best_info.Px, best_info.Py,
				best_info.Pz);
		}
	}

} // namespace domain_decomposition_3D
