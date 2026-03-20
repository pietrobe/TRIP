#ifndef DOMAIN_DECOMPOSITION_3D_HPP
#define DOMAIN_DECOMPOSITION_3D_HPP

// Require C++17 or later
static_assert(__cplusplus >= 201703L, "domain_decomposition_3D requires C++17 or later");

#include <mpi.h>
#include <deque>
#include <tuple>
#include <vector>

namespace domain_decomposition_3D
{

	// ---------------------------------------------------------------------------
	// Utility: prime factorization
	// ---------------------------------------------------------------------------

	/// Returns all prime factors of \p n (with repetition), in ascending order.
	std::vector<int>
	prime_factors(int n);

	// ---------------------------------------------------------------------------
	// DomainInfo — plain data holder for a single rank's decomposition result
	// ---------------------------------------------------------------------------

	void
	get_1d_decomposition(int n_global, int p_dims, int p_coord, int &n_local, int &start) noexcept;

	struct DomainInfo
	{
		// Processor grid dimensions
		int Px{1}, Py{1}, Pz{1};

		double score_ = 0.0; // Aspect ratio score for this decomposition (lower is better)

		// This rank's coordinate in the processor grid
		int ix{0}, iy{0}, iz{0};

		// Global domain extents
		int N_x{0}, N_y{0}, N_z{0};

		// Local subdomain extents
		int N_local_x{0}, N_local_y{0}, N_local_z{0};

		// Global start index of this rank's subdomain
		int local_start_x{0}, local_start_y{0}, local_start_z{0};

		// Rank layout: rank = ix*(Py*Pz) + iy*Pz + iz
		int
		rank_of(int px, int py, int pz) const noexcept;

		/// Returns the MPI rank of the neighbor at grid offset (dix, diy, diz)
		/// relative to this rank, or MPI_PROC_NULL at domain boundaries.
		int
		neighbor_rank(int dix, int diy, int diz) const noexcept;

		// Returns the number of cells in this rank's local subdomain, i.e. N_local_x * N_local_y * N_local_z.
		int
		number_of_cells() const noexcept
		{
			return N_local_x * N_local_y * N_local_z;
		}

		bool
		operator==(const DomainInfo &other) const noexcept
		{
			return this->score_ == other.score_;
		}

		bool
		operator<(const DomainInfo &other) const noexcept
		{
			return this->score_ < other.score_;
		}

		bool
		operator>(const DomainInfo &other) const noexcept
		{
			return this->score_ > other.score_;
		}

		bool
		operator!=(const DomainInfo &other) const noexcept
		{
			return !(*this == other);
		}

		bool
		operator<=(const DomainInfo &other) const noexcept
		{
			return !(*this > other);
		}

		bool
		operator>=(const DomainInfo &other) const noexcept
		{
			return !(*this < other);
		}

		inline void
		set_rank(int rank) noexcept
		{
			ix = rank / (Py * Pz);
			iy = (rank % (Py * Pz)) / Pz;
			iz = rank % Pz;
		}

		inline double
		load_imbalance() const noexcept
		{
			const double average_cells = static_cast<double>(N_x * N_y * N_z) / (Px * Py * Pz);
			return static_cast<double>(this->number_of_cells()) / average_cells;
		}

		inline double
		aspect_ratio_score() const noexcept
		{
			const double lx = static_cast<double>(N_x) / Px;
			const double ly = static_cast<double>(N_y) / Py;
			const double lz = static_cast<double>(N_z) / Pz;

			const double max_len = std::max(lx, std::max(ly, lz));
			const double min_len = std::min(lx, std::min(ly, lz));

			return max_len / min_len;
		}
	};

	// ---------------------------------------------------------------------------
	// DomainDecomposition — performs and owns the 3-D decomposition
	// ---------------------------------------------------------------------------

	class DomainDecomposition
	{
	   public:
		/// Construct and immediately perform the decomposition.
		///
		/// \param mpi_rank  This process's MPI rank.
		/// \param mpi_size  Total number of MPI processes.
		/// \param N_x, N_y, N_z  Global domain extents.
		DomainDecomposition(int mpi_rank, int mpi_size, int N_x, int N_y, int N_z);

		// Non-copyable, movable
		DomainDecomposition(const DomainDecomposition &) = delete;
		DomainDecomposition &
		operator=(const DomainDecomposition &) = delete;

		DomainDecomposition(DomainDecomposition &&) = default;
		DomainDecomposition &
		operator=(DomainDecomposition &&) = default;

		/// Read-only access to the decomposition result for this rank.
		const DomainInfo &
		info() const noexcept
		{
			return info_;
		}

		/// Builds decomposition info for any valid MPI rank and returns it by value.
		DomainInfo
		info_for_rank(int rank) const;

		std::deque<DomainInfo>
		best_infos(int rank) const;

		/// Convenience forwarder — neighbor rank at grid offset (dix, diy, diz).
		int
		neighbor_rank(int dix, int diy, int diz) const noexcept
		{
			return info_.neighbor_rank(dix, diy, diz);
		}

	   private:
		int		   mpi_rank_;
		int		   mpi_size_;
		DomainInfo info_;

		std::deque<DomainInfo> factorization_cache_;

		/// Chooses (Px, Py, Pz) that minimises the aspect-ratio score.
		void
		factorize_procs(int size, int Nx, int Ny, int Nz, int &Px, int &Py, int &Pz);

		/// Splits \p n_global cells across \p p_dims ranks;
		/// returns (n_local, start) for rank \p p_coord.
	};

	// ---------------------------------------------------------------------------
	// Stand-alone demo (mirrors demo_decompose_domain_3d)
	// ---------------------------------------------------------------------------

	/// Result from select_best_decomposition
	struct DecompositionSelection
	{
		DomainInfo best_info;
		double	   aspect_ratio	  = 0.0;
		double	   score		  = 0.0;
		double	   grid_ratio	  = 0.0;
		double	   load_imbalance = 0.0;
		bool	   found		  = false;
	};

	void
	analyze_decomposition(const DomainInfo &info, int mpi_size);

	/// Selects the best decomposition from a list of candidates based on:
	/// load imbalance, aspect ratio, score, and grid ratio (max/min processor count).
	/// If verbose is true, prints evaluation details for each candidate.
	DecompositionSelection
	select_best_decomposition(const std::deque<DomainInfo> &candidates, bool verbose = true);

	void
	demo_decompose_domain_3d();

} // namespace domain_decomposition_3D

#endif // DOMAIN_DECOMPOSITION_3D_HPP
