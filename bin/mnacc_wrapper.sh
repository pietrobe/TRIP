#!/bin/bash



## /apps/ACC/GROMACS/test/wrapper.sh

lrank=${OMPI_COMM_WORLD_LOCAL_RANK:-${SLURM_LOCALID}}
lrank2=$(( lrank % 4 ))
case $lrank2 in
        0)
		export CUDA_VISIBLE_DEVICES=0
		export mycores="0-19"
		export UCX_NET_DEVICES=mlx5_0:1
		;;
        1)
		export CUDA_VISIBLE_DEVICES=1
		export mycores="20-39"
		export UCX_NET_DEVICES=mlx5_1:1
                ;;
        2)
		export CUDA_VISIBLE_DEVICES=2
		export mycores="40-59"
		export UCX_NET_DEVICES=mlx5_4:1
                ;;
        3)
		export CUDA_VISIBLE_DEVICES=3
		export mycores="60-79"
		export UCX_NET_DEVICES=mlx5_5:1
                ;;
        *)
                echo "ERROR: MPI local rank larger than expected! -> $lrank"
                exit 1
esac
# echo "numactl --physcpubind=$mycores $@ using CUDA $CUDA_VISIBLE_DEVICES"

numactl --physcpubind=$mycores $@

