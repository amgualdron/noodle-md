#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# The SLURM job that ./launch-run submits. You never run this yourself.
#
# It is submitted as a job array: every array element runs one task of the run,
#     task id = <offset> + array index      (offset is 0 unless the run has more
#                                            tasks than the cluster's array limit)
# ./launch-run adds the job name, time limit, array range and output file.
# Arguments: <run folder> <offset>
# ─────────────────────────────────────────────────────────────────────────────
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# One thread per task (a single short chain gains nothing from more), and print
# progress straight into the output file.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1

TASK_ID=$(( $2 + SLURM_ARRAY_TASK_ID ))
echo "=== task $TASK_ID   job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}   $(hostname)   $(date) ==="

# The tool's own Python environment (made by setup.sh): no conda activate needed.
./env/bin/python src/simulate.py "$1" "$TASK_ID"
