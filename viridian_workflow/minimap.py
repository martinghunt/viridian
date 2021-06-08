"""minimap2 and samtools utilities
"""
import os
import subprocess
from viridian_workflow.utils import run_process, check_file


def run_se(outdir, ref_genome, fq, prefix=None, threads=1):
    """map single fastq ONT reads with minimap2
    """
    bam = os.path.join(outdir, "reference_mapped.bam")
    if prefix:
        bam = os.path.join(outdir, f"{prefix}-reference_mapped.bam")

    minimap_cmd = ["minimap2", "-t", str(threads), "-ax", "map-ont", ref_genome, fq]
    sort_cmd = ["samtools", "sort", "-o", bam]

    map_proc = subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE)
    sort_proc = subprocess.Popen(sort_cmd, stdin=map_proc.stdout)
    sort_proc.wait()
    check_file(bam)

    run_process(f"samtools index {bam}")
    check_file(bam + ".bai")
    return bam


def run(outdir, ref_genome, fq1, fq2, prefix=None, threads=1):
    bam = os.path.join(outdir, "reference_mapped.bam")
    if prefix:
        bam = os.path.join(outdir, f"{prefix}-reference_mapped.bam")

    minimap_cmd = ["minimap2", "-t", str(threads), "-ax", "sr", ref_genome, fq1, fq2]
    sort_cmd = ["samtools", "sort", "-o", bam]

    map_proc = subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE)
    sort_proc = subprocess.Popen(sort_cmd, stdin=map_proc.stdout)
    sort_proc.wait()
    check_file(bam)

    run_process(f"samtools index {bam}")
    check_file(bam + ".bai")
    return bam
