"""The pipeline definition
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional, Any

# import tempfile
import json

from viridian_workflow import readstore, self_qc
from viridian_workflow.subtasks import Cylon, Minimap, Varifier
from viridian_workflow.primers import AmpliconSet


def run_pipeline(
    work_dir: Path,
    platform: str,
    fqs: list[Path],
    amplicon_sets: list[AmpliconSet],
    ref: Path = Path("../covid/MN908947.fasta"),
    force_amp_scheme: Optional[AmpliconSet] = None,
    keep_intermediate: bool = False,
    keep_bam: bool = False,
    sample_name: str = "sample",
    frs_threshold: float = 0.1,
    self_qc_depth: int = 20,
    consensus_max_n_percent: int = 50,
    max_percent_amps_fail: int = 50,
    dump_tsv: bool = False,
    command_line_args: Optional[dict[str, Any]] = None,
    force_consensus: Optional[Path] = None,
    global_log: dict[str, Any] = {},  # global pipeline log dictionary (bad)
):

    results: dict[str, Any] = {}
    results["Progress"] = []
    work_dir = Path(work_dir)
    if work_dir.exists():
        raise Exception(f"Output directory {work_dir} already exists")
    work_dir.mkdir()

    # generate name-sorted bam from fastqs
    if platform == "illumina":
        fq1, fq2 = fqs
        minimap = Minimap(work_dir / "name_sorted.bam", ref, fq1, fq2=fq2, sort=False)
    elif platform == "ont":
        fq = fqs[0]
        minimap = Minimap(work_dir / "name_sorted.bam", ref, fq, sort=False)
    elif platform == "iontorrent":
        raise NotImplementedError
    else:
        print(f"Platform {platform} is not supported.", file=sys.stderr)
        exit(1)

    unsorted_bam: Path = minimap.run()
    results["Progress"].append(minimap.log)

    # add minimap task log to result log
    # log["minimap"] = minimap.log

    # pre-process input bam
    bam: readstore.Bam = readstore.Bam(unsorted_bam)
    results["Coverage"] = {
        "total_reads": 0,  # TODO
        "Total_fragments": 0,  # TODO
        "Reference_coverage": 0,  # TODO
        "Reference_length": 0,  # TODO
        "Average_amplicon_depth": 0,  # TODO
    }

    # detect amplicon set
    amplicon_set: AmpliconSet = bam.detect_amplicon_set(amplicon_sets)
    results["Amplicons"] = {
        "scheme": amplicon_set.name,
        "total_amplicons": len(amplicon_set.amplicons),
        "Successful_amplicons": 0,  # TODO
        "fragment_matches": 0,  # TODO
        "fragment_mismatches": 0,  # TODO
    }

    # construct readstore
    # this subsamples the reads
    reads = (
        readstore.ReadStore(amplicon_set, bam)
        if force_amp_scheme is None
        else readstore.ReadStore(force_amp_scheme, bam)
    )

    # log["amplicons"] = reads.summary

    # branch on whether to run cylon or use external assembly ("cuckoo mode")
    consensus: Optional[Path] = None
    if force_consensus is not None:
        # global_log["forced_consensus"] = str(force_consensus)
        consensus = Path(force_consensus)

    else:
        # save reads for cylon assembly
        amp_dir = work_dir / "amplicons"
        manifest_data = reads.make_reads_dir_for_cylon(amp_dir)

        # run cylon
        cylon = Cylon(work_dir, platform, ref, amp_dir, manifest_data, reads.cylon_json)
        consensus = cylon.run()
        results["Progress"].append(cylon.log)

    # satify type bounds and ensure the readstore was properly constructed
    assert consensus is not None
    assert reads.start_pos is not None
    assert reads.end_pos is not None

    # varifier
    varifier = Varifier(
        work_dir / "varifier",
        ref,
        consensus,
        min_coord=reads.start_pos,
        max_coord=reads.end_pos,
    )
    vcf, msa, varifier_consensus = varifier.run()
    results["Progress"].append(varifier.log)

    pileup = self_qc.Pileup(
        varifier_consensus,
        reads,
        msa=msa,
        config=self_qc.Config(frs_threshold, self_qc_depth),
    )

    # masked fasta output
    masked_fasta: str = pileup.mask()
    # log["self_qc"] = pileup.log
    # log["qc"] = pileup.summary

    results["Self_qc"] = {
        "Masked_by_assembler": 0,  # TODO
        "Total_masked_incl_self_qc": 0,  # TODO
        "Filters": {},  # TODO
    }

    results["Consensus"] = ("",)  # TODO
    results["reference_start"] = (0,)  # TODO
    results["reference_end"] = (0,)  # TODO

    with open(work_dir / "consensus.fa", "w", encoding="utf-8") as fasta_out:
        print(f">{sample_name}", file=fasta_out)
        print(masked_fasta, file=fasta_out)

    # annotate vcf
    annotated_vcf = pileup.annotate_vcf(vcf)

    # dump tsv
    if dump_tsv:
        _ = pileup.dump_tsv(work_dir / "all_stats.tsv", amplicon_set)

    with open(work_dir / "final.vcf", "w", encoding="utf-8") as vcf_out:
        header, records = annotated_vcf
        for h in header:
            print(h, file=vcf_out)
        for rec in records:
            print("\t".join(map(str, rec)), file=vcf_out)

    return results
