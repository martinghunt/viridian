import filecmp
import json
import os
import pytest
import subprocess

from viridian_workflow import sim_detect_amplicon_scheme

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "sim_detect_amplicon_scheme")


def test_amplicon_tsv_to_amplicon_coords():
    amplicon_tsv = os.path.join(data_dir, "amplicon_tsv_to_amplicon_coords.tsv")
    got = sim_detect_amplicon_scheme.amplicon_tsv_to_amplicon_coords(amplicon_tsv)
    expect = {
        "amp_1": {"ends": {430}, "starts": {25, 26}},
        "amp_2": {"ends": {728, 726}, "starts": {324}},
    }
    assert expect == got


def test_amplicon_coords_to_fasta_of_seqs():
    ref_seq = "AGTCGATCGATCG"
    scheme_name = "SCHEME"
    coords = {
        "amp_1": {"starts": {0, 1}, "ends": {9}},
        "amp_2": {"starts": {5,6}, "ends": {11,12}},
    }
    tmp_out = "tmp.amplicon_coords_to_fasta_of_seqs.fa"
    subprocess.check_output(f"rm -rf {tmp_out}", shell=True)
    sim_detect_amplicon_scheme.amplicon_coords_to_fasta_of_seqs(ref_seq, scheme_name, coords, tmp_out)
    expect = os.path.join(data_dir, "amplicon_coords_to_fasta_of_seqs.fa")
    assert filecmp.cmp(expect, tmp_out, shallow=False)
    os.unlink(tmp_out)


def test_amplicon_schemes_tsvs_to_fasta_of_seqs():
    # FIXME
    #sim_detect_amplicon_scheme.amplicon_schemes_tsvs_to_fasta_of_seqs(ref_seq, schemes_tsvs, outfile)
    pass


def test_parse_sam_file():
    # FIXME
    #sim_detect_amplicon_scheme.parse_sam_file()
    pass


def test_simulate_detect_scheme():
    # FIXME
    #sim_detect_amplicon_scheme.simulate_detect_scheme(
    pass
