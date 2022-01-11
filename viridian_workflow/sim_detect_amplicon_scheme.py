import csv
import itertools
import logging
import os

import pysam

from viridian_workflow import amplicon_schemes, detect_primers, minimap, utils


def amplicon_tsv_to_amplicon_coords(infile):
    coords = {}
    with open(infile) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for d in reader:
            pos = int(d["Position"])
            if d["Amplicon_name"] not in coords:
                coords[d["Amplicon_name"]] = {"starts": set(), "ends": set()}
            if d["Left_or_right"] == "left":
                coords[d["Amplicon_name"]]["starts"].add(pos)
            elif d["Left_or_right"] == "right":
                end = pos + len(d["Sequence"]) - 1
                coords[d["Amplicon_name"]]["ends"].add(end)
            else:
                raise Exception(f"Bad Left_or_right colum {d}")
    return coords


def amplicon_coords_to_fasta_of_seqs(ref_seq, scheme_name, coords, outfile):
    with open(outfile, "a") as f:
        for amplicon_name, d in coords.items():
            for (start, end) in itertools.product(d["starts"], d["ends"]):
                seq_name = f">{scheme_name}___{amplicon_name}___{start}___{end}"
                print(seq_name, file=f)
                print(ref_seq[start:end+1], file=f)


def amplicon_schemes_tsvs_to_fasta_of_seqs(ref_seq, schemes_tsvs, outfile):
    if os.path.exists(outfile):
        os.unlink(outfile)

    amplicon_coords = {}
    for scheme_name, scheme_tsv in schemes_tsvs.items():
        amplicon_coords[scheme_name] = amplicon_tsv_to_amplicon_coords(scheme_tsv)
        amplicon_coords_to_fasta_of_seqs(ref_seq, scheme_name, amplicon_coords[scheme_name], outfile)
    return amplicon_coords

def parse_sam_file(sam_file, amplicon_scheme_list, amplicon_coords):
    results = {}
    aln_file = pysam.AlignmentFile(sam_file, "r")
    for read in aln_file:
        print(read.query_name)
        if read.is_secondary or read.is_supplementary:
            continue
        scheme_name, amplicon_name, start, end = read.query_name.split("___")
        start = int(start)
        end = int(end)
        assert scheme_name in amplicon_coords
        assert amplicon_name in amplicon_coords[scheme_name]
        result_d = amplicon_coords[scheme_name][amplicon_name]
        if "matches" not in result_d:
            result_d["matches"] = {}
        result_key = (start, end)
        assert result_key not in result_d["matches"]
        matches = detect_primers.match_read_to_amplicons(read, amplicon_scheme_list)
        result_d["matches"][result_key] = {"amp_matches": matches}
        if matches is not None:
            for scheme_name, amplicon_list in matches.items():
                for amp in amplicon_list:
                    print(scheme_name, amp)
        print("_________________________________________________________")

def simulate_detect_scheme(
        built_in_amp_schemes,
        amp_schemes_tsv,
        ref_fasta,
        outprefix,
    ):
    logging.info("Gathering amplicons scheme files")
    if built_in_amp_schemes is None and amp_schemes_tsv is None:
        logging.info("No primer schemes provided. Using all built in schemes")
        built_in_amp_schemes = list(
            amplicon_schemes.get_built_in_schemes().keys()
        )
    (
        amplicon_scheme_name_to_tsv,
        amplicon_scheme_list,
    ) = amplicon_schemes.load_list_of_amplicon_sets(
        built_in_names_to_use=built_in_amp_schemes,
        tsv_others_to_use=amp_schemes_tsv,
    )

    logging.info("Making fasta of all amplicon sequences")
    ref_seq = utils.load_single_seq_fasta(ref_fasta)
    amplicon_seqs_fa = f"{outprefix}.amplicons.fa"
    amplicon_coords = amplicon_schemes_tsvs_to_fasta_of_seqs(ref_seq, amplicon_scheme_name_to_tsv, amplicon_seqs_fa)

    sam_file = f"{outprefix}.sam"
    logging.info("Mapping amplicons to reference")
    minimap.run(
        sam_file,
        ref_fasta,
        amplicon_seqs_fa,
        sort=False,
    )
    parse_sam_file(sam_file, amplicon_scheme_list, amplicon_coords)
    for scheme_name, results in amplicon_coords.items():
        for amp_name, d in results.items():
            for (start, end), match_dict in d["matches"].items():
                hits = sorted(list(match_dict["amp_matches"].keys()))
                print(scheme_name, amp_name, start, end, hits)
