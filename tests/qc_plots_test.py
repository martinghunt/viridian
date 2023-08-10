import copy
import os
import pytest

from viridian import qc_plots, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "qc_plots")


def test_qc_dict_line_to_basecall():
    default_d = {
        "Ref_nt": "A",
        "Cons_nt": "A",
        "Cons_depth": 0,
        "Mask": ".",
    }
    default_d.update({x: 0 for x in ["A", "C", "G", "T", "a", "c", "g", "t"]})

    d = copy.copy(default_d)
    d["Ref_nt"] = "-"
    assert qc_plots.qc_dict_line_to_basecall(d) == qc_plots.Basecall.INDEL
    d["Cons_nt"] = "-"
    assert qc_plots.qc_dict_line_to_basecall(d) == qc_plots.Basecall.INDEL
    d["Ref_nt"] = "A"
    assert qc_plots.qc_dict_line_to_basecall(d) == qc_plots.Basecall.INDEL

    d = copy.copy(default_d)
    d["Cons_nt"] = "N"
    assert qc_plots.qc_dict_line_to_basecall(d) == qc_plots.Basecall.N

    d = copy.copy(default_d)
    d["a"] = 100
    d["A"] = 100
    d["Mask"] = "PASS"
    assert qc_plots.qc_dict_line_to_basecall(d) == qc_plots.Basecall.ACGT_GOOD
    d["Mask"] = "DEPTH"
    assert qc_plots.qc_dict_line_to_basecall(d) == qc_plots.Basecall.ACGT_DP
    # when ref is A/C/G/T and cons is A/C/G/T, but consensus is better supported
    # by a different base, mask is FRS
    d["Mask"] = "FRS"
    assert qc_plots.qc_dict_line_to_basecall(d) == qc_plots.Basecall.ACGT_GOOD
    # Make more Cs than consensus call of A, and should get ACGT_BAD
    d["c"] = 101
    d["C"] = 101
    assert qc_plots.qc_dict_line_to_basecall(d) == qc_plots.Basecall.ACGT_BAD

    d = copy.copy(default_d)
    # IUPAC R is A or G
    d["Cons_nt"] = "R"
    d["Mask"] = "PASS"
    d["a"] = 50
    d["A"] = 50
    d["g"] = 60
    d["G"] = 60
    assert qc_plots.qc_dict_line_to_basecall(d) == qc_plots.Basecall.HET_GOOD
    d["Mask"] = "DEPTH"
    assert qc_plots.qc_dict_line_to_basecall(d) == qc_plots.Basecall.HET_DP
    d["Mask"] = "HET"
    assert qc_plots.qc_dict_line_to_basecall(d) == qc_plots.Basecall.HET_GOOD
    d["c"] = 70
    d["C"] = 70
    assert qc_plots.qc_dict_line_to_basecall(d) == qc_plots.Basecall.HET_BAD
    d["c"] = 55
    d["C"] = 55
    assert qc_plots.qc_dict_line_to_basecall(d) == qc_plots.Basecall.HET_GOOD


def test_load_qc_tsv():
    tmp_file = os.path.join(data_dir, "load_qc_tsv.tsv")
    assert qc_plots.load_qc_tsv(tmp_file) == [
        qc_plots.Basecall.N,
        qc_plots.Basecall.INDEL,
        qc_plots.Basecall.ACGT_GOOD,
    ]


def test_GenomeCalls():
    tmp_dir = "tmp.test_qc_plots.GenomeCalls"
    utils.syscall(f"rm -rf {tmp_dir}")
    os.mkdir(tmp_dir)
    tsv_fofn = os.path.join(tmp_dir, "tsvs.fofn")
    with open(tsv_fofn, "w") as f:
        for i in range(1, 4, 1):
            print(os.path.join(data_dir, f"genome_calls.qc.{i}.tsv"), file=f)
    gc = qc_plots.GenomeCalls()
    gc.add_qc_tsvs(tsv_fofn, cpus=2)
    expect_calls = [
        [0, 0, 0, 3, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 2],
        [3, 0, 0, 0, 0, 0, 0, 0],
        [2, 0, 0, 1, 0, 0, 0, 0],
        [3, 0, 0, 0, 0, 0, 0, 0],
        [2, 0, 0, 0, 0, 0, 0, 1],
    ]
    assert gc.calls == expect_calls

    calls_file = os.path.join(tmp_dir, "calls.txt.gz")
    gc.to_text_file(calls_file)
    assert os.path.exists(calls_file)
    gc = qc_plots.GenomeCalls()
    gc.from_text_file(calls_file)
    assert gc.calls == expect_calls

    calls_file = os.path.join(tmp_dir, "calls.pickle")
    gc.to_pickle(calls_file)
    assert os.path.exists(calls_file)
    gc = qc_plots.GenomeCalls()
    gc.from_pickle(calls_file)
    assert gc.calls == expect_calls

    other_gc = qc_plots.GenomeCalls()
    other_gc.calls = [
        [0, 1, 0, 0, 0, 0, 0, 1],
        [0, 2, 0, 0, 0, 0, 0, 2],
        [0, 3, 0, 0, 0, 0, 0, 0],
        [4, 0, 4, 0, 0, 0, 0, 0],
        [0, 0, 0, 5, 0, 0, 0, 0],
        [0, 0, 0, 0, 6, 1, 2, 0],
    ]
    gc.add_other_genome_calls(other_gc)
    expect_calls = [
        [0, 1, 0, 3, 0, 0, 0, 1],
        [1, 2, 0, 0, 0, 0, 0, 4],
        [3, 3, 0, 0, 0, 0, 0, 0],
        [6, 0, 4, 1, 0, 0, 0, 0],
        [3, 0, 0, 5, 0, 0, 0, 0],
        [2, 0, 0, 0, 6, 1, 2, 1],
    ]
    assert gc.calls == expect_calls

    y_vals = gc.get_plot_y_vals()
    # not worth checking this in too much detail, as y value order etc will
    # probably change. Just check dimensions are correct
    assert len(y_vals) == len(qc_plots.Basecall)
    for l in y_vals:
        assert len(l) == len(gc)

    utils.syscall(f"rm -r {tmp_dir}")


def test_Plots():
    tmp_dir = "tmp.test_qc_plots.Plots"
    utils.syscall(f"rm -rf {tmp_dir}")
    os.mkdir(tmp_dir)
    index_tsv = os.path.join(tmp_dir, "index.tsv")
    genome_calls_tsvs = [
        os.path.join(data_dir, f"genome_calls.qc.{i}.tsv") for i in [1, 2, 3]
    ]
    datasets = {
        "set1": {
            "fofn": "qc_names_file.1.txt",
            "tsv_files": genome_calls_tsvs,
        },
        "set2": {
            "fofn": "qc_names_file.2.txt",
            "tsv_files": genome_calls_tsvs[:2],
        },
    }
    call_sets = {}
    with open(index_tsv, "w") as f:
        print("name", "qc_names_file", sep="\t", file=f)
        for name, d in datasets.items():
            print(name, os.path.join(tmp_dir, d["fofn"]), sep="\t", file=f)
            fofn = os.path.join(tmp_dir, d["fofn"])
            with open(fofn, "w") as f2:
                print(*d["tsv_files"], sep="\n", file=f2)
            call_sets[name] = qc_plots.GenomeCalls()
            call_sets[name].add_qc_tsvs(fofn)

    plots = qc_plots.Plots()
    pickle_file = os.path.join(tmp_dir, "data.pickle")
    plots.from_index_file(index_tsv)
    plots.to_pickle(pickle_file)
    assert plots.genome_calls == call_sets

    pickles_fofn = os.path.join(tmp_dir, "data.pickles.fofn")
    with open(pickles_fofn, "w") as f:
        print(pickle_file, file=f)
        print(pickle_file, file=f)
    plots = qc_plots.Plots()
    pickle_file_combined = pickle_file + ".combined"
    plots.combine_pickles(pickles_fofn, pickle_file_combined)
    assert os.path.exists(pickle_file_combined)
    call_sets_combined = copy.deepcopy(call_sets)
    for k in call_sets:
        call_sets_combined[k].add_other_genome_calls(call_sets[k])
    assert plots.genome_calls == call_sets_combined

    plots = qc_plots.Plots()
    plots.from_pickle(pickle_file)
    assert plots.genome_calls == call_sets

    plot_tsvs_dir = os.path.join(tmp_dir, "Plot_tsvs")
    plots.write_plot_data_tsvs(plot_tsvs_dir)
    assert os.path.exists(plot_tsvs_dir)
    got_index = utils.load_json(os.path.join(plot_tsvs_dir, "index.json"))
    assert got_index == {"set1": "0.tsv.gz", "set2": "1.tsv.gz"}
    for filename in got_index.values():
        assert os.path.exists(os.path.join(plot_tsvs_dir, filename))

    plots_outdir = os.path.join(tmp_dir, "Plots")
    plots.plot(plots_outdir)
    assert os.path.exists(plots_outdir)
    for ext in ["svg", "png", "pdf"]:
        assert os.path.exists(os.path.join(plots_outdir, f"plot.{ext}"))
    utils.syscall(f"rm -r {tmp_dir}")
