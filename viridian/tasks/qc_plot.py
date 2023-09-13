import logging

from viridian import qc_plots, utils


def gather_data(options):
    p = qc_plots.Plots()
    p.from_index_file(options.index_tsv, cpus=options.cpus)
    p.to_pickle(options.outfile)
    if options.tsv_dir is not None:
        logging.info(
            f"--tsv_dir option used, also writing data to TSV files in {options.tsv_dir}"
        )
        p.write_plot_data_tsvs(options.tsv_dir)


def combine_data(options):
    p = qc_plots.Plots()
    p.combine_pickles(options.infile, options.outfile)


def plot(options):
    plots = qc_plots.Plots()
    plots.from_pickle(options.infile)
    if options.force:
        logging.info(f"--force option used, so deleting {options.outdir} if it exists")
        utils.syscall(f"rm -rf {options.outdir}")
    plots.plot(
        options.outdir,
        x_start=None if options.range is None else options.range[0] - 1,
        x_end=None if options.range is None else options.range[1] - 1,
        colours=None,
        x_window=None,
        plot_width=1000,
        dataset_height=120,
        genes_track_height=20,
        y_gap=45,
        x_tick_step=None,
        amp_scheme=options.amp_scheme,
        plot_amp_names=options.add_amp_names,
        plot_amp_number=options.add_amp_number,
        plot_genes=options.gene_track,
        plot_primers=options.add_primers,
        title=options.title,
        diff_track=options.plot_diff,
        datasets_to_plot=options.dataset,
        hist_datasets=options.hist_diffs,
        stats_tracks=options.stat,
        stats_tracks_bin=options.stat_bin_width,
    )
