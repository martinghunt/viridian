import logging

from viridian import qc_plots, utils


def gather_data(options):
    p = qc_plots.Plots()
    p.to_pickle(options.index_tsv, options.outfile, cpus=options.cpus)
    if options.tsv_dir is not None:
        logging.info(
            f"--tsv_dir option used, also writing data to TSV files in {options.tsv_dir}"
        )
        p.write_plot_data_tsvs(options.tsv_dir)


def plot(options):
    plots = qc_plots.Plots()
    plots.from_pickle(options.infile)
    if options.force:
        logging.info(f"--force option used, so deleting {options.outdir} if it exists")
        utils.syscall(f"rm -rf {options.outdir}")
    plots.plot(
        options.outdir,
        x_start=None,
        x_end=None,
        colours=None,
        x_window=None,
        plot_width=1000,
        dataset_height=200,
        gene_track_height=20,
        y_gap=45,
        x_tick_step=None,
    )
