from viridian_workflow import sim_detect_amplicon_scheme

def run(options):
    sim_detect_amplicon_scheme.simulate_detect_scheme(
        options.built_in_amp_schemes,
        options.amp_schemes_tsv,
        options.ref_fasta,
        options.outprefix,
    )

