import csv
from enum import Enum
import gzip
import logging
import math
import multiprocessing
import os

from viridian import constants, utils


class Basecall(Enum):
    ACGT_GOOD = 0  # cons is A/C/G/T, reads support the call
    ACGT_DP = 1  # cons is A/C/G/T, position has DEPTH mask
    ACGT_BAD = 2  # cons is A/C/G/T, but reads support different base (non-N)
    INDEL = 3  # indel (ref or cons are "-")
    HET_GOOD = 4  # cons is het, reads support the call
    HET_DP = 5  # cons is het, position has DEPTH mask
    HET_BAD = 6  # cons is het, but reads support something else (non-N)
    N = 7  # cons is N


PLOT_DEFAULT_COLOURS = {
    "ACGT_GOOD": "green",
    "ACGT_DP": "pink",
    "ACGT_BAD": "red",
    "INDEL": "gray",
    "HET_GOOD": "lightgreen",
    "HET_DP": "lightblue",
    "HET_BAD": "blue",
    "N": "lightgray",
}


# from bottom to top
PLOT_BASECALL_ORDER = [
    "ACGT_GOOD",
    "HET_GOOD",
    "INDEL",
    "N",
    "HET_DP",
    "HET_BAD",
    "ACGT_DP",
    "ACGT_BAD",
]


def svg_export(infile, outfile):
    utils.syscall(f"inkscape {infile} -o {outfile}")


def svg_header_lines(width, height):
    return [
        '<?xml version="1.0" standalone="no"?>',
        '<!DOCTYPE svg PUBLIC " -//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">',
        f'<svg width="{width}" height="{height}">',
    ]


def svg_line(x1, y1, x2, y2, colour, stroke_width=1):
    return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{colour}" stroke-width="{stroke_width}" />'


def svg_polygon(
    fill_colour,
    border_colour,
    xy_coords=None,
    x_coords=None,
    y_coords=None,
    border_width=1,
    opacity=-1,
):
    coord_to_join = xy_coords if xy_coords is not None else zip(x_coords, y_coords)
    l = [
        '<polygon points="',
        " ".join([f"{x[0]},{x[1]}" for x in coord_to_join]) + '"',
        f'fill="{fill_colour}"',
    ]

    if opacity != -1:
        l.append(f'fill-opacity="{opacity}"')

    l.extend(
        [
            f'stroke="{border_colour}"',
            f'stroke-width="{border_width}"',
            "/>",
        ]
    )
    return " ".join(l)


def svg_text(
    x_coord,
    y_coord,
    text,
    colour="black",
    vertical=False,
    font_size=12,
    h_right=False,
    h_center=False,
    v_center=False,
    font_style="normal",
    font_weight="normal",
):
    l = [
        f'<text x="{x_coord}" y="{y_coord+1}"',
        f'font-size="{font_size}"',
        f'font-style="{font_style}"',
        f'font-weight="{font_weight}"',
        f'fill="{colour}"',
    ]
    if h_center:
        l.append('text-anchor="middle"')
    elif h_right:
        l.append('text-anchor="end"')
    if v_center:
        l.append('dominant-baseline="middle"')
    if vertical:
        l.append(f'transform="rotate(-90,{x_coord},{y_coord})"')
    return " ".join(l) + f">{text}</text>"


def tick_positions_from_range(minimum, maximum, step=None):
    if step is None:
        magnitude = math.floor(math.log(maximum - minimum, 10))
        step = 10**magnitude
        while step * 3 > maximum - minimum:
            step = int(step / 2)
        step = max(step, 1)
    first = step * (minimum // step)
    return [x for x in range(first, maximum + 1, step) if minimum <= x <= maximum]


def x_axis(
    left,
    right,
    y,
    tick_labels,
    tick_positions,
    axis_colour="dimgray",
    tick_colour="dimgray",
    text_colour="dimgray",
):
    axis_thickness = 1
    tick_thickness = 0.5
    tick_height = 3
    lines = [
        svg_line(left, y, right, y, colour=axis_colour, stroke_width=axis_thickness)
    ]

    for xpos, xlab in zip(tick_positions, tick_labels):
        lines.extend(
            [
                svg_line(
                    xpos,
                    y,
                    xpos,
                    y + tick_height,
                    colour=tick_colour,
                    stroke_width=tick_thickness,
                ),
                svg_text(
                    xpos,
                    y + tick_height + 9,
                    xlab,
                    font_size=10,
                    h_center=True,
                    colour=text_colour,
                ),
            ]
        )
    return lines


def y_axis(
    top,
    bottom,
    x,
    tick_labels,
    tick_positions,
    axis_colour="dimgray",
    tick_colour="dimgray",
    text_colour="dimgray",
):
    axis_thickness = 1
    tick_thickness = 0.5
    tick_height = 3
    lines = [
        svg_line(x, top, x, bottom, colour=axis_colour, stroke_width=axis_thickness)
    ]

    for ypos, ylab in zip(tick_positions, tick_labels):
        lines.extend(
            [
                svg_line(
                    x - tick_height,
                    ypos,
                    x,
                    ypos,
                    colour=tick_colour,
                    stroke_width=tick_thickness,
                ),
                svg_text(
                    x - tick_height - 2,
                    ypos,
                    ylab,
                    font_size=10,
                    h_right=True,
                    v_center=True,
                    colour=text_colour,
                ),
            ]
        )
    return lines


def qc_dict_to_best_other_depths(d, exclude_base):
    return max(
        int(d[c]) + int(d[c.lower()]) for c in constants.ACGT.difference(exclude_base)
    )


def qc_dict_line_to_basecall(d):
    if d["Ref_nt"] == "-" or d["Cons_nt"] == "-":
        return Basecall.INDEL

    cons = d["Cons_nt"]

    if cons == "N":
        return Basecall.N
    elif cons in constants.ACGT:
        if d["Mask"] == "PASS":
            return Basecall.ACGT_GOOD
        elif "DEPTH" in d["Mask"]:
            return Basecall.ACGT_DP
        elif int(d["Cons_depth"]) >= qc_dict_to_best_other_depths(d, cons):
            return Basecall.ACGT_GOOD
        else:
            return Basecall.ACGT_BAD
    else:
        assert cons in constants.IUPAC
        if "DEPTH" in d["Mask"]:
            return Basecall.HET_DP
        elif d["Mask"] == "PASS":
            return Basecall.HET_GOOD

        cons_bases = set(constants.IUPAC[cons])
        for c in cons_bases:
            c_depth = int(d[c]) + int(d[c.lower()])
            if c_depth >= qc_dict_to_best_other_depths(d, c):
                return Basecall.HET_GOOD
        return Basecall.HET_BAD


def load_qc_tsv(filename):
    results = []
    previous_ref_pos = None
    of = gzip.open if filename.endswith(".gz") else open
    with of(filename, "rt") as f:
        for d in csv.DictReader(f, delimiter="\t"):
            # if we see the same ref position more than once, it's an indel
            if previous_ref_pos is not None and previous_ref_pos == d["Ref_pos"]:
                results[-1] = Basecall.INDEL
            else:
                previous_ref_pos = d["Ref_pos"]
                results.append(qc_dict_line_to_basecall(d))
    return results


class GenomeCalls:
    def __init__(self):
        self.calls = []

    def __eq__(self, other):
        return type(other) is type(self) and self.__dict__ == other.__dict__

    def __len__(self):
        return len(self.calls)

    def add_one_genome(self, new_results):
        if len(self.calls) == 0:
            for call in new_results:
                self.calls.append([0] * len(Basecall))
                self.calls[-1][call.value] += 1
        else:
            assert len(self.calls) == len(new_results)
            for i, call in enumerate(new_results):
                self.calls[i][call.value] += 1

    def add_list_of_genomes(self, new_genomes):
        for new_genome in new_genomes:
            self.add_one_genome(new_genome)

    def to_text_file(self, outfile):
        logging.info(f"Writing data to file {outfile}")
        of = gzip.open if outfile.endswith(".gz") else open
        with of(outfile, "wt") as f:
            print(*[x.name for x in Basecall], sep="\t", file=f)
            for l in self.calls:
                print(*l, sep="\t", file=f)

    def from_text_file(self, infile):
        logging.debug(f"Loading file {infile}")
        self.calls = []
        first = True
        of = gzip.open if infile.endswith(".gz") else open
        with of(infile, "rt") as f:
            for line in f:
                if first:
                    assert len(line.split()) == len(Basecall)
                    first = False
                else:
                    self.calls.append([int(i) for i in line.rstrip().split()])
        logging.debug(f"Loaded file {infile}")

    def to_pickle(self, outfile):
        utils.write_pickle(outfile, self.calls)

    def from_pickle(self, infile):
        self.calls = utils.load_pickle(infile)

    def add_qc_tsvs(self, file_of_qc_filenames, cpus=1):
        logging.info(f"Getting QC files to load from {file_of_qc_filenames}")
        with open(file_of_qc_filenames) as f:
            all_files = [x.rstrip() for x in f]
        logging.info(f"{len(all_files)} QC files to load from {file_of_qc_filenames}")
        next_debug_i = 0.1 * len(all_files)

        for i in range(0, len(all_files), cpus):
            with multiprocessing.Pool(processes=cpus) as p:
                new_results = p.map(
                    load_qc_tsv,
                    all_files[i : i + cpus],
                )
            self.add_list_of_genomes(new_results)
            if i >= next_debug_i:
                logging.info(
                    f"Loaded {min(i+cpus, len(all_files))}/{len(all_files)} files from {file_of_qc_filenames}"
                )
                next_debug_i += 0.1 * len(all_files)

        logging.info(f"Finished loading files from {file_of_qc_filenames}")

    def get_plot_y_vals(self, x_start=None, x_end=None):
        y_vals = [[] for _ in range(len(Basecall))]
        if x_end is None:
            x_end = len(self) - 1
        for counts in self.calls[x_start : x_end + 1]:
            total = 0
            for name in PLOT_BASECALL_ORDER:
                total += counts[Basecall[name].value]
                y_vals[Basecall[name].value].append(total)
        return y_vals


class Plots:
    def __init__(self):
        self.genome_calls = {}
        self.ref_length = None

    def add_genome_call_set(self, name, genome_calls):
        assert name not in self.genome_calls
        self.genome_calls[name] = genome_calls

    def add_genomes(self, index_file, cpus=1):
        with open(index_file) as f:
            for d in csv.DictReader(f, delimiter="\t"):
                call_set = GenomeCalls()
                call_set.add_qc_tsvs(d["qc_names_file"], cpus=cpus)
                if self.ref_length is None:
                    self.ref_length = len(call_set)
                else:
                    assert self.ref_length == len(call_set)
                self.add_genome_call_set(d["name"], call_set)
                logging.info(
                    f"Loaded data set {d['name']} from file {d['qc_names_file']}"
                )

    def to_pickle(self, index_file, outfile, cpus=1):
        self.add_genomes(index_file, cpus=cpus)
        utils.write_pickle(outfile, self.genome_calls)
        logging.info(f"Data written to file {outfile}")

    def from_pickle(self, infile):
        self.genome_calls = utils.load_pickle(infile)
        genome_lengths = set(len(x) for x in self.genome_calls.values())
        assert len(genome_lengths) == 1
        self.ref_length = genome_lengths.pop()
        logging.info(f"Data loaded from file {infile}")

    def write_plot_data_tsvs(self, outdir):
        os.mkdir(outdir)
        file_lookup = {}
        for name, call_set in self.genome_calls.items():
            file_lookup[name] = f"{len(file_lookup)}.tsv.gz"
            call_set.to_text_file(os.path.join(outdir, file_lookup[name]))
            logging.info(f"Finished data set{name}, filename: {file_lookup[name]}")
        utils.write_json(os.path.join(outdir, "index.json"), file_lookup)
        logging.info("All done")

    def plot(
        self,
        outdir,
        x_start=None,
        x_end=None,
        colours=None,
        x_window=None,
        plot_width=1000,
        dataset_height=200,
        gene_track_height=20,
        y_gap=50,
        x_tick_step=None,
        y_tick_step=None,
    ):
        bottom_gap = 50
        if colours is None:
            colours = {}
        assert len(self.genome_calls) > 0
        if x_start is None:
            x_start = 0
        if x_end is None:
            x_end = self.ref_length - 1
        assert 0 <= x_start < x_end < self.ref_length
        colours = {
            Basecall[k].name: colours.get(k, v) for k, v in PLOT_DEFAULT_COLOURS.items()
        }
        total_height = (dataset_height + y_gap) * len(self.genome_calls) + bottom_gap
        plot_rect_left_x = 100
        plot_rect_right_x = plot_width - 100
        rect_width = plot_rect_right_x - plot_rect_left_x
        rect_y_top = y_gap
        rect_y_bottom = rect_y_top + dataset_height
        rect_x_scale = rect_width / (x_end - x_start)
        rect_middle = 0.5 * (plot_rect_left_x + plot_rect_right_x)
        x_ticks_abs = tick_positions_from_range(x_start, x_end, x_tick_step)
        x_ticks_pos = [plot_rect_left_x + x * rect_x_scale for x in x_ticks_abs]
        print("x_ticks_abs", x_ticks_abs)
        print("x_ticks_pos", x_ticks_pos)

        os.mkdir(outdir)
        svg_out = os.path.join(outdir, "plot.svg")
        with open(svg_out, "w") as f:
            print(*svg_header_lines(plot_width, total_height), sep="\n", file=f)
            x_coords = [
                plot_rect_left_x + x * rect_x_scale
                for x in [0] + list(range(x_start, x_end + 1)) + [x_end, 0]
            ]

            for set_name, calls in self.genome_calls.items():
                print(
                    svg_text(rect_middle, rect_y_top - 8, set_name, h_center=True),
                    file=f,
                )
                y_vals = calls.get_plot_y_vals(x_start=x_start, x_end=x_end)
                max_y = max(y_vals[Basecall[PLOT_BASECALL_ORDER[-1]].value])

                for base_type in reversed(PLOT_BASECALL_ORDER):
                    y_coords = [rect_y_bottom] + [
                        rect_y_bottom - (y * dataset_height / max_y)
                        for y in y_vals[Basecall[base_type].value]
                    ]
                    y_coords.extend([rect_y_bottom, rect_y_bottom])
                    print(
                        svg_polygon(
                            colours[base_type],
                            colours[base_type],
                            x_coords=x_coords,
                            y_coords=y_coords,
                            border_width=0,
                        ),
                        file=f,
                    )

                y_ticks_abs = tick_positions_from_range(0, max_y, y_tick_step)
                y_ticks_pos = [
                    rect_y_bottom - (y * dataset_height / max_y) for y in y_ticks_abs
                ]
                print(
                    *y_axis(
                        rect_y_top,
                        rect_y_bottom,
                        plot_rect_left_x,
                        y_ticks_abs,
                        y_ticks_pos,
                    ),
                    sep="\n",
                    file=f,
                )
                print(
                    *x_axis(
                        plot_rect_left_x,
                        plot_rect_right_x,
                        rect_y_bottom,
                        x_ticks_abs,
                        x_ticks_pos,
                    ),
                    sep="\n",
                    file=f,
                )
                rect_y_top = rect_y_bottom + y_gap
                rect_y_bottom = rect_y_top + dataset_height

            last_rect_bottom = rect_y_top - y_gap
            print(
                svg_text(
                    rect_middle,
                    last_rect_bottom + 20,
                    "Position in genome",
                    h_center=True,
                    v_center=True,
                    colour="dimgrey",
                ),
                file=f,
            )
            print(
                svg_text(
                    plot_rect_left_x - 30,
                    0.5 * (last_rect_bottom + y_gap),
                    "Number of samples",
                    v_center=True,
                    h_center=True,
                    colour="dimgrey",
                    vertical=True,
                ),
                file=f,
            )

            print("</svg>", file=f)

        svg_export(svg_out, os.path.join(outdir, "plot.pdf"))
        svg_export(svg_out, os.path.join(outdir, "plot.png"))
