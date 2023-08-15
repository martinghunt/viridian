import csv
from enum import Enum
import gzip
import logging
import math
import multiprocessing
from operator import itemgetter
import os

from viridian import amplicon_schemes, constants, scheme_id, utils


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
    "ACGT_GOOD": "#40826D",
    "ACGT_DP": "pink",
    "ACGT_BAD": "red",
    "INDEL": "lemonchiffon",
    "HET_GOOD": "lightgreen",
    "HET_DP": "lightblue",
    "HET_BAD": "cornflowerblue",
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


COVID_GENES = [
    {"start": 0, "end": 264, "name": "5'UTR"},
    {"start": 265, "end": 21554, "name": "ORF1ab"},
    {"start": 21562, "end": 25383, "name": "S"},
    {"start": 25392, "end": 26219, "name": "ORF3a"},
    {"start": 26244, "end": 26471, "name": "E"},
    {"start": 26522, "end": 27190, "name": "M"},
    {"start": 27201, "end": 27386, "name": "ORF6"},
    {"start": 27393, "end": 27758, "name": "ORF7a"},
    {"start": 27893, "end": 28258, "name": "ORF8"},
    {"start": 28273, "end": 29532, "name": "N"},
    {"start": 29557, "end": 29673, "name": "ORF10"},
    {"start": 29674, "end": 29902, "name": "3'UTR"},
]


def svg_export(infile, outfile):
    if outfile.endswith(".png"):
        opts = "-b 'rgb(255, 255, 255)' -d 400"
    else:
        opts = ""
    utils.syscall(f"inkscape {opts} {infile} -o {outfile}")


def svg_header_lines(width, height):
    return [
        '<?xml version="1.0" standalone="no"?>',
        '<!DOCTYPE svg PUBLIC " -//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">',
        f'<svg width="{width}" height="{height}">',
    ]


def svg_line(x1, y1, x2, y2, colour, stroke_width=1, linecap="butt", dasharray=None):
    dasharray = "" if dasharray is None else f'stroke-dasharray="{dasharray}"'
    return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{colour}" stroke-width="{stroke_width}" stroke-linecap="{linecap}" {dasharray} />'


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
    axis_colour="black",
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


def grid_v_lines(x_positions, top, bottom, x_min, x_max):
    lines = []
    step = x_positions[1] - x_positions[0]
    opt_big = {"colour": "cornflowerblue", "stroke_width": 0.75, "dasharray": "8 2"}
    opt_small = {"colour": "cornflowerblue", "stroke_width": 0.5, "dasharray": "5 2"}
    first_small = x_positions[0] - 0.5 * step
    if x_min <= first_small:
        lines.append(svg_line(first_small, top, first_small, bottom, **opt_small))

    for i, x in enumerate(x_positions):
        lines.append(svg_line(x, top, x, bottom, **opt_big))
        small_x = x + 0.5 * step
        if small_x <= x_max:
            lines.append(svg_line(small_x, top, small_x, bottom, **opt_small))
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


def amp_primer_track(
    genome_start,
    genome_end,
    x_left,
    x_right,
    y_top,
    y_bottom,
    scheme_name,
    plot_primers=True,
    plot_amp_names=True,
    plot_amp_number=False,
):
    logging.info(f"Making amplicon/primer track for scheme {scheme_name}")
    schemes = amplicon_schemes.get_built_in_schemes()
    assert scheme_name in schemes
    scheme = scheme_id.Scheme(tsv_file=schemes[scheme_name])
    logging.info(f"Loaded scheme info {scheme_name}")
    lines = []
    at_top = False
    x_scale = (x_right - x_left) / (genome_end - genome_start)
    x_trans = lambda x: x_left + (x - genome_start) * x_scale
    text_opts = {
        "h_center": True,
        "v_center": True,
        "font_size": 8,
        "colour": "dimgrey",
    }
    primer_opts = {"stroke_width": 2, "linecap": "round"}

    for amp_number, amp in enumerate(sorted(scheme.amplicons, key=itemgetter("start"))):
        if genome_start > amp["end"] or amp["start"] > genome_end:
            continue
        at_top = not at_top
        x_start = max(x_left, x_trans(amp["start"]))
        x_end = min(x_right, x_trans(amp["end"]))
        y_middle = 0.5 * (y_top + y_bottom)
        y = y_top if at_top else y_bottom
        lines.append(svg_line(x_start, y, x_end, y, "dimgrey", stroke_width=5))

        if plot_amp_number:
            lines.append(
                svg_text(0.5 * (x_start + x_end), y_middle, amp_number + 1, **text_opts)
            )
        elif plot_amp_names:
            lines.append(
                svg_text(0.5 * (x_start + x_end), y_middle, amp["name"], **text_opts)
            )

        if not plot_primers:
            continue

        for l_r in amp["primers"]:
            y = y_top + 4 if at_top else y_bottom - 4
            previous_start = previous_end = None
            for start, end in sorted(amp["primers"][l_r], reverse=l_r == "right"):
                if not genome_start <= start <= end <= genome_end:
                    continue

                if (
                    previous_start is not None
                    and previous_end is not None
                    and previous_start <= end
                    and start <= previous_end
                ):
                    if at_top:
                        y += 5
                    else:
                        y -= 5
                previous_start = start
                previous_end = end
                y_add = 2 if at_top else -2
                s = x_trans(start) + 1
                e = x_trans(end) - 1
                lines.append(svg_line(s, y, e, y, "red", **primer_opts))

                if e - s < 3:
                    continue

                if l_r == "right":
                    lines.append(svg_line(s, y, s + 3, y + y_add, "red", **primer_opts))
                else:
                    lines.append(svg_line(e, y, e - 3, y + y_add, "red", **primer_opts))

    return lines


def genes_track(
    genome_start,
    genome_end,
    x_left,
    x_right,
    y_top,
    y_bottom,
):
    logging.info(f"Making genes track")
    print("x_left, x_right", x_left, x_right)
    x_scale = (x_right - x_left) / (genome_end - genome_start)
    x_trans = lambda x: x_left + (x - genome_start) * x_scale
    y = [y_top, y_top, y_bottom, y_bottom]
    lines = []
    font_size = 10
    label_above = True

    for d in COVID_GENES:
        if genome_start > d["end"] or d["start"] > genome_end:
            continue
        s = x_trans(max(d["start"], genome_start))
        e = x_trans(min(d["end"], genome_end))
        x_middle = 0.5 * (s + e)
        x_coords = [s, e, e, s]

        lines.append(
            svg_polygon(
                "mintcream",
                "dimgrey",
                x_coords=x_coords,
                y_coords=y,
                border_width=1,
            )
        )

        if 7 * len(d["name"]) > e - s:  # label too big for gene rectangle
            if d["name"] == "ORF8" and genome_end - genome_start > 20000:
                extra = font_size
            else:
                extra = 0

            if label_above:
                y_pos = y_top - font_size - extra
                line_y1 = y_top - 0.4 * font_size - extra
                line_y2 = y_top
            else:
                y_pos = y_bottom + font_size + extra
                line_y1 = y_bottom + 0.4 * font_size + extra
                line_y2 = y_bottom

            lines.append(
                svg_line(
                    x_middle, line_y1, x_middle, line_y2, "dimgrey", stroke_width=1
                )
            )
            label_above = not label_above
        else:
            y_pos = 0.5 * (y_top + y_bottom)
            label_above = True

        lines.append(
            svg_text(
                x_middle,
                y_pos,
                d["name"],
                colour="dimgrey",
                font_size=font_size,
                v_center=True,
                h_center=True,
            )
        )

    return lines


def legend(x_left, y_middle, colours, square_size=11, y_gap=5, font_size=11):
    total_height = (
        len(PLOT_BASECALL_ORDER) * square_size + (len(PLOT_BASECALL_ORDER) - 1) * y_gap
    )
    y_top = y_middle - 0.5 * total_height
    y = y_top
    x_r = x_left + square_size
    text_x = x_r + 3
    x_coords = [x_left, x_left, x_r, x_r]
    text_opts = {"font_size": font_size, "colour": "black", "v_center": True}
    lines = []
    for key in reversed(PLOT_BASECALL_ORDER):
        lines.append(svg_text(text_x, y + 0.5 * square_size, key, **text_opts))
        y_coords = [y, y + square_size, y + square_size, y]
        lines.append(
            svg_polygon(colours[key], "black", x_coords=x_coords, y_coords=y_coords)
        )
        y += square_size + y_gap
    return lines

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

    def add_other_genome_calls(self, other):
        assert len(self) == len(other)
        for this_list, other_list in zip(self.calls, other.calls):
            for i in range(len(this_list)):
                this_list[i] += other_list[i]

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
        of = gzip.open if file_of_qc_filenames.endswith(".gz") else open
        with of(file_of_qc_filenames, "rt") as f:
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

    def from_index_file(self, index_file, cpus=1):
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

    def to_pickle(self, outfile, cpus=1):
        utils.write_pickle(outfile, self.genome_calls)
        logging.info(f"Data written to file {outfile}")

    def from_pickle(self, infile):
        self.genome_calls = utils.load_pickle(infile)
        genome_lengths = set(len(x) for x in self.genome_calls.values())
        assert len(genome_lengths) == 1
        self.ref_length = genome_lengths.pop()
        logging.info(f"Data loaded from file {infile}")

    def add_other_plot(self, other):
        for other_name, other_calls in other.genome_calls.items():
            if other_name not in self.genome_calls:
                logging.info(f"Adding new data set {other_name}")
                self.add_genome_call_set(other_name, other_calls)
            else:
                logging.info(f"Updating existing data set {other_name}")
                self.genome_calls[other_name].add_other_genome_calls(other_calls)

    def combine_pickles(self, pickles_fofn, outfile):
        logging.info(f"Combining all data from files in {pickles_fofn}")
        with open(pickles_fofn) as f:
            for filename in map(str.rstrip, f):
                logging.info(f"Adding data from file {filename}")
                new_plots = Plots()
                new_plots.from_pickle(filename)
                self.add_other_plot(new_plots)
        logging.info(f"Finished combining data from {pickles_fofn}")
        self.to_pickle(outfile)

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
        dataset_height=100,
        genes_track_height=10,
        amp_track_height=20,
        y_gap=50,
        x_tick_step=None,
        y_tick_step=None,
        amp_scheme=None,
        plot_amp_names=False,
        plot_amp_number=False,
        plot_primers=True,
        plot_genes=False,
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
        logging.info(f"Making plots for genome coords {x_start+1}-{x_end+1}")
        colours = {
            Basecall[k].name: colours.get(k, v) for k, v in PLOT_DEFAULT_COLOURS.items()
        }
        plot_rect_left_x = 100
        plot_rect_right_x = plot_width - 100
        rect_width = plot_rect_right_x - plot_rect_left_x
        rect_y_top = y_gap
        grid_top = rect_y_top
        rect_y_bottom = rect_y_top + dataset_height
        rect_x_scale = rect_width / (x_end - x_start)
        rect_middle = 0.5 * (plot_rect_left_x + plot_rect_right_x)
        x_ticks_abs = tick_positions_from_range(x_start, x_end, x_tick_step)
        x_ticks_pos = [
            plot_rect_left_x + (x - x_start) * rect_x_scale for x in x_ticks_abs
        ]

        os.mkdir(outdir)
        svg_lines = []
        x_coords = [
            plot_rect_left_x + (x - x_start) * rect_x_scale
            for x in [x_start] + list(range(x_start, x_end + 1)) + [x_end, x_start]
        ]

        for set_name, calls in self.genome_calls.items():
            svg_lines.append(
                svg_text(rect_middle, rect_y_top - 8, set_name, h_center=True)
            )
            y_vals = calls.get_plot_y_vals(x_start=x_start, x_end=x_end)
            max_y = max(y_vals[Basecall[PLOT_BASECALL_ORDER[-1]].value])

            for base_type in reversed(PLOT_BASECALL_ORDER):
                y_coords = [rect_y_bottom] + [
                    rect_y_bottom - (y * dataset_height / max_y)
                    for y in y_vals[Basecall[base_type].value]
                ]
                y_coords.extend([rect_y_bottom, rect_y_bottom])
                svg_lines.append(
                    svg_polygon(
                        colours[base_type],
                        colours[base_type],
                        x_coords=x_coords,
                        y_coords=y_coords,
                        border_width=0,
                    )
                )

            y_ticks_abs = tick_positions_from_range(0, max_y, y_tick_step)
            y_ticks_pos = [
                rect_y_bottom - (y * dataset_height / max_y) for y in y_ticks_abs
            ]
            svg_lines.extend(
                y_axis(
                    rect_y_top,
                    rect_y_bottom,
                    plot_rect_left_x,
                    y_ticks_abs,
                    y_ticks_pos,
                )
            )
            rect_y_top = rect_y_bottom + y_gap
            rect_y_bottom = rect_y_top + dataset_height
            logging.info(f"Plotted data set '{set_name}'")

        last_rect_bottom = rect_y_top - y_gap
        grid_bottom = last_rect_bottom
        svg_lines.append(
            svg_text(
                plot_rect_left_x - 30,
                0.5 * (last_rect_bottom + y_gap),
                "Number of samples",
                v_center=True,
                h_center=True,
                colour="black",
                vertical=True,
            )
        )

        if amp_scheme is not None:
            primer_top = last_rect_bottom + 0.5 * y_gap
            primer_bottom = primer_top + amp_track_height
            primer_middle = 0.5 * (primer_top + primer_bottom)
            svg_lines.extend(
                amp_primer_track(
                    x_start,
                    x_end,
                    plot_rect_left_x,
                    plot_rect_right_x,
                    primer_top,
                    primer_bottom,
                    amp_scheme,
                    plot_primers=plot_primers,
                    plot_amp_names=plot_amp_names,
                    plot_amp_number=plot_amp_number,
                )
            )
            svg_lines.append(
                svg_text(
                    plot_rect_left_x - 7,
                    primer_middle,
                    "Amplicons",
                    colour="black",
                    h_right=True,
                    v_center=True,
                )
            )
            grid_bottom = primer_bottom + y_gap

        if plot_genes:
            genes_top = primer_bottom + 0.75 * y_gap
            genes_bottom = genes_top + genes_track_height
            svg_lines.extend(
                genes_track(
                    x_start,
                    x_end,
                    plot_rect_left_x,
                    plot_rect_right_x,
                    genes_top,
                    genes_bottom,
                )
            )
            svg_lines.append(
                svg_text(
                    plot_rect_left_x - 7,
                    0.5 * (genes_top + genes_bottom),
                    "Genes",
                    colour="black",
                    h_right=True,
                    v_center=True,
                )
            )
            grid_bottom = genes_bottom + 30

        svg_lines.extend(
            x_axis(
                plot_rect_left_x,
                plot_rect_right_x,
                grid_bottom,
                x_ticks_abs,
                x_ticks_pos,
            )
        )
        svg_lines.append(
            svg_text(
                rect_middle,
                grid_bottom + 25,
                "Position in genome",
                h_center=True,
                v_center=True,
                colour="black",
            )
        )

        v_lines = grid_v_lines(
            x_ticks_pos, grid_top - 4, grid_bottom, plot_rect_left_x, plot_rect_right_x
        )
        svg_out = os.path.join(outdir, "plot.svg")
        with open(svg_out, "w") as f:
            print(*svg_header_lines(plot_width, grid_bottom + 40), sep="\n", file=f)
            print(*v_lines, sep="\n", file=f)
            print(*svg_lines, sep="\n", file=f)
            print(
                *legend(
                    plot_rect_right_x + 10, 0.5 * (grid_top + last_rect_bottom), colours
                ),
                sep="\n",
                file=f,
            )
            print("</svg>", file=f)

        svg_export(svg_out, os.path.join(outdir, "plot.pdf"))
        svg_export(svg_out, os.path.join(outdir, "plot.png"))
