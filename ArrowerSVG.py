#! /usr/bin/env python3

######################################################################
#                                                                    #
#           PLOT ARROWS FOR GENE CLUSTER GIVEN A GenBank FILE        #
#                           Peter Cimermancic                        #
#                               April 2010                           #
#                                                                    #
######################################################################

import argparse
import csv
from xml.sax.saxutils import escape

from Bio import SeqIO

# Default fill for genes that have no annotation / category.
DEFAULT_COLOR = "#cccccc"

# Qualitative palette (ColorBrewer "Paired", 12 colours). Categories beyond
# the palette length wrap around, so very domain-rich clusters may reuse hues.
PALETTE = [
    "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
    "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
    "#cab2d6", "#6a3d9a", "#ffff99", "#b15928",
]


# --- Draw arrow for gene
def arrow(X, Y, L, H, strand, h, l, color, title=""):
    '''
    SVG code for arrow:
        - (X,Y) ... upper left (+) or right (-) corner of the arrow
        - L ... arrow length
        - H ... arrow height
        - strand
        - h ... arrow head edge width
        - l ... arrow head length
        - color ... fill colour (hex string, e.g. "#1f78b4")
        - title ... hover text (gene name / category)
    the edges are ABCDEFG starting from (X,Y)
    '''

    if strand == '+':
        points = [
            [X, Y], [X + L - l, Y], [X + L - l, Y - h], [X + L, Y + H / 2],
            [X + L - l, Y + H + h], [X + L - l, Y + H], [X, Y + H],
        ]
        if L < l:
            # squize arrow if length shorter than head length
            points = [
                [X, Y], [X, Y], [X, Y - h], [X + L, Y + H / 2],
                [X, Y + H + h], [X, Y + H], [X, Y + H],
            ]

    elif strand == '-':
        points = [
            [X + L, Y], [X + l, Y], [X + l, Y - h], [X, Y + H / 2],
            [X + l, Y + H + h], [X + l, Y + H], [X + L, Y + H],
        ]
        if L < l:
            # squize arrow if length shorter than head length
            points = [
                [X + L, Y], [X + L, Y], [X + L, Y - h], [X, Y + H / 2],
                [X + L, Y + H + h], [X + L, Y + H], [X + L, Y + H],
            ]

    else:
        return ""

    points_str = " ".join("%i,%i" % (p[0], p[1]) for p in points)
    title_tag = "<title>%s</title>" % escape(title) if title else ""

    return """ <polygon points="%s"
                   style="fill:%s;fill-opacity:1.0;
                   stroke:#000000;stroke-width:2">%s</polygon>""" % (
        points_str, color, title_tag)


def line(X, Y, L):
    '''
    Draw a line below genes
    '''

    line = """<line x1="%i" y1="%i" x2="%i" y2="%i"
              style="stroke:rgb(99,99,99);stroke-width:2"/>""" % (X, Y, X + L, Y)
    return line


def text(text, X, Y, F):
    '''
    Write gene name
    '''
    line = """<text x="%i"  y="%i"
          style="font-family: Arial;
                 font-size  : %i;
                 font-style : italic;
                "
          >%s</text>""" % (X, Y, F, escape(text))
    return line


def legend(color_map, X, Y, font):
    '''
    Render a legend (colour swatch + category label) for the colour map.
    Returns the SVG fragment and the total height it occupies (in px).
    '''
    box = font
    gap = 6
    row_h = box + gap

    parts = []
    for i, category in enumerate(sorted(color_map)):
        color = color_map[category]
        y = Y + i * row_h
        parts.append(
            """<rect x="%i" y="%i" width="%i" height="%i"
               style="fill:%s;stroke:#000000;stroke-width:1"/>""" % (
                X, y, box, box, color))
        parts.append(
            """<text x="%i" y="%i"
               style="font-family: Arial; font-size: %i;"
               >%s</text>""" % (X + box + gap, y + box - 2, font, escape(category)))

    return "\n".join(parts), len(color_map) * row_h


def RGBToHTMLColor(rgb_tuple):
    '''
    Convert RGB color to HTML color format
    '''
    hexcolor = '#%02x%02x%02x' % rgb_tuple

    return hexcolor.strip('-')


# Qualifiers that 'auto' tries, in order, to identify a CDS feature.
ID_SOURCES = ('auto', 'gene', 'locus_tag', 'protein_id')


def feature_identifier(feature, id_source='auto'):
    '''
    Resolve the identifier of a CDS feature. This is both the label drawn under
    the arrow and the key matched against the annotation file.

    id_source selects which GenBank qualifier to use; 'auto' tries 'gene',
    then 'locus_tag', then 'protein_id'. Returns 'unknown' if none is present.
    '''
    keys = ('gene', 'locus_tag', 'protein_id') if id_source == 'auto' else (id_source,)
    for key in keys:
        if key in feature.qualifiers:
            return feature.qualifiers[key][0]
    return 'unknown'


def load_annotations(path):
    '''
    Load a gene -> category mapping (and optional category -> colour map) from
    a TSV/CSV file.

    Format: 2 or 3 columns, "gene_id <delimiter> category [<delimiter> colour]".
    The delimiter (tab, comma or semicolon) is auto-detected. The optional 3rd
    column sets an explicit colour for that category (hex like "#ff0000" or an
    SVG colour name like "red"); the first colour seen for a category wins.
    Blank lines and lines starting with '#' (e.g. a header) are ignored.

    Returns (mapping, category_colors).
    '''
    mapping = {}
    category_colors = {}
    with open(path, newline='') as fh:
        sample = fh.read(2048)
        fh.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters="\t,;")
        except csv.Error:
            dialect = csv.excel_tab
        for row in csv.reader(fh, dialect):
            if not row or row[0].lstrip().startswith('#'):
                continue
            if len(row) < 2:
                continue
            gene, category = row[0].strip(), row[1].strip()
            if not (gene and category):
                continue
            mapping[gene] = category
            color = row[2].strip() if len(row) >= 3 else ''
            if color and category not in category_colors:
                category_colors[category] = color
    return mapping, category_colors


def build_color_map(categories, overrides=None):
    '''
    Assign a colour to each category. Categories listed in `overrides` keep
    their explicit colour; the rest take the next palette colour. Order is
    preserved so that the same input always yields the same colours.
    '''
    overrides = overrides or {}
    color_map = {}
    palette_i = 0
    for category in categories:
        if category in color_map:
            continue
        if category in overrides:
            color_map[category] = overrides[category]
        else:
            color_map[category] = PALETTE[palette_i % len(PALETTE)]
            palette_i += 1
    return color_map


def SVG(GenBankFile, annotations=None, category_colors=None, id_source='auto',
        ArrowHeight=20, HeadEdge=8, HeadLength=10,
        marginX=100, marginY=30, scaling=100.0, font=14):
    '''
    Create the main SVG document:
        - read in GenBank document
        - find genes, start and stop positions, and strands
        - colour arrows by category (from the annotation mapping); when no
          annotations are given every arrow is grey and no legend is drawn
        - draw a legend and return the SVG as a string
    '''
    annotations = annotations or {}

    # --- read in GenBank file
    with open(GenBankFile, 'r') as handle:
        seq_records = list(SeqIO.parse(handle, "genbank"))

    # --- collect the categories actually present, in genomic order, and map
    #     them to colours
    ordered_categories = []
    for seq_record in seq_records:
        for feature in seq_record.features:
            if feature.type == 'CDS':
                category = annotations.get(feature_identifier(feature, id_source))
                if category and category not in ordered_categories:
                    ordered_categories.append(category)
    color_map = build_color_map(ordered_categories, category_colors)

    # --- build the body, tracking the extent so the canvas can be sized
    body = ""
    max_right = marginX
    content_bottom = marginY + ArrowHeight

    for seq_record in seq_records:

        # draw a line that corresponds to cluster size
        ClusterSize = len(seq_record.seq)
        cluster_px = ClusterSize / scaling
        max_right = max(max_right, marginX + cluster_px)
        body += line(marginX, marginY + ArrowHeight / 2, cluster_px)

        previousStop = 0  # keep track of previous stop to adjust text location
        previousLevel = 1
        for feature in seq_record.features:
            if feature.type == 'CDS':

                GeneName = feature_identifier(feature, id_source)
                category = annotations.get(GeneName)
                color = color_map.get(category, DEFAULT_COLOR)

                strand = feature.location.strand
                if strand == -1:
                    strand = '-'
                elif strand == 1:
                    strand = '+'

                start = str(feature.location.start)
                if '>' in start or '<' in start:
                    start = start.replace('>', '').replace('<', '')
                start = int(start) / scaling

                stop = str(feature.location.end)
                if '>' in stop or '<' in stop:
                    stop = stop.replace('>', '').replace('<', '')
                stop = int(stop) / scaling

                # write arrow to SVG file
                title = GeneName if category is None else "%s (%s)" % (GeneName, category)
                body += arrow(start + marginX, marginY, stop - start, ArrowHeight,
                              strand, HeadEdge, HeadLength, color, title=title)

                # annotate genes
                if previousLevel == 4:
                    previousLevel = 0
                if previousStop == 0:
                    text_y = marginY + 50
                else:
                    if ((start + (stop - start) / 2 - len(GeneName) * font / 4.) - previousStop) < 10:
                        text_y = marginY + 50 + (font + 6) * previousLevel
                        previousLevel += 1
                    else:
                        previousLevel = 0
                        text_y = marginY + 50 + (font + 6) * previousLevel
                text_x = marginX + start + (stop - start) / 2 - len(GeneName) * font / 4.
                body += text(GeneName, text_x, text_y, font)
                content_bottom = max(content_bottom, text_y)
                previousStop = stop

    # --- legend (only when there is something to show)
    if color_map:
        legend_y = content_bottom + 2 * font
        legend_svg, legend_h = legend(color_map, marginX, legend_y, font)
        body += "\n" + legend_svg
        content_bottom = legend_y + legend_h

    # --- size the canvas to the content
    width = int(max_right + marginX)
    height = int(content_bottom + marginY)

    header = """<?xml version="1.0" standalone="no"?>
    <!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"
    "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">

    <svg width="%i" height="%i" xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink'>
    """ % (width, height)

    return header + body + '\n</svg>'


def parse_args(argv=None):
    '''
    Parse command line arguments.
    '''
    parser = argparse.ArgumentParser(
        description="Plot arrows for a gene cluster given a GenBank file, "
                    "producing an SVG that can be further edited in e.g. Adobe Illustrator.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "genbank",
        metavar="GENBANK_FILE",
        help="input GenBank file describing the gene cluster",
    )
    parser.add_argument(
        "-o", "--output",
        metavar="SVG_FILE",
        default=None,
        help="write the SVG to this file (default: print to stdout)",
    )
    parser.add_argument(
        "-a", "--annotations",
        metavar="TSV_FILE",
        default=None,
        help="optional TSV/CSV file with columns 'gene_id, category[, colour]'. "
             "Arrows are coloured by category and a legend is drawn. The optional "
             "3rd column sets an explicit colour for that category (e.g. '#ff0000' "
             "or 'red'), otherwise a palette colour is used. gene_id is matched "
             "against the qualifier chosen by --id-source. Without this option all "
             "arrows are grey and no legend is drawn.",
    )
    parser.add_argument(
        "--id-source",
        choices=ID_SOURCES,
        default="auto",
        help="which GenBank CDS qualifier identifies a gene, used both as the "
             "arrow label and to match the annotation file. 'auto' tries 'gene', "
             "then 'locus_tag', then 'protein_id'",
    )
    parser.add_argument(
        "-H", "--arrow-height",
        type=int, default=20,
        help="width of the arrow central part, in pixels",
    )
    parser.add_argument(
        "-E", "--head-edge",
        type=int, default=8,
        help="additional width of the arrow head, in pixels",
    )
    parser.add_argument(
        "-l", "--head-length",
        type=int, default=10,
        help="arrow head length, in pixels",
    )
    parser.add_argument(
        "-X", "--margin-x",
        type=int, default=100,
        help="left-side margin, in pixels",
    )
    parser.add_argument(
        "-Y", "--margin-y",
        type=int, default=30,
        help="top-side margin, in pixels",
    )
    parser.add_argument(
        "-S", "--scaling",
        type=float, default=100.0,
        help="scaling of bp per px (100 means 100 bp/px)",
    )
    parser.add_argument(
        "-F", "--font-size",
        type=int, default=14,
        help="gene annotation font size",
    )

    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    if args.annotations:
        annotations, category_colors = load_annotations(args.annotations)
    else:
        annotations, category_colors = {}, {}

    svg = SVG(
        args.genbank,
        annotations=annotations,
        category_colors=category_colors,
        id_source=args.id_source,
        ArrowHeight=args.arrow_height,
        HeadEdge=args.head_edge,
        HeadLength=args.head_length,
        marginX=args.margin_x,
        marginY=args.margin_y,
        scaling=args.scaling,
        font=args.font_size,
    )

    if args.output:
        with open(args.output, 'w') as out:
            out.write(svg)
    else:
        print(svg)


if __name__ == '__main__':
    main()
