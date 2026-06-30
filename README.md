# Arrower

Plot arrows for a gene cluster from a GenBank file and render them as SVG.

## Requirements

```
pip install -r requirements.txt
```

(Requires Python 3 and Biopython.)

## Usage

```
python ArrowerSVG.py GENBANK_FILE [options]
```

The SVG is printed to stdout by default. Use `-o` to write a file directly:

```
python ArrowerSVG.py cluster.gbk -o output.svg
```

The output SVG can then be edited in e.g. Adobe Illustrator.

## Coloring by category (e.g. Pfam domains)

Pass a `-a` annotation file mapping each gene to a category. Arrows are colored
per category and a legend is drawn automatically; each arrow also shows
`gene (category)` on hover.

```
python ArrowerSVG.py cluster.gbk -a annotations.tsv -o output.svg
```

The annotation file is a two-column TSV/CSV (delimiter auto-detected). Blank
lines and lines starting with `#` are ignored:

```
geneA	Kinase
geneB	Transporter
geneC	Kinase
geneD	Hydrolase
```

Genes absent from the file are drawn grey. With no `-a` file, all arrows are grey
and no legend is drawn. The canvas size adapts to the cluster length and legend.

## Options

| Flag | Long form         | Default | Description                                      |
|------|-------------------|---------|--------------------------------------------------|
| `-o` | `--output`        | stdout  | Write the SVG to this file                       |
| `-a` | `--annotations`   | none    | TSV/CSV mapping gene id → category (for coloring) |
| `-H` | `--arrow-height`  | 20      | Width of the arrow central part (px)             |
| `-E` | `--head-edge`     | 8       | Additional width of the arrow head (px)          |
| `-l` | `--head-length`   | 10      | Arrow head length (px)                           |
| `-X` | `--margin-x`      | 100     | Left-side margin (px)                            |
| `-Y` | `--margin-y`      | 30      | Top-side margin (px)                             |
| `-S` | `--scaling`       | 100.0   | Scaling of bp per px (100 = 100 bp/px)           |
| `-F` | `--font-size`     | 14      | Gene annotation font size                        |

Run `python ArrowerSVG.py --help` for the full list.
