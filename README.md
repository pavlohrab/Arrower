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

Coloring is optional. With **no** `-a` file, all arrows are grey and no legend is
drawn. Pass a `-a` annotation file to color arrows per category and draw a legend
automatically; each arrow also shows `gene (category)` on hover.

```
python ArrowerSVG.py cluster.gbk -a annotations.tsv -o output.svg
```

The annotation file is a TSV/CSV (delimiter auto-detected) with **2 or 3 columns**:
`gene_id, category[, color]`. Blank lines and lines starting with `#` are ignored.

```
geneA	Kinase	    #e31a1c
geneB	Transporter
geneC	Kinase	    #e31a1c
geneD	Hydrolase
```

- **Column 1** (`gene_id`) is matched against the GenBank qualifier chosen by
  `--id-source` (see below).
- **Column 2** (`category`) groups genes; equal categories share a color.
- **Column 3** (`color`, optional) sets an explicit color for that category —
  a hex value like `#e31a1c` or an SVG color name like `red`. Categories without
  an explicit color get one from the built-in palette. (The first color seen for
  a category wins.)

Genes absent from the file are drawn grey. The canvas size adapts to the cluster
length and legend.

### Choosing the gene identifier

`--id-source` controls which GenBank CDS qualifier identifies a gene — this is
both the arrow label and the key matched against column 1 of the annotation file:

- `auto` (default) — tries `gene`, then `locus_tag`, then `protein_id`.
- `gene`, `locus_tag`, `protein_id` — force that specific qualifier.

`locus_tag` is often the most reliable choice, since it is usually present and
unique for every CDS:

```
python ArrowerSVG.py cluster.gbk -a annotations.tsv --id-source locus_tag
```

## Options

| Flag | Long form         | Default | Description                                      |
|------|-------------------|---------|--------------------------------------------------|
| `-o` | `--output`        | stdout  | Write the SVG to this file                       |
| `-a` | `--annotations`   | none    | TSV/CSV: `gene_id, category[, color]` (enables coloring) |
|      | `--id-source`     | auto    | GenBank qualifier to identify genes: `auto`/`gene`/`locus_tag`/`protein_id` |
| `-H` | `--arrow-height`  | 20      | Width of the arrow central part (px)             |
| `-E` | `--head-edge`     | 8       | Additional width of the arrow head (px)          |
| `-l` | `--head-length`   | 10      | Arrow head length (px)                           |
| `-X` | `--margin-x`      | 100     | Left-side margin (px)                            |
| `-Y` | `--margin-y`      | 30      | Top-side margin (px)                             |
| `-S` | `--scaling`       | 100.0   | Scaling of bp per px (100 = 100 bp/px)           |
| `-F` | `--font-size`     | 14      | Gene annotation font size                        |

Run `python ArrowerSVG.py --help` for the full list.
