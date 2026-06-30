#! /usr/bin/env python3

######################################################################
#                                                                    #
#           PLOT ARROWS FOR GENE CLUSTER GIVEN A GenBank FILE        #
#                           Peter Cimermancic                        #
#                               April 2010                           #
#                                                                    #
######################################################################

import argparse

from Bio import SeqIO

# --- Draw arrow for gene
def arrow(X,Y,L,H,strand,h,l,color):
    '''
    SVG code for arrow:
        - (X,Y) ... upper left (+) or right (-) corner of the arrow
        - L ... arrow length
        - H ... arrow height
        - strand
        - h ... arrow head edge width
        - l ... arrow head length
        - color
        - strand
    the edges are ABCDEFG starting from (X,Y)     
    '''
    
    if strand == '+':
        
        A = [X,Y]
        B = [X+L-l,Y]
        C = [X+L-l,Y-h]
        D = [X+L,Y+H/2]
        E = [X+L-l,Y+H+h]
        F = [X+L-l,Y+H]
        G = [X,Y+H]

        if L < l:
            # squize arrow if length shorter than head length
            B = [X,Y]
            C = [X,Y-h]
            D = [X+L,Y+H/2]
            E = [X,Y+H+h]
            F = [X,Y+H]

        line = """ <polygon id="bla" points="%i,%i %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i"
                   style="fill:#cccccc;fill-opacity:1.0;
                   stroke:#000000;stroke-width:2">
                   </polygon>""" % (A[0],A[1],B[0],B[1],C[0],C[1],D[0],D[1],E[0],E[1],F[0],F[1],G[0],G[1])
        
        return line

    elif strand == '-':
        
        A = [X+L,Y]
        B = [X+l,Y]
        C = [X+l,Y-h]
        D = [X,Y+H/2]
        E = [X+l,Y+H+h]
        F = [X+l,Y+H]
        G = [X+L,Y+H]

        if L < l:
            # squize arrow if length shorter than head length
            B = [X+L,Y]
            C = [X+L,Y-h]
            D = [X,Y+H/2]
            E = [X+L,Y+H+h]
            F = [X+L,Y+H]

        line = """ <polygon id="bla" points="%i,%i %i,%i %i,%i %i,%i %i,%i %i,%i %i,%i"
                   style="fill:#cccccc;fill-opacity:1.0;
                   stroke:#000000;stroke-width:2">
                   </polygon>""" % (A[0],A[1],B[0],B[1],C[0],C[1],D[0],D[1],E[0],E[1],F[0],F[1],G[0],G[1])
        
        return line
    
    else: return 0

def line(X,Y,L):
    '''
    Draw a line below genes
    '''
    
    line = """<line x1="%i" y1="%i" x2="%i" y2="%i"
              style="stroke:rgb(99,99,99);stroke-width:2"/>""" % (X,Y,X+L,Y)
    return line


def text(text,X,Y,F):
    '''
    Write gene name
    '''
    line = """<text x="%i"  y="%i"
          style="font-family: Arial;
                 font-size  : %i;
                 font-style : italic;
                "
          >%s</text>""" % (X,Y,F,text)
    return line
    

def RGBToHTMLColor(rgb_tuple):
    '''
    Convert RGB color to HTML color format
    '''
    hexcolor = '#%02x%02x%02x' % rgb_tuple
    
    return hexcolor.strip('-')


def SVG(GenBankFile,ArrowHeight=20,HeadEdge=8,HeadLength=10,marginX=100,marginY=30,scaling=100.0,font=14):
    '''
    Create the main SVG document:
        - read in GenBank documnet
        - find genes, start and stop positions, and strands
        - write the SVG files
    '''
    
    # --- create SVG header
    ALL_TEXT = ""

    header = """<?xml version="1.0" standalone="no"?>
    <!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" 
    "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">

    <svg width="870" height="300" xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink'
    onload='Init(evt)'>
    """ #% Width

    ALL_TEXT += header

    # --- read in GenBank file
    with open(GenBankFile, 'r') as handle:
        seq_records = list(SeqIO.parse(handle, "genbank"))

    for seq_record in seq_records:
        
        # draw a line that corespond to cluster size
        ClusterSize = len(seq_record.seq)
        ALL_TEXT += line(marginX,marginY+ArrowHeight/2,ClusterSize / scaling)
        
        previousStop = 0 # keep track of previous stop to adjust text location
        previousLevel = 1
        for feature in seq_record.features:
            if feature.type == 'CDS':
                
                GeneName = feature.qualifiers['gene'][0]

                strand = feature.strand
                if strand == -1:
                    strand = '-'
                elif strand == 1:
                    strand = '+'

                start = str(feature.location.start)
                if '>' in start or '<' in start:
                    start = start.replace('>','').replace('<','')
                start = int(start) / scaling
                
                stop = str(feature.location.end) 
                if '>' in stop or '<' in stop:
                    stop = stop.replace('>','').replace('<','')
                stop = int(stop) / scaling
                
                # write arrow to SVG file
                ALL_TEXT += arrow(start+marginX,marginY,stop-start,ArrowHeight,strand,HeadEdge,HeadLength,GeneName)

                # annotate genes
                if previousLevel  == 4: previousLevel = 0
                if previousStop == 0:
                    ALL_TEXT += text(GeneName,marginX+start+(stop-start)/2-len(GeneName)*font/4.,marginY+50,font)
                else:
                    if ((start+(stop-start)/2-len(GeneName)*font/4.) - previousStop) < 10:
                        ALL_TEXT += text(GeneName,marginX+start+(stop-start)/2-len(GeneName)*font/4.,marginY+50+(font+6)*previousLevel,font)
                        previousLevel += 1
                    else:
                        previousLevel = 0
                        ALL_TEXT += text(GeneName,marginX+start+(stop-start)/2-len(GeneName)*font/4.,marginY+50+(font+6)*previousLevel,font)
                previousStop = stop


    ALL_TEXT += '</svg>'

    return ALL_TEXT

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

    svg = SVG(
        args.genbank,
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
