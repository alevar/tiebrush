#!/usr/bin/env python

import numpy as np
import argparse
import sys
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.gridspec import GridSpec
from adjustText import adjust_text
mpl.use('Agg')


class TX:
    def __init__(self):
        self.seqid = None
        self.strand = None
        self.exons = []
        self.orf = []

        self.tid = None

    def parse_from_gtf(self, gtf_lines):
        for line in gtf_lines.split("\n"):
            lcs = line.strip().split("\t")
            if not len(lcs) == 9:
                continue

            assert self.seqid is None or self.seqid == lcs[0], "variable values in 1st column (seqid)"
            self.seqid = lcs[0]
            assert self.strand is None or self.strand == lcs[6], "variable values in 7th column (strand)"
            self.strand = lcs[6]

            txid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]
            if self.tid is None:
                self.tid = txid
            assert self.tid == txid, "transcript IDs do not match. function accepts a single transcript entry only"

            if lcs[2].lower() == "transcript" or lcs[2].lower() == "mrna":
                continue

            if lcs[2].lower() == "exon":
                self.exons.append((int(lcs[3]), int(lcs[4])))
            if lcs[2].lower() == "cds":
                self.orf.append((int(lcs[3]), int(lcs[4])))

        # sort exons and orf
        self.exons.sort(key=lambda l: l[0])
        self.orf.sort(key=lambda l: l[0])

        assert self.tid is not None, "tid not set"
        assert len(self.exons) > 0, "exon chain is empty"

    def nume(self):
        return len(self.exons)

    def numc(self):
        return len(self.orf)

    def get_introns(self):
        if len(self.exons) > 1:
            for i in range(len(self.exons) - 1):
                yield self.exons[i][1], self.exons[i + 1][0]

    def get_exons(self):
        for e in self.exons:
            yield e

    def get_cds(self):
        for c in self.orf:
            yield c

    def get_tid(self):
        return self.tid

    def get_strand(self):
        return self.strand

    def get_seqid(self):
        return self.seqid

    def get_start(self):
        return self.exons[0][0]

    def get_end(self):
        return self.exons[-1][1]

    def get_cstart(self):
        return self.orf[0][0]

    def get_cend(self):
        return self.orf[-1][1]

    def print_contents(self):
        print(self.tid)
        print(self.seqid + self.strand + ":" + str(self.get_start()) + "-" + str(self.get_end()))
        print(self.exons)
        print(self.orf)


class Locus:
    def __init__(self):
        self.txs = list()
        self.seqid = None
        self.strand = None

        self.intervals = []  # union of all exons in the locus (minus the introns)
        self.introns = dict()
        self.intron_cov_lst = list()

        self.exon_starts = []
        self.exon_ends = []

        self.graphcoords = None
        self.graphToGene = None
        self.covx_lst = list()
        self.cov_lst = list()

        self.cov_full_lst = list()

        self.settings = None

        self.num_cov_tracks = 0
        self.num_sj_tracks = 0

    @staticmethod
    def union(intervals):
        res = []
        for s, e in sorted(intervals):
            if res and res[-1][1] >= s - 1:
                res[-1][1] = max(res[-1][1], e)
            else:
                res.append([s, e])

        return [(x[0], x[1]) for x in res]

    @staticmethod
    def cubic_bezier(pts, t):
        """
        Get points in a cubic bezier.
        """
        p0, p1, p2, p3 = pts
        p0 = np.array(p0)
        p1 = np.array(p1)
        p2 = np.array(p2)
        p3 = np.array(p3)
        return p0 * (1 - t) ** 3 + 3 * t * p1 * (1 - t) ** 2 + \
               3 * t ** 2 * (1 - t) * p2 + t ** 3 * p3

    def add_tx(self, tx):
        assert self.seqid is None or self.seqid == tx.get_seqid(), "mismatching seqids"
        assert self.strand is None or self.strand == tx.get_strand(), "mismatching strands"

        self.seqid = tx.get_seqid()
        self.strand = tx.get_strand()

        self.intervals = Locus.union(self.intervals + tx.exons)

        intron = [-1, -1]
        for s, e in tx.exons:
            self.exon_starts.append(s)
            self.exon_ends.append(e)

            intron[1] = s
            if not intron[0] == -1:
                self.introns[tuple(intron)] = 0
            intron[0] = e

        self.txs.append(tx)

    def set_scaling(self):
        # get graphcoords
        if self.graphcoords is None:
            self.graphcoords, self.graphToGene = self.getScaling(self.settings["intron_scale"], self.settings["exon_scale"],
                                                                 self.settings["reverse_minus"])
            
    def get_start(self):
        return self.intervals[0][0]

    def get_end(self):
        return self.intervals[-1][1]

    def add_introns(self, sj_fname):
        assert os.path.exists(sj_fname),"Splice Junction track does not exist"

        self.intron_cov_lst.append(dict())

        with open(sj_fname, "r") as inFP:
            for line in inFP:
                lcs = line.strip().split("\t")
                if not len(lcs) == 6:
                    continue

                if not lcs[0] == self.seqid:
                    continue

                if not lcs[5] == self.strand:
                    continue

                intron = (int(lcs[1]), int(lcs[2]) + 1)
                if intron in self.introns:
                    self.intron_cov_lst[-1][intron] = int(lcs[4])

        self.num_sj_tracks+=1

    def getScaling(self, intron_scale, exon_scale, reverse_minus):
        """
        Compute the scaling factor across various genic regions.
        """

        tx_start = self.get_start()
        tx_end = self.get_end()

        exoncoords = np.zeros((tx_end - tx_start + 1))
        for i in range(len(self.exon_starts)):
            exoncoords[self.exon_starts[i] - tx_start: self.exon_ends[i] - tx_start] = 1

        graphToGene = {}
        graphcoords = np.zeros((tx_end - tx_start + 1), dtype='f')
        x = 0
        if self.strand == '+' or not reverse_minus:
            for i in range(tx_end - tx_start + 1):
                graphcoords[i] = x
                graphToGene[int(x)] = i + tx_start
                if exoncoords[i] == 1:
                    x += 1. / exon_scale
                else:
                    x += 1. / intron_scale
        else:
            for i in range(tx_end - tx_start + 1):
                graphcoords[-(i + 1)] = x
                graphToGene[int(x)] = tx_end - i + 1
                if exoncoords[-(i + 1)] == 1:
                    x += 1. / exon_scale
                else:
                    x += 1. / intron_scale
        return graphcoords, graphToGene

    def set_settings(self, settings):
        self.settings = settings

    def compress_intervals(self, vals, graphcoords):  # intervals with optional values if more
        compressed_x = []
        compressed_wiggle = []
        prevx = graphcoords[0]
        tmpval = []
        for i in range(len(graphcoords)):
            tmpval.append(vals[i])
            if abs(graphcoords[i] - prevx) > self.settings["resolution"]:
                compressed_wiggle.append(np.mean(tmpval))
                compressed_x.append(prevx)
                prevx = graphcoords[i]
                tmpval = []

        return compressed_x, compressed_wiggle

    def add_coverage(self, cov_fname):
        assert os.path.exists(cov_fname),"Coverage track does not exist: "+cov_fname
        assert os.path.exists(cov_fname),"Coverage track does not exist: "+cov_fname

        # use the graphcoords to perform interval compression below
        self.cov_full_lst.append(list())
        self.cov_lst.append(list())
        self.covx_lst.append(list())

        self.cov_full_lst[-1] = [0 for i in range(self.get_start(), self.get_end() + 1, 1)]
        with open(cov_fname, "r") as inFP:
            for line in inFP:
                lcs = line.strip().split("\t")
                if not len(lcs) == 4:
                    continue

                if not lcs[0] == self.seqid:
                    continue

                if int(lcs[2]) <= self.get_start() or int(lcs[1]) > self.get_end():
                    continue

                # process coverage
                for v in range(int(lcs[1]), min(self.get_end(),int(lcs[2])), 1):
                    self.cov_full_lst[-1][v - self.get_start()] = int(lcs[3])

        # compress the vals
        self.covx_lst[-1], self.cov_lst[-1] = self.compress_intervals(self.cov_full_lst[-1], self.graphcoords)
        self.num_cov_tracks+=1

    def get_coords_str(self):
        return self.seqid + self.strand + ":" + str(self.get_start()) + "-" + str(self.get_end())

    def plot(self,out_fname):
        assert self.num_sj_tracks==self.num_cov_tracks or self.num_sj_tracks==0,"incompatible number of splice junciton and coverage tracks - the numebrs should either be the same or no splice junction tracks provided at all"

        color_dens = "#ffb703"
        color_spines = "#fb8500"
        color_exon = "#023047"
        color_cds = "#219ebc"

        hrs = [4]*self.num_cov_tracks+[1 for i in range(len(self.txs))]
        gs1hs = 1
        gs2hs = 0.3

        fig = plt.figure(figsize=(self.settings["fig_width"],self.settings["fig_height"]))

        gs1 = GridSpec(len(self.txs)+self.num_cov_tracks, 1, height_ratios=hrs)
        if self.num_cov_tracks>0:
            gs1.update(hspace=gs1hs)

        for c in range(self.num_cov_tracks):
            final_plot = c==self.num_cov_tracks-1
            ax = plt.subplot(gs1[c,:])

            ax.fill_between(self.covx_lst[c], self.cov_lst[c],y2=0, color=color_dens, lw=0)

            maxheight = max(self.cov_lst[c])
            ymax = 1.1 * maxheight
            ymin = -.5 * ymax

            annotations = []

            if self.num_sj_tracks>0:
                for jxn,val in self.intron_cov_lst[c].items():
                    leftss, rightss = jxn

                    ss1, ss2 = [self.graphcoords[leftss - self.get_start() - 1],
                                self.graphcoords[rightss - self.get_start()]]

                    mid = (ss1 + ss2) / 2
                    h = -3 * ymin / 4

                    leftdens = self.cov_full_lst[c][int(ss1)]
                    rightdens = self.cov_full_lst[c][int(ss2)]

                    pts = [(ss1, leftdens),
                           (ss1, leftdens + h),
                           (ss2, rightdens + h),
                           (ss2, rightdens)]

                    midpt = Locus.cubic_bezier(pts, .5)

                    if self.settings["number_junctions"]:
                        annotations.append(ax.annotate('%s'%(val), xy=(midpt[0], midpt[1]), xytext=(midpt[0], midpt[1]+.3),fontsize=self.settings["font_size"]))

                    pp1 = PathPatch(Path(pts,[Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]),
                                    ec=color_spines, lw=np.log(val + 1) / np.log(10), fc='none')

                    ax.add_patch(pp1)

                adjust_text(annotations, autoalign='y', expand_objects=(0.1, 1),
                            only_move={'points':'', 'text':'y', 'objects':'y'}, force_text=0.75, force_objects=0.1)

            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
            max_graphcoords = max(self.graphcoords) - 1
            ax.set_xlim(0, max(self.graphcoords))

            if not final_plot:
                ax.spines['bottom'].set_color('none')
                ax.set_xticks([])
                ax.set_xticks([],minor=True)
            else:
                ax.xaxis.set_ticks_position('bottom')
                ax.set_xlabel("Genomic coordinates : "+self.get_coords_str(),fontsize=self.settings["font_size"])

                coords_fontsize = self.settings["font_size"] - (self.settings["font_size"] * 0.2)
                ax.set_xticks(np.linspace(0, max_graphcoords, self.settings["nxticks"]),
                              [self.graphToGene[int(x)] for x in \
                               np.linspace(0, max_graphcoords, self.settings["nxticks"])],
                              fontsize=coords_fontsize)


            ax.set_ylabel("Coverage",fontsize=self.settings["font_size"])
            ax.spines["left"].set_bounds(0, max(self.cov_full_lst[c]))
            ax.tick_params(axis='y',labelsize=self.settings["font_size"])
            ax.set_ybound(lower=ax.get_ybound()[0], upper=max(self.cov_full_lst[c]))
            ax.yaxis.set_ticks_position('left')

        exonwidth = .3
        narrows = 50

        locus_start = self.get_start()

        gs2 = GridSpec(len(self.txs)+self.num_cov_tracks, 1,height_ratios=hrs)
        gs2.update(hspace=gs2hs)
        for i,tx in enumerate(self.txs):
            ax2 = plt.subplot(gs2[i+self.num_cov_tracks,:])
            ax2.set_xlabel(tx.get_tid(),fontsize=self.settings["font_size"])

            for s, e in tx.orf:
                s = s - locus_start
                e = e - locus_start
                x = [self.graphcoords[s], self.graphcoords[e], self.graphcoords[e], self.graphcoords[s]]
                y = [-exonwidth / 6, -exonwidth / 6, exonwidth / 6, exonwidth / 6]
                ax2.fill(x, y,color=color_cds, lw=.5, zorder=30)

            for s, e in tx.exons:
                s = s - locus_start
                e = e - locus_start
                x = [self.graphcoords[s], self.graphcoords[e], self.graphcoords[e], self.graphcoords[s]]
                y = [-exonwidth / 8, -exonwidth / 8, exonwidth / 8, exonwidth / 8]
                ax2.fill(x, y, color=color_exon, lw=.5, zorder=20)

            # Draw intron.
            tx_start = tx.get_start() - locus_start
            tx_end = tx.get_end() - locus_start

            hline_left = self.graphcoords[tx_start]/max(self.graphcoords)
            hline_right = self.graphcoords[tx_end]/max(self.graphcoords)
            ax2.axhline(0,xmin=hline_left,xmax=hline_right, color=color_exon, lw=2)

            # Draw intron arrows.
            spread = .2 * max(self.graphcoords) / narrows
            for i in range(narrows):
                loc = float(i) * max(self.graphcoords) / narrows
                if tx.get_strand == '+' or self.settings["reverse_minus"]:
                    x = [loc - spread, loc, loc - spread]
                else:
                    x = [loc + spread, loc, loc + spread]
                y = [-exonwidth / 20, 0, exonwidth / 20]
                if x[0]>=self.graphcoords[tx_start] and x[0]<=self.graphcoords[tx_end]:
                    ax2.plot(x, y, lw=2, color=color_exon)

            ax2.set_xlim(0, max(self.graphcoords))
            plt.box(on=False)
            ax2.set_xticks([])
            ax2.set_yticks([])

        plt.subplots_adjust(hspace=.5, wspace=.7)
        plt.savefig(out_fname)


def sashimi(args):
    assert os.path.exists(args.gtf), "GTF does not exist: " + args.gtf
    # assert os.path.exists(args.cov), "Coverage file does not exist: " + args.cov
    # assert os.path.exists(args.sj), "Splice Junction file does not exist: " + args.sj

    settings = {"intron_scale": args.intron_scale,
                "exon_scale": args.exon_scale,
                "logged": args.logged,
                "ymax": args.ymax,
                "number_junctions": args.number_junctions,
                "resolution": args.resolution,
                "fig_width": args.fig_width,
                "fig_height": args.fig_height,
                "junction_log_base": args.junction_log_base,
                "reverse_minus": args.reverse_minus,
                "font_size": args.font_size,
                "nyticks": args.nyticks,
                "nxticks": args.nxticks,
                "show_ylabel": args.show_ylabel,
                "show_xlabel": args.show_xlabel,
                "sans_serif": args.sans_serif,
                "bar_color": args.bar_color}

    tids_seen = set()

    locus = Locus()
    locus.set_settings(settings)

    with open(args.gtf, "r") as inFP:
        cur_tid = None
        cur_tid_lines = ""
        for line in inFP:
            lcs = line.split("\t")
            if not len(lcs) == 9:
                continue

            tid = lcs[8].split("transcript_id \"", 1)[1].split("\"", 1)[0]
            if cur_tid is None:
                cur_tid = tid
                cur_tid_lines = ""

            if not cur_tid == tid:
                assert tid not in tids_seen, "mixed tids"
                tx = TX()
                tx.parse_from_gtf(cur_tid_lines)
                locus.add_tx(tx)
                cur_tid = tid
                cur_tid_lines = line

            else:
                cur_tid_lines += line

        tx = TX()
        tx.parse_from_gtf(cur_tid_lines)
        locus.add_tx(tx)

    locus.set_scaling()

    # read in only values for which the transcriptome has been constructed
    is_cov_lst_file = args.cov is not None
    if args.cov is not None:
        # check if it's a file listing a set of files with introns
        with open(args.cov,"r") as inFP:
            for line in inFP:
                tmp = line.strip()
                if not os.path.exists(tmp) and len(tmp)>0:
                    is_cov_lst_file = False
                    break

        if is_cov_lst_file:
            with open(args.cov,"r") as inFP:
                for line in inFP:
                    tmp = line.strip()
                    if len(tmp)>0:
                        locus.add_coverage(tmp)
        else:
            locus.add_coverage(args.cov)

    # add coverage
    is_sj_lst_file = True
    if args.sj is not None:
        with open(args.sj,"r") as inFP:
            for line in inFP:
                tmp = line.strip()
                if not os.path.exists(tmp) and len(tmp)>0:
                    is_sj_lst_file = False
                    break

        if is_sj_lst_file:
            assert is_cov_lst_file,"can not add splice junction tracks as list without coverage tracks provided as list as well"
            with open(args.sj,"r") as inFP:
                for line in inFP:
                    tmp = line.strip()
                    if len(tmp)>0:
                        locus.add_introns(tmp)
        else:
            locus.add_introns(args.sj)

    locus.plot(args.output)


def main(args):
    parser = argparse.ArgumentParser(description='''Help Page''')
    parser.add_argument("--gtf",
                        required=True,
                        help="annotation in a GFF/GTF format")
    parser.add_argument("--cov",
                        required=False,
                        help="coverage in bedgraph format or a file containing a list of filenames with coverage in bedgraph for multiple samples. If a list is provided - the files should be in the same order as the splice junctions below (if provided)")
    parser.add_argument("--sj",
                        required=False,
                        help="splice junctions in bed format or a file containing a list of filenames with splice junctions in bed format for multiple samples. If a list is provided - the files should be in the same order as the coverage tracks.")
    parser.add_argument("-o",
                        "--output",
                        required=True,
                        help="output basename")
    parser.add_argument("--intron_scale",
                        required=False,
                        type=int,
                        default=20,
                        help="intron_scale")
    parser.add_argument("--exon_scale",
                        required=False,
                        type=int,
                        default=1,
                        help="exon_scale")
    parser.add_argument("--resolution",
                        required=False,
                        type=int,
                        default=6,
                        help="resolution")
    parser.add_argument("--fig_width",
                        required=False,
                        type=int,
                        default=20,
                        help="fig_width")
    parser.add_argument("--fig_height",
                        required=False,
                        type=int,
                        default=10,
                        help="fig_height")
    parser.add_argument("--junction_log_base",
                        required=False,
                        type=float,
                        default=10.,
                        help="junction_log_base")
    parser.add_argument("--font_size",
                        required=False,
                        type=int,
                        default=18,
                        help="fig_height")
    parser.add_argument("--nyticks",
                        required=False,
                        type=int,
                        default=3,
                        help="nyticks")
    parser.add_argument("--nxticks",
                        required=False,
                        type=int,
                        default=4,
                        help="nxticks")
    parser.add_argument("--ymax",
                        required=False,
                        type=int,
                        default=None,
                        help="ymax")
    parser.add_argument("--logged",
                        action="store_true",
                        required=False,
                        help="logged - False by default")
    parser.add_argument("--number_junctions",
                        action="store_false",
                        required=False,
                        help="number_junctions - True by default")
    parser.add_argument("--reverse_minus",
                        required=False,
                        action="store_true",
                        help="reverse_minus - False by default")
    parser.add_argument("--show_ylabel",
                        required=False,
                        action="store_false",
                        help="show_ylabel - True by default")
    parser.add_argument("--show_xlabel",
                        required=False,
                        action="store_false",
                        help="show_xlabel - True by default")
    parser.add_argument("--sans_serif",
                        required=False,
                        action="store_true",
                        help="sans_serif - False by default")
    parser.add_argument("--bar_color",
                        required=False,
                        type=str,
                        default="k",
                        help="bar_color")

    parser.set_defaults(func=sashimi)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
