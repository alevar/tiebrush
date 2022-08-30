#!/usr/bin/env python

import numpy as np
import argparse
import pickle
import sys
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import Patch
from matplotlib.patches import PathPatch
from matplotlib.gridspec import GridSpec
from adjustText import adjust_text

mpl.use('Agg')

def intersect(s1,s2):
    res = [0,-1,0]
    tis = max(s1[0],s2[0])
    tie = min(s1[1],s2[1])
    if(tis<=tie):
        res[0] = tis
        res[1] = tie
        return (tie-tis)+1,res
    return 0,res

def split(s1,s2):
    left  = [0,-1,-1]
    right = [0,-1,-1]

    il,inter = intersect(s1,s2)
    if il>0:
        if inter[0]>s1[0]:
            left[0] = s1[0]
            left[1] = inter[0]-1
            left[2] = s1[2]
        if inter[1]<s1[1]:
            right[0] = inter[1]+1
            right[1] = s1[1]
            right[2] = s1[2]
    else:
        if s1[0]<s2[0]:
            left = s1
        else:
            right = s1

    return left,inter,right

def slen(s):
    return (s[1]-s[0])+1

def clen(chain):
    res = 0
    for c in chain:
        res+=slen(c)
    return res

def compare(i1,i2):
    intervals = []
    for i in i1:
        intervals.append([i[0],i[1]])
        intervals[-1].append(-1)
    for i in i2:
        intervals.append([i[0],i[1]])
        intervals[-1].append(1)
    intervals.sort()

    if len(i1)==0 and len(i2)==0:
        return []

    stack = []
    stack.append(intervals[0])
    for i in intervals[1:]:

        left,inter,right = split(stack[-1],i)
        if slen(right)>0:
            assert slen(inter)==slen(i) # must be intirely contained within
        else:
            tmp,inter2,right = split(i,stack[-1])
            if(slen(tmp)>0):
                t2 = stack[-1]
                stack[-1]=tmp
                stack.append(t2)

            else:
                assert slen(tmp)<=0,str(tmp)+","+str(inter2)+","+str(right)
            assert inter==inter2

        stack.pop()

        if slen(left)>0:
            stack.append(left)
        if slen(inter)>0:
            inter[2] = 0
            stack.append(inter)
        if slen(right)>0:
            stack.append(right)

    return stack

# runs compare() funciton and labels all matches as in and out of frame accordingly
def compare_label_frame(chain1,chain2,strand):
    if chain2 is np.nan or len(chain2)==0:
        [[x[0],x[1],-1] for x in chain1]
    if chain1 is np.nan or len(chain1)==0:
        [[x[0],x[1],1] for x in chain2]

    mod_chain = compare(chain1,chain2)

    if strand=="-":
        mod_chain.reverse()

    t_frame = 0
    q_frame = 0

    for mc in mod_chain:
        if (mc[2] == -1): # extra positions in the query
            q_frame += slen(mc)
        elif (mc[2] == 1): # template positions missing from the query
            t_frame += slen(mc)
        elif (mc[2] == 0): # matching positions between query and template
            if (q_frame % 3 == t_frame % 3):
                mc[2] = 100 # inframe
            else:
                mc[2] = -100 # outframe
        else:
            print("wrong code: "+str(mc[2]))
            return

    return mod_chain

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
        self.ref_tx = None
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

    def add_tx(self, tx, ref=False):
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

        if ref:
            self.ref_tx = len(self.txs)-1

    def set_scaling(self):
        # get graphcoords
        if self.graphcoords is None:
            self.graphcoords, self.graphToGene = self.getScaling(self.settings["intron_scale"], self.settings["exon_scale"],
                                                                 self.settings["reverse"])
            
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
                if intron in self.introns or self.settings["all_junctions"]:
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

    def plot(self,out_fname,title,save_pickle,compare):
        assert self.num_sj_tracks==self.num_cov_tracks or self.num_sj_tracks==0,"incompatible number of splice junciton and coverage tracks - the numebrs should either be the same or no splice junction tracks provided at all"

        # sns.color_palette("colorblind").as_hex()[3]
        color_dens = "#ffb703"
        color_spines = "#fb8500"
        colors_compare = {-1:("#029e73","Missing From Reference"),
                          1:("#949494","Extra In Reference"),
                          100:("#56b4e9","Matching In Frame"),
                          -100:("#d55e00","Matching Out Of Frame"),
                          0:("#023047","Non-Coding Positions")}
        colors_non_compare = {100:("#56b4e9","Coding Positions"),
                              0:("#023047","Non-Coding Positions")}

        hrs = [4]*self.num_cov_tracks+[1 for i in range(len(self.txs))]
        gs1hs = 1
        gs2hs = 0.4

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

        # gs2 = GridSpec(len(self.txs)+self.num_cov_tracks, 1,height_ratios=hrs)
        gs1.update(hspace=gs2hs)
        for i,tx in enumerate(self.txs):
            ax2 = plt.subplot(gs1[i+self.num_cov_tracks,:])
            ax2.set_xlabel(tx.get_tid(),fontsize=self.settings["font_size"])

            if compare and self.ref_tx is not None:
                stack = compare_label_frame(tx.orf,self.txs[self.ref_tx].orf,self.strand)

                for s, e, l in stack:
                    s = s - locus_start
                    e = e - locus_start
                    x = [self.graphcoords[s], self.graphcoords[e], self.graphcoords[e], self.graphcoords[s]]
                    y = [-exonwidth / 6, -exonwidth / 6, exonwidth / 6, exonwidth / 6]
                    if l==1:
                        ax2.fill(x, y,linestyle="-",color=colors_compare[l][0],lw=2,zorder=30,fill=False)
                    else:
                        ax2.fill(x, y,color=colors_compare[l][0], lw=.5, zorder=30)

            else:
                for s, e in tx.orf:
                    s = s - locus_start
                    e = e - locus_start
                    x = [self.graphcoords[s], self.graphcoords[e], self.graphcoords[e], self.graphcoords[s]]
                    y = [-exonwidth / 6, -exonwidth / 6, exonwidth / 6, exonwidth / 6]
                    ax2.fill(x, y,color=colors_non_compare[100], lw=.5, zorder=30)

            for s, e in tx.exons:
                s = s - locus_start
                e = e - locus_start
                x = [self.graphcoords[s], self.graphcoords[e], self.graphcoords[e], self.graphcoords[s]]
                y = [-exonwidth / 8, -exonwidth / 8, exonwidth / 8, exonwidth / 8]
                ax2.fill(x, y, color=colors_non_compare[0][0], lw=.5, zorder=20)

            # Draw intron.
            tx_start = tx.get_start() - locus_start
            tx_end = tx.get_end() - locus_start

            hline_left = self.graphcoords[tx_start]/max(self.graphcoords)
            hline_right = self.graphcoords[tx_end]/max(self.graphcoords)
            ax2.axhline(0,xmin=hline_left,xmax=hline_right, color=colors_non_compare[0][0], lw=2)

            # Draw intron arrows.
            spread = .2 * max(self.graphcoords) / narrows
            for i in range(narrows):
                loc = float(i) * max(self.graphcoords) / narrows
                if tx.get_strand() == '+' or self.settings["reverse"]:
                    x = [loc - spread, loc, loc - spread]
                else:
                    x = [loc + spread, loc, loc + spread]
                y = [-exonwidth / 20, 0, exonwidth / 20]
                if x[0]>=self.graphcoords[tx_start] and x[0]<=self.graphcoords[tx_end]:
                    ax2.plot(x, y, lw=2, color=colors_non_compare[0][0])

            ax2.set_xlim(0, max(self.graphcoords))
            plt.box(on=False)
            ax2.set_xticks([])
            ax2.set_yticks([])


        if not title is None:
            plt.suptitle(title,wrap=True)

        plt.subplots_adjust(hspace=.5, wspace=.7)

        handles = [Patch(color=c[0],label=c[1]) for i,c in colors_non_compare.items()]
        if compare:
            handles = [Patch(color=c[0],label=c[1]) for i,c in colors_compare.items()]
            handles[1] = Patch(edgecolor=colors_compare[1][0],facecolor=None,fill=False,linestyle="-",linewidth=3,label=colors_compare[1][1])

        lgd = plt.legend(handles=handles,
                       fontsize=self.settings["font_size"],
                       loc="lower left",
                       bbox_to_anchor=(0., -2.02, 1., .102),
                       ncol=2,
                       mode="expand",
                       borderaxespad=0.)

        if save_pickle:
            with open(out_fname+".pickle","wb") as outFP:
                pickle.dump((fig,gs1),outFP)

        plt.savefig(out_fname,bbox_extra_artists=(lgd,), bbox_inches='tight')


def sashimi(args):
    assert os.path.exists(args.gtf), "GTF does not exist: " + args.gtf
    # assert os.path.exists(args.cov), "Coverage file does not exist: " + args.cov
    # assert os.path.exists(args.sj), "Splice Junction file does not exist: " + args.sj

    settings = {"intron_scale": args.intron_scale,
                "exon_scale": args.exon_scale,
                "number_junctions": args.number_junctions,
                "resolution": args.resolution,
                "fig_width": args.fig_width,
                "fig_height": args.fig_height,
                "reverse": args.reverse,
                "font_size": args.font_size,
                "nxticks": args.nxticks,
                "title": args.title,
                "pickle": args.pickle,
                "compare": args.compare,
                "all_junctions": args.all_junctions}

    tids_seen = set()

    locus = Locus()
    locus.set_settings(settings)

    found_ref = False

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

                is_ref = args.compare is not None and tx.get_tid()==args.compare
                found_ref = found_ref or is_ref
                locus.add_tx(tx,is_ref)

                cur_tid = tid
                cur_tid_lines = line

            else:
                cur_tid_lines += line

        tx = TX()
        tx.parse_from_gtf(cur_tid_lines)
        is_ref = args.compare is not None and tx.get_tid()==args.compare
        found_ref = found_ref or is_ref
        locus.add_tx(tx,is_ref)

    locus.set_scaling()

    if compare and not found_ref:
        print("could not find the reference transcript for comparison: "+args.compare)
        exit(1)

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

    title = None
    if not args.title is None:
        title = " ".join(args.title)
    locus.plot(args.output,title,args.pickle,args.compare)


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
                        help="Filename for the output figure. The format (png,svg, ...) will be automatically deduced based on the extension.")
    parser.add_argument("--intron_scale",
                        required=False,
                        type=int,
                        default=20,
                        help="Parameter regulating the scaling of the introns (Default: 20). Decreasing the integer value will scale introns down in size compared to exons.")
    parser.add_argument("--exon_scale",
                        required=False,
                        type=int,
                        default=1,
                        help="Parameter regulating the scaling of the exons (Default: 1). Increasing the integer value will scale exons down in size compared to introns.")
    parser.add_argument("--resolution",
                        required=False,
                        type=int,
                        default=6,
                        help="Parameter regulates the smoothing factor of the coverage track (Default: 6). Increasing the value will increasing the smoothing by reducing the number of points on the coverage track.")
    parser.add_argument("--fig_width",
                        required=False,
                        type=int,
                        default=20,
                        help="Width of the figure in inches (Default: 20).")
    parser.add_argument("--fig_height",
                        required=False,
                        type=int,
                        default=10,
                        help="Height of the figure in inches (Default: 10).")
    parser.add_argument("--font_size",
                        required=False,
                        type=int,
                        default=18,
                        help="Size of the font (Default: 18)")
    parser.add_argument("--nxticks",
                        required=False,
                        type=int,
                        default=4,
                        help="Number of positional markers to include on the x-axis with labels (Default: 4).")
    parser.add_argument("--number_junctions",
                        action="store_false",
                        required=False,
                        help="Disables labels idicating coverage of splice junctions")
    parser.add_argument("--reverse",
                        required=False,
                        action="store_true",
                        help="Flips image horizontally, which is equivalent to setting strand to the opposite value.")
    parser.add_argument("--title",
                        required=False,
                        nargs='+',
                        default=None,
                        help="Title of the figure.")
    parser.add_argument("--pickle",
                        required=False,
                        action="store_true",
                        help="Save a pickle alongside the figure which can be loaded into a separate instance of matplotlib for modification.")
    parser.add_argument("--compare",
                        required=False,
                        type=str,
                        help="Users can specify one of the input transcripts to serve as a reference. If set, all transcripts in the input will be compared to the reference and plotted using a dedicated color pallete. The comparison will visualize in-frame and out-of-frame positions as well as any intervals missing and extra between the reference and each query transcript")
    parser.add_argument("--all-junctions",
                        required=False,
                        action="store_true",
                        help="Will force the script to display all junctions, including those not present in the GTF")

    parser.set_defaults(func=sashimi)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
