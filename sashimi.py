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

def cut(chain,s,e):
    dup = [c for c in chain]
    if s<=dup[0][0] and e>=dup[-1][1]: # nothing to cut
        return dup
        
    tmp = list()
    for i,c in enumerate(dup): # not reference - the copy is being modified
        if c[0]<=s and c[1]>=s:
            dup[i][0]=s
            c[0]=s
        if c[1]>=e:
            dup[i][1]=e
            c[1]=e
        if c[1]<s or c[0]>e:
            continue
            
        tmp.append(c)
        
    return tmp

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
        self.dummy = False
        self.store_cds = True

        self.attrs = dict()

    def clear(self):
        self.seqid = None
        self.strand = None
        self.exons = []
        self.orf = []

        self.tid = None
        self.dummy = False
        self.store_cds = True

        self.attrs = dict()

    def set_nocds(self):
        self.store_cds = False
    def set_seqid(self,seqid):
        self.seqid = seqid
    def set_strand(self,strand):
        self.strand = strand
    def set_exons(self,exons):
        self.exons = exons
    def set_orf(self,orf):
        self.orf = orf
    def set_tid(self,tid):
        self.tid = tid
    def set_dummy(self,dummy):
        self.dummy=dummy

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
                # extract attributes
                for av in lcs[8].strip().rstrip(";").split(";"):
                    # print(lcs[8].strip().rstrip(";"),av)
                    k,v = av.strip().rstrip("\"").split(" \"")
                    assert k not in self.attrs,"duplicate attribute key: "+k+" in: "+lcs[8].strip()
                    self.attrs[k]=v
                continue

            if lcs[2].lower() == "exon":
                self.exons.append([int(lcs[3]), int(lcs[4])])
            if self.store_cds and lcs[2].lower() == "cds":
                self.orf.append([int(lcs[3]), int(lcs[4])])

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

    def get_attr(self,attr):
        return self.attrs.get(attr,"")

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
        self.track_names = list()

        self.exon_starts = []
        self.exon_ends = []

        self.graphcoords = None
        self.graphToGene = None
        
        self.covx_lst = list()
        self.cov_lst = list()
        self.cov_full_lst = list()
        
        # these contain values for the zoomed-in view if requested
        self.covx_zoom_lst = list()
        self.cov_zoom_lst = list()
        self.cov_full_zoom_lst = list()

        self.settings = None

        self.num_cov_tracks = 0
        self.num_sj_tracks = 0
        
        self.zoom_ratio = 1
        
        self.color_dens = "#ffb703"
        self.color_spines = "#fb8500"
        self.colors_compare = {-1:("#029e73","Missing From Reference"),
                          1:("#949494","Extra In Reference"),
                          100:("#56b4e9","Matching In Frame"),
                          -100:("#d55e00","Matching Out Of Frame"),
                          0:("#023047","Non-Coding Positions")}
        self.colors_exon_compare = {-1:("#029e73","Missing From Reference"),
                          1:("#949494","Extra In Reference"),
                          100:("#023047","Matching")}
        self.colors_non_compare = {100:("#56b4e9","Coding Positions"),
                              0:("#023047","Non-Coding Positions")}

    @staticmethod
    def union(intervals):
        res = []
        for s, e in sorted(intervals):
            if res and res[-1][1] >= s - 1:
                res[-1][1] = max(res[-1][1], e)
            else:
                res.append([s, e])

        return [[x[0], x[1]] for x in res]

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
        assert self.seqid is None or self.seqid == tx.get_seqid(), "mismatching seqids: "+tx.get_tid()
        assert self.strand is None or self.strand == tx.get_strand(), "mismatching strands: "+tx.get_tid()

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
                                                                 
        # get ratios
        if self.settings["zoom"]:
            zoom_start_transform = self.settings["zoom_start"]-self.get_start()
            zoom_end_transform = self.settings["zoom_end"]-self.get_start()
            self.zoom_ratio = 1 if not self.settings["zoom"] else (zoom_end_transform-zoom_start_transform)/(self.get_end()-self.get_start())

    def get_start(self):
        return self.intervals[0][0]

    def get_end(self):
        return self.intervals[-1][1]

    def add_track_names(self,name_fname):
        if os.path.exists(name_fname):
            with open(name_fname, "r") as inFP:
                for line in inFP:
                    self.track_names.append(line.strip())
        else:
            self.track_names.append(name_fname)

    def add_introns(self,sj_fname):
        assert os.path.exists(sj_fname),"Splice Junction track does not exist"

        self.intron_cov_lst.append(dict())

        total_sj_cov = 0
        max_sj_cov = 0
        min_sj_cov = sys.maxsize

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
                    self.intron_cov_lst[-1][intron] = [int(lcs[4])]
                    total_sj_cov += int(lcs[4])
                    max_sj_cov = max(max_sj_cov,int(lcs[4]))
                    min_sj_cov = min(min_sj_cov,int(lcs[4]))

        factor = 0.00001
        try:
            factor = total_sj_cov/len(self.intron_cov_lst[-1])
        except:
            pass
        for k,v in self.intron_cov_lst[-1].items():
            self.intron_cov_lst[-1][k].append(round(self.intron_cov_lst[-1][k][0]/factor,2))
            top = 99*self.intron_cov_lst[-1][k][0]
            bottom = max(0.1,max_sj_cov)
            self.intron_cov_lst[-1][k].append((top/bottom)+1)

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

    def compress_intervals(self, vals, graphcoords,resolution):  # intervals with optional values if more
        compressed_x = []
        compressed_wiggle = []
        prevx = graphcoords[0]
        tmpval = []
        for i in range(len(graphcoords)):
            tmpval.append(vals[i])
            if abs(graphcoords[i] - prevx) > resolution:
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
        self.covx_lst[-1], self.cov_lst[-1] = self.compress_intervals(self.cov_full_lst[-1], self.graphcoords,self.settings["resolution"])
        self.num_cov_tracks+=1
               
        if self.settings["zoom"]:            
            self.covx_zoom_lst.append(list())
            self.cov_zoom_lst.append(list())
            
            zoom_resolution = max(1,int(self.settings["resolution"]*self.zoom_ratio))
            
            self.covx_zoom_lst[-1],self.cov_zoom_lst[-1] = self.compress_intervals(self.cov_full_lst[-1], self.graphcoords,zoom_resolution)
            

    def get_coords_str(self):
        return self.seqid + self.strand + ":" + str(self.get_start()) + "-" + str(self.get_end())

    @staticmethod
    def scale_val(v, min_new, max_new, min_cur, max_cur):
        return (max_new - min_new) * (v - min_cur) / (max_cur - min_cur) + min_new;
        
    # build general gridspace
    # build general grid and populate axis
    # populate data for each
    # separate funciton to draw zoom boxes
    
    def build_gridspace(self,fig):
        gs_main = None
        if self.settings["zoom"]:
            gs_main = GridSpec(2, 1, figure=fig, height_ratios=[2,1], hspace=0.2)
        else:
            gs_main = GridSpec(1, 1, figure=fig)
        
        return gs_main
        
    def build_cov_tx_grid(self,fig,gs,zoom_view):
        height_ratio_cov2tx = 0
        try:
            if zoom_view:
                height_ratio_cov2tx = (self.settings["cov_height"])/self.settings["tx_height"]
            else:
                height_ratio_cov2tx = self.settings["cov_height"]/self.settings["tx_height"]
        except:
            pass
        
        # first separate into cov and tx sections
        gs_sub_cov = None
        gs_sub_tx = None
        if self.num_cov_tracks>0:
            cov_hr = height_ratio_cov2tx*self.num_cov_tracks
            tx_hr = len(self.txs)
            gs_sub = gs.subgridspec(2, 1, height_ratios=[cov_hr,tx_hr],hspace=0.2)
            gs_sub_cov = gs_sub[0].subgridspec(self.num_cov_tracks,1,hspace=1)
            gs_sub_tx = gs_sub[1].subgridspec(len(self.txs),1,hspace=0.4)
        else:
            gs_sub_tx = gs.subgridspec(len(self.txs), 1,hspace=0.4)
            
        return gs_sub_cov,gs_sub_tx
    
    @staticmethod
    def get_belly_arc_coords(start,end,leftdens,rightdens,h,thickness):
        pts = [(start, leftdens),
               (start, (leftdens+h)+thickness),
               (end,   (rightdens+h)+thickness),
               (end,   rightdens),
               (end,   (rightdens+h)-thickness),
               (start, (leftdens+h)-thickness),
               (start, leftdens),
               (start, leftdens)]
        
        codes = [Path.MOVETO,
                 Path.CURVE4,
                 Path.CURVE4,
                 Path.CURVE4,
                 Path.CURVE4,
                 Path.CURVE4,
                 Path.CURVE4,
                 Path.CLOSEPOLY]
        
        midpt = Locus.cubic_bezier([(start, leftdens),
                                    (start, leftdens+h),
                                    (end,   rightdens+h),
                                    (end,   rightdens)], .5)
        return pts,codes,midpt
    
    def plot_coverage(self,fig,gs,title,compare,text_attr,rel,zoom_view):    
        axes = []

        for c in range(self.num_cov_tracks):
            covx = self.covx_lst[c]
            cov = self.cov_lst[c]
            if zoom_view:
                covx = self.covx_zoom_lst[c]
                cov = self.cov_zoom_lst[c]
                

            final_plot = c==self.num_cov_tracks-1
            ax = fig.add_subplot(gs[c,:])
            axes.append(ax)

            if c==0 and not title is None:
                ax.set_title(title,wrap=True,fontsize=self.settings["font_size"]*1.5)

            if not zoom_view:
                ax.fill_between(covx,cov,y2=0, color=self.color_dens, lw=0)
            else:
                ax.fill_between(covx,cov,y2=0, color=self.color_dens, lw=0,step="post")

            annotations = []
            
            y_limits = [0,max(self.cov_full_lst[c])]
            
            min_height = 0.25*max(self.cov_full_lst[c])
            max_height = 0.4*max(self.cov_full_lst[c])

            if self.num_sj_tracks>0:
                for jxn,val in self.intron_cov_lst[c].items():
                    leftss, rightss = jxn

                    ss1, ss2 = [self.graphcoords[leftss - self.get_start() - 1],
                                self.graphcoords[rightss - self.get_start()]]

                    leftdens = self.cov_full_lst[c][leftss - self.get_start()-1]
                    rightdens = self.cov_full_lst[c][rightss - self.get_start()]
                    maxdens = max(leftdens,rightdens)
                    
                    h = min(maxdens*0.75,max_height)
                    h = max(h,min_height)
                    
                    # for the super small values - we need to make the height higher
                    # for super tall values - need smaller
                    # else standardized
                    
                    thickness = min(max(self.cov_full_lst[c])*0.1,val[0]*0.1)
                    y_limits[1] = max(y_limits[1],maxdens+h+(thickness/2))
                    
                    pts,codes,midpt = self.get_belly_arc_coords(ss1,ss2,leftdens,rightdens,h,thickness)

                    if self.settings["number_junctions"]:
                        sj_v = val[1] if rel else val[0]
                        annotations.append(ax.annotate('%s'%(sj_v), xy=(midpt[0], midpt[1]), xytext=(midpt[0], midpt[1]+.3),fontsize=self.settings["font_size"]))

                    pp1 = PathPatch(Path(pts,codes),ec=self.color_spines, fc=self.color_spines,alpha=0.75,lw=2/3)

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


            ax.set_ylim(y_limits[0],y_limits[1])
            ax.set_ybound(lower=y_limits[0], upper=y_limits[1])
            ax.yaxis.set_ticks_position('left')
            if len(self.track_names)>0:
                ax.set_xlabel(self.track_names[c],fontsize=self.settings["font_size"])
                
        return axes
    
    def plot_txs(self,fig,gs,title,compare,text_attr,rel,zoom_view):
        
        axes = []
        
        exonwidth = .3
        narrows = 50
        if zoom_view:
            narrows /= self.zoom_ratio
            narrows = int(narrows)

        locus_start = self.get_start()

        for i,tx in enumerate(self.txs):
            if tx.dummy: # skip any dummies
                continue

            ax = fig.add_subplot(gs[i,:])
            axes.append(ax)

            if i==0 and not title is None and self.num_cov_tracks==0:
                ax.set_title(title,wrap=True,fontsize=self.settings["font_size"]*1.5)

            xlabel = tx.get_tid()
            if text_attr != "transcript_id":
                ta = tx.get_attr(text_attr)
                if len(ta)>0:
                    xlabel = tx.get_attr(text_attr)

            ax.set_xlabel(xlabel,fontsize=self.settings["font_size"])

            if compare and self.ref_tx is not None:
                stack = compare_label_frame(tx.orf,self.txs[self.ref_tx].orf,self.strand)
                for s, e, l in stack:                        
                    s = s - locus_start
                    e = e - locus_start
                    x = [self.graphcoords[s], self.graphcoords[e], self.graphcoords[e], self.graphcoords[s]]
                    y = [-exonwidth / 6, -exonwidth / 6, exonwidth / 6, exonwidth / 6]
                    
                    ax.fill(x, y,linestyle="-",color=self.colors_compare[l][0],lw=2,zorder=30,fill=False)
                    if not l==1:
                        ax.fill(x, y,color=self.colors_compare[l][0], zorder=30)

                if self.ref_tx == i:
                    ax.set_facecolor((1,0,0,0.1))
                    #x = [self.graphcoords[0], self.graphcoords[self.get_end()-self.get_start()], self.graphcoords[self.get_end()-self.get_start()], self.graphcoords[0]]
                    #y = [-exonwidth / 5, -exonwidth / 5, exonwidth / 5, exonwidth / 5]
                    #ax.fill(x, y,linestyle="-",color="xkcd:salmon",alpha=0.15,lw=2,zorder=30,fill=True)

            else:
                cur_orf = [o for o in tx.orf]
                for s, e in cur_orf:
                    s = s - locus_start
                    e = e - locus_start
                    x = [self.graphcoords[s], self.graphcoords[e], self.graphcoords[e], self.graphcoords[s]]
                    y = [-exonwidth / 6, -exonwidth / 6, exonwidth / 6, exonwidth / 6]
                    ax.fill(x, y,linestyle="-",color=self.colors_non_compare[100][0],lw=2,zorder=30,fill=False)
                    ax.fill(x, y,color=self.colors_non_compare[100][0], lw=.5, zorder=30)

            cur_exons = [e for e in tx.exons]  
            for s, e in cur_exons:
                s = s - locus_start
                e = e - locus_start
                x = [self.graphcoords[s], self.graphcoords[e], self.graphcoords[e], self.graphcoords[s]]
                y = [-exonwidth / 8, -exonwidth / 8, exonwidth / 8, exonwidth / 8]
                ax.fill(x, y, color=self.colors_non_compare[0][0], lw=.5, zorder=20)

            # Draw intron.
            tx_start = tx.get_start() - locus_start
            tx_end = tx.get_end() - locus_start
            max_val = max(self.graphcoords)

            hline_left = self.graphcoords[tx_start]/max_val
            hline_right = self.graphcoords[tx_end]/max_val
            ax.axhline(0,xmin=hline_left,xmax=hline_right, color=self.colors_non_compare[0][0], lw=2)

            # Draw intron arrows.
            spread = .2 * max_val / narrows
                
            for im in range(narrows):
                loc = float(im) * max_val / narrows
                if tx.get_strand() == '+' or self.settings["reverse"]:
                    x = [loc - spread, loc, loc - spread]
                else:
                    x = [loc + spread, loc, loc + spread]
                y = [-exonwidth / 20, 0, exonwidth / 20]
                if x[0]>=self.graphcoords[tx_start] and x[0]<=self.graphcoords[tx_end]:
                    ax.plot(x, y, lw=2, color=self.colors_non_compare[0][0])

            ax.set_xlim(0, max_val)
            #plt.box(on=False)
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
            ax.spines['bottom'].set_color('none')
            ax.spines['left'].set_color('none')
            ax.set_xticks([])
            ax.set_yticks([])

                
        return axes
    
    def build_zoom(self,fig,axes_norm,axes_zoom):
        # first need coordinates for the zoom region
        zoom_start_transform = self.graphcoords[self.settings["zoom_start"]-self.get_start()]
        zoom_end_transform = self.graphcoords[self.settings["zoom_end"]-self.get_start()]
        
        num_tracks = len(self.txs)+self.num_cov_tracks
        
        scaled_start_zoom = self.scale_val(zoom_start_transform,
                                      axes_norm[0].get_position().get_points()[0][0],
                                      axes_norm[0].get_position().get_points()[1][0],
                                      axes_norm[0].get_xlim()[0],
                                      axes_norm[0].get_xlim()[1])
        scaled_end_zoom = self.scale_val(zoom_end_transform,
                                      axes_norm[0].get_position().get_points()[0][0],
                                      axes_norm[0].get_position().get_points()[1][0],
                                      axes_norm[0].get_xlim()[0],
                                      axes_norm[0].get_xlim()[1])
        
        scaled_start_01 = self.scale_val(zoom_start_transform,
                                      0,
                                      1,
                                      axes_norm[0].get_xlim()[0],
                                      axes_norm[0].get_xlim()[1])
        scaled_end_01 = self.scale_val(zoom_end_transform,
                                      0,
                                      1,
                                      axes_norm[0].get_xlim()[0],
                                      axes_norm[0].get_xlim()[1])
        
        
        top_gs = GridSpec(2, 1, figure=fig)
        top_gs.update(top=axes_norm[0].get_position().get_points()[1][1],
                        bottom=axes_norm[-1].get_position().get_points()[0][1],
                        left=scaled_start_zoom,
                        right=scaled_end_zoom)
        top_ax = fig.add_subplot(top_gs[0:])
        top_ax.tick_params(axis='both',which='both',bottom=0,left=0, labelbottom=0, labelleft=0)
        top_ax.set_facecolor("grey")
        top_ax.patch.set_alpha(0.1)
        
        mid_gs = GridSpec(2, 1, figure=fig)
        mid_gs.update(top=axes_norm[-1].get_position().get_points()[0][1],
                        bottom=axes_zoom[0].get_position().get_points()[1][1],
                        left=axes_norm[0].get_position().get_points()[0][0],
                        right=axes_norm[0].get_position().get_points()[1][0])
        mid_ax = fig.add_subplot(mid_gs[0:])
        mid_ax.tick_params(axis='both',which='both',bottom=0,left=0,labelbottom=0,labelleft=0)
        p = plt.Polygon([[scaled_start_01,1],[0,0],[1,0],[scaled_end_01,1]],closed=True,color="grey",alpha=0.1)
        mid_ax.add_patch(p)
        mid_ax.set_axis_off()

        bot_gs = GridSpec(2, 1, figure=fig)
        bot_gs.update(top=axes_zoom[0].get_position().get_points()[1][1],
                        bottom=axes_zoom[-1].get_position().get_points()[0][1],
                        left=axes_norm[0].get_position().get_points()[0][0],
                        right=axes_norm[0].get_position().get_points()[1][0])
        bot_ax = fig.add_subplot(bot_gs[0:])
        bot_ax.tick_params(axis='both',which='both',bottom=0,left=0,labelbottom=0,labelleft=0)
        bot_ax.set_facecolor("grey")
        bot_ax.patch.set_alpha(0.1)
        
        
        
    def build_fig(self,out_fname,title,save_pickle,compare,legend,text_attr,rel):
        assert self.num_sj_tracks==self.num_cov_tracks or self.num_sj_tracks==0,"incompatible number of splice junciton and coverage tracks - the numebrs should either be the same or no splice junction tracks provided at all"
        
        fig_height = (self.settings["tx_height"]*len(self.txs)) + (self.settings["cov_height"]*self.num_cov_tracks)
        if self.settings["zoom"]:
            fig_height *= 1.75
            
        fig = plt.figure(figsize=(self.settings["fig_width"],fig_height))
        
        gs_main = self.build_gridspace(fig)
        gs_subs = []
        
        for nr in range(gs_main.nrows):
            gs_sub_cov,gs_sub_tx = self.build_cov_tx_grid(fig,gs_main[nr],nr==1)
            gs_subs.append([(gs_sub_cov,gs_sub_tx),[]])
            
            # add data
            axes = self.plot_coverage(fig,gs_sub_cov,title,compare,text_attr,rel,nr==1)
            gs_subs[-1][1].extend(axes)
            axes = self.plot_txs(fig,gs_sub_tx,title,compare,text_attr,rel,nr==1)
            gs_subs[-1][1].extend(axes)
            
            if nr==1:
                # remove axis labels, etc
                for ax in gs_subs[-1][1]:
                    ax.set_xlabel("")
                    ax.set_title("")
                    zoom_start_transform = self.graphcoords[self.settings["zoom_start"]-self.get_start()]
                    zoom_end_transform = self.graphcoords[self.settings["zoom_end"]-self.get_start()]
                    ax.set_xlim(zoom_start_transform,zoom_end_transform)
                    
            
        if self.settings["zoom"]:
            self.build_zoom(fig,gs_subs[0][1],gs_subs[1][1])
            
        return fig
            

    def plot(self,out_fname,title,save_pickle,compare,legend,text_attr,rel):                              
        fig = self.build_fig(out_fname,title,save_pickle,compare,legend,text_attr,rel)
        
        if legend:
            handles = [Patch(color=c[0],label=c[1]) for i,c in self.colors_non_compare.items()]
            if compare:
                handles = [Patch(color=c[0],label=c[1]) for i,c in self.colors_compare.items()]
                handles[1] = Patch(edgecolor=self.colors_compare[1][0],facecolor=None,fill=False,linestyle="-",linewidth=3,label=self.colors_compare[1][1])

            lgd = plt.legend(handles=handles,
                             fontsize=self.settings["font_size"],
                             loc="lower left",
                             bbox_to_anchor=(0., -0.4, 1., .102),
                             ncol=2,
                             mode="expand",
                             borderaxespad=0.)

            if save_pickle:
                with open(out_fname+".pickle","wb") as outFP:
                    pickle.dump((fig,gs1),outFP)

            plt.savefig(out_fname,bbox_extra_artists=(lgd,), bbox_inches='tight')
        else:
            plt.savefig(out_fname, bbox_inches='tight')

# gets minimum and maximum values from BED/bedgraph files
def get_MinMax_val(fname,seqid,strand):
    sj_fname_lst = []
    is_sj_lst_file = True
    if fname is not None:
        with open(fname,"r") as inFP:
            for line in inFP:
                tmp = line.strip()
                if not os.path.exists(tmp) and len(tmp)>0:
                    is_sj_lst_file = False
                    break

        if is_sj_lst_file:
            with open(fname,"r") as inFP:
                for line in inFP:
                    tmp = line.strip()
                    if len(tmp)>0:
                        sj_fname_lst.append(tmp)
        else:
            sj_fname_lst.append(fname)


        intron_starts = []
        intron_ends = []
        for sf in sj_fname_lst:
            with open(sf, "r") as inFP:
                for line in inFP:
                    lcs = line.strip().split("\t")
                    if not len(lcs) == 6:
                        continue

                    if not lcs[0] == seqid:
                        continue

                    if not lcs[5] == strand:
                        continue

                    intron_starts.append(int(lcs[1]))
                    intron_ends.append(int(lcs[2])+1)
    return min(intron_starts),max(intron_ends)

def sashimi(args):
    assert os.path.exists(args.gtf), "GTF does not exist: " + args.gtf
    # assert os.path.exists(args.cov), "Coverage file does not exist: " + args.cov
    # assert os.path.exists(args.sj), "Splice Junction file does not exist: " + args.sj

    settings = {"intron_scale": args.intron_scale,
                "exon_scale": args.exon_scale,
                "number_junctions": args.number_junctions,
                "resolution": args.resolution,
                "fig_width": args.fig_width,
                "tx_height": args.tx_height,
                "cov_height": args.cov_height,
                "reverse": args.reverse,
                "font_size": args.font_size,
                "nxticks": args.nxticks,
                "title": args.title,
                "pickle": args.pickle,
                "compare": args.compare,
                "all_junctions": args.all_junctions,
                "zoom_start":args.zoom_start,
                "zoom_end":args.zoom_end,
                "zoom":args.zoom_start is not None and args.zoom_end is not None}

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
                if args.nocds:
                    tx.set_nocds()
                tx.parse_from_gtf(cur_tid_lines)

                is_ref = args.compare is not None and tx.get_tid()==args.compare
                found_ref = found_ref or is_ref
                locus.add_tx(tx,is_ref)

                cur_tid = tid
                cur_tid_lines = line

            else:
                cur_tid_lines += line

        tx = TX()
        if args.nocds:
            tx.set_nocds()
        tx.parse_from_gtf(cur_tid_lines)
        is_ref = args.compare is not None and tx.get_tid()==args.compare
        found_ref = found_ref or is_ref
        locus.add_tx(tx,is_ref)

    # if requested - read junctions and create a dummy transcript spanning min to max junc/transcript values
    if args.all_junctions and args.sj is not None:
        min_val,max_val = get_MinMax_val(args.sj,locus.seqid,locus.strand)
        dummy_start = min(min_val-1,locus.get_start())
        dummy_end = max(max_val+1,locus.get_end())
        tx = TX()
        tx.set_seqid(locus.seqid)
        tx.set_strand(locus.strand)
        tx.set_exons([(dummy_start,dummy_end)])
        tx.set_orf([])
        tx.set_tid("dummy")
        tx.set_dummy(True)
        locus.add_tx(tx,False)

    locus.set_scaling()

    if args.compare is not None and not found_ref:
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

    # add junctions
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

    # add track names
    if args.tn is not None:
        locus.add_track_names(args.tn)

    title = None
    if not args.title is None:
        title = " ".join(args.title)
    locus.plot(args.output,title,args.pickle,args.compare,args.legend,args.text_attr,args.rel)


def main(args):
    parser = argparse.ArgumentParser(description='''Help Page''')
    parser.add_argument("--gtf",
                        required=True,
                        type=str,
                        help="Annotation in a GFF/GTF format")
    parser.add_argument("--cov",
                        required=False,
                        type=str,
                        help="Coverage in bedgraph format or a file containing a list of filenames with coverage in bedgraph for multiple samples. If a list is provided - the files should be in the same order as the splice junctions below (if provided)")
    parser.add_argument("--sj",
                        required=False,
                        type=str,
                        help="Splice junctions in bed format or a file containing a list of filenames with splice junctions in bed format for multiple samples. If a list is provided - the files should be in the same order as the coverage tracks.")
    parser.add_argument("--tn",
                        required=False,
                        type=str,
                        help="Names for the coverage/junction tracks. If a list is provided - the files should be in the same order as the coverage tracks.")
    parser.add_argument("--rel",
                        action="store_true",
                        required=False,
                        help="Display junciton coverage valus as relative usage of the junction compared to average junciton usage.")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        required=True,
                        help="Filename for the output figure. The format (png,svg, ...) will be automatically deduced based on the extension.")
    parser.add_argument("-c",
                        "--nocds",
                        action="store_true",
                        help="If enabled, will display no CDS features in the output plots.")
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
                        default=5,
                        help="Parameter regulates the smoothing factor of the coverage track (Default: 5). Increasing the value will increasing the smoothing by reducing the number of points on the coverage track.")
    parser.add_argument("--fig_width",
                        required=False,
                        type=int,
                        default=20,
                        help="Width of the figure in inches (Default: 20).")
    parser.add_argument("--tx_height",
                        required=False,
                        type=int,
                        default=2,
                        help="Height of the transcript elements in the figure in inches (Default: 2).")
    parser.add_argument("--cov_height",
                        required=False,
                        type=int,
                        default=3,
                        help="Height of the coverage elements in the figure in inches (Default: 3).")
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
    parser.add_argument("--zoom_start",
                        required=False,
                        type=int,
                        help="Start coordinate to show in the zoomed-in view")
    parser.add_argument("--zoom_end",
                        required=False,
                        type=int,
                        help="End coordinate to show in the zoomed-in view")
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
    parser.add_argument("--text_attr",
                        required=False,
                        type=str,
                        default="transcript_id",
                        help="If specified, the value will determine which attribute key will be used as text displayed for each transcript along with the transcript ID.")
    parser.add_argument("--legend",
                        action="store_true",
                        help="Add legend to the plot")

    parser.set_defaults(func=sashimi)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
    

# TODO: add transcript-level comparison
