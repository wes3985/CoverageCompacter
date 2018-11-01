# -*- coding: utf-8 -*-
"""
@author: w.woollard

"""

import sys
from collections import OrderedDict
import math
import datetime

class Loci(object):

    def __init__(self, start=1, end=1):
        self.start = start                      # start of loci
        self.end = end                          # end of loci
        self.N_basesCovered = 0                 # total number of bases with coverage
        self.depthSum = 0                       # sum of all depth
        self.LastCoveredBase = None             # the last base with coverage >0
        self.FirstCoveredBase = None            # first base covered >0
        
        self.size = 1                           # size of loci
        self.meanCoverage = 0.0                 # the mean coverage of the loci
        self.coverageFraction = 0.0             # fraction of the loci that has coverage >0
        
    def update_metrics(self):
        if self.end < self.start:
            self.end = self.start
        self.size = (self.end - self.start) +1
        self.start = int(self.start)
        self.end = int(self.end)
        try:
            self.meanCoverage = self.depthSum / self.size
            self.coverageFraction = self.N_basesCovered / self.size
        except ZeroDivisionError:
            self.meanCoverage = 0.0
            self.coverageFraction = 0.0

    def increment_depth(self, depth, currentPosition):
        # this method is invoked whenever a covered base is read
        self.LastCoveredBase = currentPosition
        self.N_basesCovered += 1
        self.depthSum += depth  
            
        
       
class Chromosome(object):

    def __init__(self, sample, binSize, CHR, NoCov=0):
        self.sample = sample                            # name of the sample that this loci belongs to
        self.binSize = binSize                          # the total size that a covered region may extend to OR the max size of non-covered regions between covered regions
        self.binEdge = math.ceil(binSize / 2)           # the possible extended range of the amplicon in bp beyond the last covered base
        self.CHR = CHR                                  # name of chromosom that the loci is on 
        self.NoCov = NoCov                              # The minimum depth to consider as having acceptable coverage
        self.currentPosition = 1                        # the current position in the chromosome that the generator is pointing to    
        self.LociList = []                              # list to hold all loci in the chromosome
        self.current_loci = self.initiate_loci(1,1)     # current chunk of chromosome 
        
    def initiate_loci(self, start, end):
        return Loci(start, end)      
    
    # called at the end of the chromsome to finish the current loci
    def make_tidy(self):
        self.current_loci.end -=1                          
        self.current_loci.update_metrics()
        self.LociList.append(self.current_loci)    
        
    def merge_loci(self):
        self.current_loci.start = self.LociList[-1].start
        if self.LociList[-1].N_basesCovered != None:
            self.current_loci.N_basesCovered += self.LociList[-1].N_basesCovered
        self.current_loci.depthSum += self.LociList[-1].depthSum
        if self.LociList[-1].FirstCoveredBase != None:
            self.current_loci.FirstCoveredBase = self.LociList[-1].FirstCoveredBase
            self.current_loci.LastCoveredBase = self.LociList[-1].LastCoveredBase
        self.LociList.pop()
        if len(self.LociList) != 0:
            if self.LociList[-1].LastCoveredBase != None:
                if ((self.currentPosition-1) - self.LociList[-1].LastCoveredBase) <= self.binSize:
                    # This nested merge is a failsafe incase merging has failed to comlete properly, 
                    # it doesn't appear to be called in any simulations
                    self.current_loci = self.merge_loci
                    print('WARNING: nested merge occured at:',self.CHR,self.currentPosition)
        return self.current_loci
    

    
        
    def update_current_loci(self, depth):
       
        if self.current_loci.N_basesCovered == 0 and depth > self.NoCov:            
           
            if self.current_loci.end - self.current_loci.start >= self.binEdge:                
                
                if len(self.LociList) > 0 and \
                    ((self.currentPosition-1) - self.LociList[-1].LastCoveredBase) <= self.binSize:
                    # Merge bins if current and previous edges overlap 
                    self.current_loci = self.merge_loci()
                    self.current_loci.increment_depth(depth, self.currentPosition)         
                    
                else:
                # if current and previous bin edges do not overlap then a loci with no coverage
                # must be added to the list that sits between 2 loci with coverage
                    self.current_loci.end = (self.currentPosition - (self.binEdge+1))                                       # reset this loci end to the amplicon boundary
                    self.current_loci.update_metrics()
                    self.LociList.append(self.current_loci)                                                                 # add the empty loci to the list                                         
                    self.current_loci = self.initiate_loci((self.currentPosition - self.binEdge), self.currentPosition)     # intitate a new loci with bin edges adjusted backwards
                    self.current_loci.FirstCoveredBase = self.currentPosition
                    if len(self.LociList) > 0:
                        if ((self.currentPosition) - self.LociList[-1].start) <= self.binEdge:
                            # final check just incase loci can be merged
                            self.current_loci = self.merge_loci()
                    self.current_loci.increment_depth(depth, self.currentPosition) 

            else:
                # Current loci has no coverage so far but is small enough to be within range of the bin edge
                if len(self.LociList) > 0 and \
                    ((self.currentPosition-1) - self.LociList[-1].LastCoveredBase) <= self.binSize:
                    # merge bins if current and previous edges overlap 
                    self.current_loci = self.merge_loci()
                    self.current_loci.increment_depth(depth, self.currentPosition)
          
                     
                else:
                # if there is no previous bin, or the last covered base is to far away to be
                # within range of the binSize and so can't be merged then add coverage to this loci
                    self.current_loci.FirstCoveredBase = self.currentPosition
                    self.current_loci.increment_depth(depth, self.currentPosition)                    
                
        elif self.current_loci.N_basesCovered != 0 and depth <= self.NoCov:
            # if the loci has accumilated a large enough region of no coverage
            # so that the bin range has been reached then append to the list and initiate a new loci
            if self.currentPosition - self.current_loci.LastCoveredBase > self.binEdge:
                self.current_loci.end = self.currentPosition -1
                self.current_loci.update_metrics()
                self.LociList.append(self.current_loci)  # add the empty loci to the list
                self.current_loci = self.initiate_loci(self.LociList[-1].end +1,self.LociList[-1].end +1)     # intitate a new loci
                
            else:
            # the bin range hasn't been reached yet so just add a non-covered base to the current loci
                pass
                
        elif self.current_loci.N_basesCovered != 0 and depth > self.NoCov:
            # current loci has coverage and will be extended by new coverage
            self.current_loci.increment_depth(depth, self.currentPosition)           
            
        else:
        # if there is no reason to start a new loci then add an empty base
            pass

        self.currentPosition += 1    
        self.current_loci.end = self.currentPosition 

# MAIN FUNCTION
def CoverageCompacter(depth_file='', outfilePrefix='', samples='', CHROM='', binSize=10000, NoCov=0):
    binSize=int(binSize)
    NoCov=int(NoCov)
    SAMPLES=samples.split(',')
    s_dict = OrderedDict()
    c_dict = OrderedDict()
    for i,s in enumerate(SAMPLES,2):
        s_dict[s]=i                                                   # store position in each row for each sample in an orderedict
        c_dict[i]=Chromosome(s,binSize,CHROM, NoCov)                  # initiate chromosomes for each sample
    
    # read in file and update chromsomes for each sample per line of depth file read
    with open(depth_file, 'r') as fh:
        depth_lines = get_data(fh)
        while True:
            try:
                d=next(depth_lines)
                for s in c_dict.keys():
                    c_dict[s].update_current_loci(int(d[s]))          # pass the depth data to be processed for each sample                               
            except:
                StopIteration  
                # at end of file remove the last base and update_metrics and append to LociList
                for s in c_dict:
                    c_dict[s].make_tidy()           
                break 
    
    for i,s in zip(c_dict,s_dict):
        with open(outfilePrefix+'_'+s+'.txt','w') as w:
            w.write('chr\tstart\tend\tsize\tfirstCoveredBase\tlastCoveredBase\tmeanCoverage\tNBasesCovered\tDepthSum\tcoverageFraction\n')
            for bed in c_dict[i].LociList:
                w.write('\t'.join([str(x) for x in [c_dict[i].CHR, bed.start, bed.end, bed.size, bed.FirstCoveredBase, bed.LastCoveredBase, bed.meanCoverage, bed.N_basesCovered, bed.depthSum, bed.coverageFraction]])+'\n')        
        
    return

# GENERATOR FUNCTION
def get_data(file_handle):
    for line in file_handle:
        yield line.split()     # line structure from samtools depth (files have no header) [CHR, position, sample1, .....sampleN]


if __name__ == '__main__':
    print("start: ",datetime.datetime.now())
    CoverageCompacter(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    print("finished: ",datetime.datetime.now())    
