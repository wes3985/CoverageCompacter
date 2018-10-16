# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 11:15:36 2018

@author: rmhawwo
"""

import unittest
import timeit

from CoverageCompacter import Loci
from CoverageCompacter import Chromosome
from CoverageCompacter import CoverageCompacter


class Test_CoverageCompacter(unittest.TestCase):
    
    def test_Loci_standard_init(self):
        print("test_Loci_standard_init")
        l = Loci()
        self.assertEqual(l.start, 1)
        self.assertEqual(l.end, 1)
        
    def test_Loci_update_metrics_size(self):
        print("test_Loci_update_metrics_size")
        l = Loci(1,10000)
        l.update_metrics()
        self.assertEqual(l.size, 10000)
        l = Loci(100,10000)
        l.update_metrics()
        self.assertEqual(l.size, 9901)
        
    def test_Loci_update_metrics_meanCoverage(self):
        print("test_Loci_update_metrics_meanCoverage")
        l = Loci(1,10000)
        l.depthSum = 10000
        l.update_metrics()
        self.assertEqual(l.meanCoverage, 1)
        l.depthSum = 20000
        l.update_metrics()
        self.assertEqual(l.meanCoverage, 2)
        l.depthSum = 15000
        l.update_metrics()
        self.assertEqual(l.meanCoverage, 1.5)
        l.depthSum = 0
        l.update_metrics()
        self.assertEqual(l.meanCoverage, 0)
        
    def test_Loci_update_metrics_coverageFraction(self):
        print("test_Loci_update_metrics_coverageFraction")
        l = Loci(1,10000)
        l.N_basesCovered = 10000
        l.depthSum = 10000
        l.update_metrics()
        self.assertEqual(l.coverageFraction, 1.0)
        l.N_basesCovered = 5000
        l.depthSum = 10000
        l.update_metrics()
        self.assertEqual(l.coverageFraction, 0.5)
        l.N_basesCovered = 0
        l.depthSum = 0
        l.update_metrics()
        self.assertEqual(l.coverageFraction, 0)
        l.N_basesCovered = 5000
        l.depthSum = 0
        l.update_metrics()
        self.assertEqual(l.coverageFraction, 0.5)            # This test checks that depthSum doesn't affect the coverageFraction, though depthSum should not be 0 if N_basesCovered >0

    
    def test_Chromosome(self):
        print("test_Chromosome")
        c = Chromosome('testSample', 10000, 'chr1')            
        self.assertEqual(c.sample,"testSample")
        self.assertEqual(c.CHR, "chr1")
        self.assertEqual(c.binEdge,5000)
        self.assertEqual(c.current_loci.start, 1)
        self.assertEqual(c.current_loci.end, 1)
        
    def test_Chromosome_update_current_loci(self):
        print("test_Chromosome_update_current_loci")
        c = Chromosome('testSample', 10, 'chr1')
        c.update_current_loci(1)                            # add a base with depth 1
        self.assertEqual(c.current_loci.start, 1)
        self.assertEqual(c.current_loci.end, 2)
        c.update_current_loci(5)                            # add a base with depth 5
        self.assertEqual(c.current_loci.start, 1)
        self.assertEqual(c.current_loci.end, 3)
        c.update_current_loci(0)                            # add a base with depth 0
        self.assertEqual(c.current_loci.start, 1)
        self.assertEqual(c.current_loci.end, 4)
        
    def test_Chromosome_zeroBin(self):
        print("test_Chromosome_zeroBin")
        c = Chromosome('testSample', 0, 'chr1')            
        self.assertEqual(c.sample,"testSample")
        self.assertEqual(c.CHR, "chr1")
        self.assertEqual(c.binEdge,0)
        self.assertEqual(c.current_loci.start, 1)
        self.assertEqual(c.current_loci.end, 1)
        
    def test_Chromosome_makes_empty_loci(self):
        print("Chromosome_makes_empty_loci")
        c = Chromosome('testSample', 0, 'chr1')
        c.update_current_loci(0)  
        self.assertEqual(c.current_loci.start, 1)
        self.assertEqual(c.current_loci.end, 2)
        c.update_current_loci(0)  
        self.assertEqual(c.current_loci.start, 1)
        self.assertEqual(c.current_loci.end, 3)
        c.update_current_loci(0)  
        self.assertEqual(c.current_loci.start, 1)
        self.assertEqual(c.current_loci.end, 4)
        c.update_current_loci(1)  
        self.assertEqual(c.LociList[-1].end, 3)             # check the last loci end has been adjusted backwards
        self.assertEqual(c.current_loci.start, 4)
        self.assertEqual(c.current_loci.end, 5)
        
    def test_Chromosome_zeroBin_update_current_loci(self):
        print("test_Chromosome_zeroBin_update_current_loci")
        c = Chromosome('testSample', 0, 'chr1')
        c.update_current_loci(1)                            # add a base with depth 1
        self.assertEqual(c.current_loci.start, 1)
        self.assertEqual(c.current_loci.end, 2)
        c.update_current_loci(5)                            # add a base with depth 5
        self.assertEqual(c.current_loci.start, 1)
        self.assertEqual(c.current_loci.end, 3)
        c.update_current_loci(0)                            # add a base with depth 0
        self.assertEqual(c.LociList[-1].end, 2)             # check the last loci end has been adjusted backwards
        self.assertEqual(c.current_loci.start, 3)
        self.assertEqual(c.current_loci.end, 4)
        c.update_current_loci(0)                            # add another base with depth 0
        self.assertEqual(c.current_loci.start, 3)
        self.assertEqual(c.current_loci.end, 5)   
        
    def test_Chromosome_binsize1_update_current_loci(self):
        print("test_Chromosome_binsize1_update_current_loci")
        c = Chromosome('testSample', 1, 'chr1')
        c.update_current_loci(0)   
        c.update_current_loci(1)
        self.assertEqual(c.current_loci.start, 1)
        self.assertEqual(c.current_loci.end, 3)     
        
    def test_Chromosome_binsize1_binBoundaries(self):
        print("test_Chromosome_binsize1_binBoundaries")
        c = Chromosome('testSample', 1, 'chr1')
        for x in range(5):
            c.update_current_loci(2)
        c.update_current_loci(0)
        for x in range(5):
            c.update_current_loci(2)  
        self.assertEqual(c.current_loci.start, 1)
        self.assertEqual(c.current_loci.end, 12)  
        c.make_tidy()
        self.assertEqual(c.current_loci.end, 11)
        
    def test_Chromosome_binsize2_binBoundaries(self):
        print("test_Chromosome_binsize2_binBoundaries")
        c = Chromosome('testSample', 2, 'chr1')
        self.assertEqual(c.binEdge, 1)
        for x in [1,0,0,1]:
            #print(c.current_loci.start, c.current_loci.end)
            c.update_current_loci(x)
        for p in c.LociList:
            print(p.start, p.end, p.LastCoveredBase)
        print(c.current_loci.start, c.current_loci.end, c.current_loci.LastCoveredBase)
        self.assertEqual(c.current_loci.end, 5)
#        self.assertEqual(c.
        self.assertEqual(c.current_loci.start, 1)
        self.assertEqual(len(c.LociList), 0)
        
        
    def test_Chromosome_binsize3_binBoundaries(self):
        print("test_Chromosome_binsize3_binBoundaries")
        c = Chromosome('testSample', 3, 'chr1')
        self.assertEqual(len(c.LociList), 0)
        self.assertEqual(c.binEdge, 2)
        for y,x in enumerate([0,1,0,0,0],1):
            print(y, x); c.update_current_loci(x)
        # CHECKS BEFORE MERGE
        self.assertEqual(c.current_loci.N_basesCovered, 0)
        self.assertEqual(c.current_loci.end, 6)
        self.assertEqual(c.current_loci.start, 5)
        self.assertEqual(len(c.LociList), 1)
        self.assertEqual(c.LociList[-1].LastCoveredBase, 2)
        c.update_current_loci(1)
#        # CHECKS AFTER MERGE
        self.assertEqual(len(c.LociList), 0)
        self.assertEqual(c.currentPosition, 7)
        self.assertEqual(c.current_loci.start, 1)
        self.assertEqual(c.current_loci.end, 7)
        self.assertEqual(c.current_loci.N_basesCovered, 2)
        self.assertEqual(c.current_loci.depthSum, 2)
        self.assertEqual(c.current_loci.FirstCoveredBase, 2)
        self.assertEqual(c.current_loci.LastCoveredBase, 6)


    def test_Chromosome_binsize3_binEdgeBegins(self):
        print("test_Chromosome_binsize3_binEdgeBegins")
        c = Chromosome('testSample', 3, 'chr1')
        for y,x in enumerate([0,0,1,1], 1):
            print(y); c.update_current_loci(x)
        self.assertEqual(len(c.LociList), 0)
        
    def test_Chromosome_binsize3_binEdgeOutOfBounds(self):
        print("test_Chromosome_binsize3_binEdgeOutOfBounds")
        c = Chromosome('testSample', 3, 'chr1')
        for y,x in enumerate([0,0,0,0,1,1], 1):
            print(y); c.update_current_loci(x)
        self.assertEqual(len(c.LociList), 1)
        self.assertEqual(c.LociList[-1].start, 1)
        self.assertEqual(c.LociList[-1].end, 2)
        self.assertEqual(c.current_loci.start, 3)
        self.assertEqual(c.current_loci.end,7)
        
    def test_Chromosome_binsize3_multibin(self):
        print("test_Chromosome_binsize3_multibin")
        c = Chromosome('testSample', 3, 'chr1')
        for y,x in enumerate([0,0,1,1,1,0,1,1,1,0,0,0,1,1,0,1,1,1,0],1):
            print(y); c.update_current_loci(x)
        self.assertEqual(len(c.LociList), 0)
        
    def test_Chromosome_reDefine_binBoundary_if_no_coverage_more_than_binEdge_but_less_than_binSize(self):
        print("test_Chromosome_reDefine_binBoundary_if_no_coverage_more_than_binEdge_but_less_than_binSize")
        c = Chromosome('testSample', 5, 'chr1')
        for y,x in enumerate([1,1,0,0,0,0,0,0,1,1],1):
            print(y); c.update_current_loci(x)
        self.assertEqual(len(c.LociList), 1)
        self.assertEqual(c.LociList[-1].start, 1)
        self.assertEqual(c.LociList[-1].end, 5)
        self.assertEqual(c.current_loci.start, 6)
        self.assertEqual(c.current_loci.end, 11)
            
    def test_CoverageCompacter(self):    
        dep_f="coverage_test_file_100k.txt"
        CoverageCompacter(dep_f, 'CBATest', 'AE,BE,HE,C', 'chr1', 0)

# TIME WHOLE FUNCTION
def time_test():
    dep_f="coverage_test_file_100k.txt"
    t=timeit.Timer(lambda: CoverageCompacter(dep_f, 'CBATest', 'AE,BE,HE,C', 'chr1', 0)).timeit(10)
    print('Time for 10 runs: ',t)
    
if __name__ == '__main__':
    unittest.main()
    #time_test()