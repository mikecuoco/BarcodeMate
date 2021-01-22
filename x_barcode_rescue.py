import os
import sys
import pysam


##Why this module?
##In NA12878, there are 865,601,723 reads in all, and 47,523,635 (~5%) do not have "BX" tag ("RX is the original"),
## out of them, 41829412 are with mapping quality >=30, and  38595972 with mapping quality >=60


#for single cell RNA-seq data, some reads do not have tag: "CB" (where "CR" is the original sequencer reported one)
# 117,104 out of 8,719,155 reads fall in this range in one test. Seems not that much

class Shift_OR():
    def __init__(self, sf_text):
        self.sf_text=sf_text
    def match(self, spattern):
        return



class XBarcodeRescue():
    def __init__(self, sf_bam, s_working_folder, n_jobs):
        self.sf_bam=sf_bam
        self.working_folder=s_working_folder
        self.n_jobs=n_jobs

    def rescue_by_region(self, record):
        return

    ##for those unique mapped reads, but no "BX" tag,
    #1) first, we collect all the barcode within 100kb region
    #2) we compare the "RX" tag again these collected barcodes
    #3) pick the best one if multiple hits found
    def rescue_barcode(self, sf_out_bam):
        return

    def shift_or(self):
        return

    ##collect all the barcodes of a region
    def collect_barcode_in_region(self, chrm, istart, iend):
        set_barcode = set()
        samfile = pysam.AlignmentFile(self.sf_bam, "rb")
        for alignmt in samfile.fetch(chrm, istart, iend):
            if alignmt.has_tag("BX"):
                s_barcode = alignmt.get_tag("BX")
                set_barcode.add(s_barcode)
        samfile.close()


####should not be aligned here
####