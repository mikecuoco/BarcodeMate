##11/22/2017
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

import os
import global_values
from x_converter import *
from x_merger import *
from x_trimmer import *
from x_hap_block import *
from x_matrix import *
from x_barcode_rescue import *
from x_statistic import *
from optparse import OptionParser

##Todo List: 1. barcode trimming algorithm (consider indels)
##           2. molecule (barcode) coverage bed file
##           3. barcode matrix
##           4. molecule length statistic ("-S" option), the right boundary need to check the cigar to get exact value

####
##parse the options
def parse_option():
    parser = OptionParser()
    parser.add_option("-C", "--convert",
                      action="store_true", dest="convert", default=False,
                      help="convert the bam to barcode indexed bam")
    parser.add_option("-X", "--molecule",
                      action="store_true", dest="molecule", default=False,
                      help="convert the bam to cell barcode indexed bam")
    parser.add_option("-K", "--rconvert",
                      action="store_true", dest="rconvert", default=False,
                      help="convert the barcode indexed bam back to normal bam")
    parser.add_option("-M", "--merge",
                      action="store_true", dest="merge", default=False,
                      help="merge 10X bams")
    parser.add_option("-S", "--statistic",
                      action="store_true", dest="statistic", default=False,
                      help="basic statistic from 10X bam")
    parser.add_option("-B", "--hap_block",
                      action="store_true", dest="block", default=False,
                      help="Haplotype block information")
    parser.add_option("-R", "--barcode_rescue",
                      action="store_true", dest="barcode_rescue", default=False,
                      help="Rescue reads without BX tags")
    parser.add_option("-G", "--matrix",
                      action="store_true", dest="matrix", default=False,
                      help="Generate barcode matrix")
    parser.add_option("-T", "--trim",
                      action="store_true", dest="trim", default=False,
                      help="Trim the clipped reads caused by errors from 10X bam")
    parser.add_option("-i", "--input", dest="input",
                      help="input file) ", metavar="FILE")
    parser.add_option("-r", "--reference", dest="reference",
                      help="The reference file) ", metavar="FILE")
    parser.add_option("-b", "--bam", dest="bam",
                      help="Input bam file", metavar="FILE")
    parser.add_option("-o", "--output", dest="output",
                      help="The output file", metavar="FILE")
    parser.add_option("-p", "--path", dest="wfolder", type="string",
                      help="Working folder")
    parser.add_option("-w", "--bin", dest="bin", type="int",
                      help="Bin size")
    parser.add_option("-k", "--nbin", dest="nbin", type="int",
                      help="# of bins per block")
    parser.add_option("-m", "--min_cnt", dest="min", type="int",
                      help="Minimum number of shared barcodes to save")
    parser.add_option("-n", "--cores", dest="cores", type="int",
                      help="Number of cores")

    (options, args) = parser.parse_args()
    return (options, args)

####
##main function
if __name__ == '__main__':
    (options, args) = parse_option()
    if options.convert:
        sf_ori_bam = options.bam
        sf_barcode_bam = options.output
        s_working_folder = options.wfolder
        n_jobs = options.cores
        xcvtr = XConverter(sf_ori_bam, sf_barcode_bam, n_jobs)
        xcvtr.set_working_folder(s_working_folder)
        xcvtr.cvt_bam_to_barcode_bam()

    if options.molecule:#index by "MI" barcode information
        global_values.change_to_molecule_barcode()#use "MI" as barcode
        sf_ori_bam = options.bam
        sf_barcode_bam = options.output
        s_working_folder = options.wfolder
        n_jobs = options.cores
        xcvtr = XConverter(sf_ori_bam, sf_barcode_bam, n_jobs)
        xcvtr.set_working_folder(s_working_folder)
        xcvtr.cvt_bam_to_barcode_bam()

    elif options.barcode_rescue:
        sf_bam = options.bam
        s_working_folder = options.wfolder
        n_jobs = options.cores
        sf_rescued_bam = options.output
        xbx_rsc = XBarcodeRescue(sf_bam, s_working_folder, n_jobs)
        xbx_rsc.rescue_barcode(sf_rescued_bam)

    elif options.statistic:#statistic for molecule length
        sf_bam = options.bam ##this is barcode (BX) based indexed bam
        s_working_folder = options.wfolder
        n_jobs = options.cores
        sf_rslt = options.output
        #algorithms
        x_molecule=MoleculeInfo(sf_bam, s_working_folder, n_jobs)
        x_molecule.get_molecule_size_BX_bam(sf_rslt)
        #x_molecule.statistic_molecule_size(sf_rslt)

    elif options.rconvert:
        sf_barcode_bam = options.bam
        sf_head_bam=options.input
        sf_ori_bam = options.output
        n_jobs = 1 ##this module is not parallelized
        xcvtr = XConverter(sf_ori_bam, sf_barcode_bam, n_jobs)
        xcvtr.cvt_barcode_bam_to_normal(sf_head_bam)

    elif options.matrix:
        sf_bam = options.bam
        n_bins_block = options.nbin #number of bins per block
        bin_size = options.bin #size of each bin in bp
        i_barcode_cutoff = 5000 #maximum number of barcode per block
        sf_working_folder = options.wfolder
        n_cores = options.cores
        i_min_freq_cutoff = options.min #minimum number of shared barcodes between two blocks
        sf_chrm_size=options.input
        sf_out = options.output
        xmatrix = XMatrix(sf_bam, n_bins_block, bin_size, i_barcode_cutoff, sf_working_folder, i_min_freq_cutoff)
        xmatrix.cnt_shared_barcodes(sf_chrm_size, n_cores, sf_out)

    elif options.merge:
        sf_file_list= options.input
        sf_head_bam=options.bam
        sf_merged=options.output
        s_tmp_folder = options.wfolder
        n_jobs=options.cores

        bam_merger=BamMerger(sf_file_list, sf_head_bam, n_jobs, s_tmp_folder, sf_merged)
        bam_merger.merge_bam()
        #bam_merger.sort_index_bam_with_sambamba(sf_merged, s_tmp_folder)

    elif options.trim:
        sf_input_bam=options.bam
        sf_ref=options.reference
        s_tmp_folder = options.wfolder
        n_jobs = options.cores
        sf_trimmed_bam=options.output
        rt = ReadTrimmer(sf_input_bam, n_jobs)
        rt.trim_reads(s_tmp_folder, sf_ref, sf_trimmed_bam)

    elif options.block:
        sf_input_bam = options.bam
        s_tmp_folder = options.wfolder
        n_jobs = options.cores
        s_out=options.output

        xblock = XHapBlock(sf_input_bam, s_tmp_folder, n_jobs)
        xblock.get_block_info(s_out)

##TDList:
##1. given a region, get out how many barcodes inside
##2. given a site, get out how many barcodes cover
##3. sambamba sort -m 50G --tmpdir=./tmp

##
##