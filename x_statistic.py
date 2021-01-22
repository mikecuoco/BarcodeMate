##10/16/2018
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu


import os
import sys
import pysam
from subprocess import *
from multiprocessing import Pool

AVE_MOLECULE_SIZE=100000 #suppose average molecule size if 100K

def unwrap_self_molecule_statistic_by_MI(arg, **kwarg):
    return MoleculeInfo.run_molecule_statistic_for_MIs(*arg, **kwarg)

def unwrap_self_molecule_statistic_by_BX_MI(arg, **kwarg):
    return MoleculeInfo.run_molecule_statistic_for_BX_MI(*arg, **kwarg)
#

class MoleculeInfo():
    def __init__(self, sf_bam, working_folder, n_jobs):
        self._sf_bam = sf_bam
        self._working_folder = working_folder
        self._n_jobs = n_jobs

    def run_molecule_statistic_for_MIs(self, record):
        istart=int(record[0][0])
        iend=int(record[0][1])
        sf_bam=record[1]
        sf_tmp_out=record[2]

        samfile = pysam.AlignmentFile(sf_bam, "rb")
        barcodes = samfile.references

        m_statistic={}
        for s_mi in barcodes[istart:iend+1]:
            m_chrs={}
            l_algnmts=[]
            for alignmt in samfile.fetch(s_mi):
                ch = "*"
                if alignmt.has_tag("ch"):
                    ch = alignmt.get_tag("ch")
                if ch not in m_chrs:
                    m_chrs[ch]=1
                else:
                    m_chrs[ch] += 1

                cp = 0
                if alignmt.has_tag("cp"):
                    cp = alignmt.get_tag("cp")
                l_algnmts.append((ch, cp))

            max_chrm, i_left_pos, i_right_pos=self.statistic_one_molecule(m_chrs, l_algnmts)
            if max_chrm=="*":
                continue
            m_statistic[s_mi]=(i_left_pos, i_right_pos)
        samfile.close()

        with open(sf_tmp_out, "w") as fout_tmp:
            for s_mi in m_statistic:
                i_size=m_statistic[s_mi][1]-m_statistic[s_mi][0]
                s_rcd="{0}\t{1}\t{2}\t{3}\n".format(s_mi, m_statistic[s_mi][0], m_statistic[s_mi][1], i_size)
                fout_tmp.write(s_rcd)


    def run_molecule_statistic_for_BX_MI(self, record):
        istart=int(record[0][0])
        iend=int(record[0][1])
        sf_bam=record[1]
        sf_tmp_out=record[2]

        samfile = pysam.AlignmentFile(sf_bam, "rb")
        barcodes = samfile.references
        #print barcodes[0], istart, iend ###################################################
        with open(sf_tmp_out, "w") as fout_tmp:
            for s_bx in barcodes[istart:iend+1]:
                m_chrs={}
                m_algnmts={}
                #m_statistic = {} #for each BX barcode
                for alignmt in samfile.fetch(s_bx):#all the reads of the same BX barcode
                    ####each BX should have several MI barcode
                    #
                    s_mi="null"
                    if alignmt.has_tag("MI"):
                        s_mi=str(alignmt.get_tag("MI"))
                    if s_mi=="null":#skip the alignments without a "MI" tag
                        continue
                    #
                    if s_mi not in m_chrs:
                        m_chrs[s_mi]={}
                    if s_mi not in m_algnmts:
                        m_algnmts[s_mi]=[]

                    ch = "*"
                    if alignmt.has_tag("ch"):
                        ch = alignmt.get_tag("ch")
                    if ch not in m_chrs[s_mi]:
                        m_chrs[s_mi][ch]=1
                    else:
                        m_chrs[s_mi][ch] += 1

                    cp = 0
                    if alignmt.has_tag("cp"):
                        cp = alignmt.get_tag("cp")
                    m_algnmts[s_mi].append((ch, cp))

                for s_mi in m_algnmts:
                    if s_mi not in m_chrs:
                        continue

                    max_chrm, i_left_pos, i_right_pos=self.statistic_one_molecule(m_chrs[s_mi], m_algnmts[s_mi])
                    if max_chrm=="*":
                        continue
                    # if s_mi not in m_statistic:
                    #     m_statistic[s_mi]=[]
                    #m_statistic[s_mi].append((s_bx, i_left_pos, i_right_pos))
                    i_size=i_right_pos-i_left_pos
                    n_reads=len(m_algnmts[s_mi])
                    s_rcd = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(s_mi, s_bx, max_chrm, i_left_pos,
                                                                         i_right_pos, i_size, n_reads)
                    fout_tmp.write(s_rcd)
                #write to file

        samfile.close()

        # with open(sf_tmp_out, "w") as fout_tmp:
        #     for s_mi in m_statistic:
        #         for rcd in m_statistic[s_mi]:
        #             s_bx=rcd[0]
        #             i_left=rcd[1]
        #             i_right=rcd[2]
        #             i_size=i_right-i_left
        #             s_rcd="{0}\t{1}\t{2}\t{3}\t{4}\n".format(s_mi, s_bx, i_left, i_right, i_size)
        #             fout_tmp.write(s_rcd)
#

    def statistic_one_molecule(self, m_chrs, l_algnmts):
        i_max = 0
        max_chrm = ""
        for tmp_chrm in m_chrs:
            if m_chrs[tmp_chrm] > i_max:
                i_max = m_chrs[tmp_chrm]
                max_chrm = tmp_chrm
        if max_chrm == "*":
            #continue
            return max_chrm, -1, -1

        # find cluster, and then the start and end position of the cluster
        l_tmp_pos = []
        for tmp_rcd in l_algnmts:
            if tmp_rcd[0] == max_chrm:
                l_tmp_pos.append(int(tmp_rcd[1]))
        l_tmp_pos.sort()

        n_size = len(l_tmp_pos)
        i_medium = n_size / 2
        medium_pos = l_tmp_pos[i_medium]
        i_left_pos = l_tmp_pos[0]
        i_right_pos = l_tmp_pos[-1]
        for i in range(n_size):
            if abs(l_tmp_pos[i] - medium_pos) < AVE_MOLECULE_SIZE:
                i_left_pos = int(l_tmp_pos[i])
                break
        for tmp_pos in reversed(l_tmp_pos):
            if abs(tmp_pos - medium_pos) < AVE_MOLECULE_SIZE:
                i_right_pos = tmp_pos
                break
        #m_statistic[s_mi] = (i_left_pos, i_right_pos)
        return max_chrm, i_left_pos, i_right_pos

    def get_molecule_size(self, sf_out):
        #first get all the MI from the head
        samfile = pysam.AlignmentFile(self._sf_bam, "rb")
        barcodes = samfile.references
        n_barcodes=len(barcodes)
        l_regions = []
        n_ave = n_barcodes / self._n_jobs
        istart = 0
        iend = 0
        while istart < (self._n_jobs * n_ave):
            iend = istart + n_ave
            l_regions.append((istart, iend - 1))
            istart = iend
        if iend < (n_barcodes - 1):
            l_regions.append((iend, n_barcodes - 1))
        samfile.close()
        #
        l_records=[]
        for rcd in l_regions:
            sf_tmp_out = self._working_folder + "{0}_{1}.molecule_stat.txt".format(rcd[0], rcd[1])
            l_records.append((rcd, self._sf_bam, sf_tmp_out))

        pool = Pool(self._n_jobs)
        pool.map(unwrap_self_molecule_statistic_by_MI, zip([self] * len(l_records), l_records), 1)
        pool.close()
        pool.join()

        #merge the tmp results:
        with open(sf_out, "w") as fout_rslt:
            fout_rslt.write("#MI\tChrom\tMolecule-Start\tMolecule-End\tMolecule-length\n")
            for rcd in l_regions:
                sf_tmp_out = self._working_folder + "{0}_{1}.molecule_stat.txt".format(rcd[0], rcd[1])
                if os.path.isfile(sf_tmp_out)==True:
                    with open(sf_tmp_out) as fin_tmp:
                        for line in fin_tmp:
                            fout_rslt.write(line)

####
    #This is assume the bam is barcode (BX) based, and already sorted
    def get_molecule_size_BX_bam(self, sf_out):
        # first get all the MI from the head
        samfile = pysam.AlignmentFile(self._sf_bam, "rb")
        barcodes = samfile.references
        n_barcodes = len(barcodes)
        l_regions = []
        n_ave = n_barcodes / self._n_jobs
        istart = 0
        iend = 0
        while istart < (self._n_jobs * n_ave):
            iend = istart + n_ave
            l_regions.append((istart, iend - 1))
            istart = iend
        if iend < (n_barcodes - 1):
            l_regions.append((iend, n_barcodes - 1))
        samfile.close()
        #
        l_records = []
        for rcd in l_regions:
            sf_tmp_out = self._working_folder + "{0}_{1}.molecule_stat.txt".format(rcd[0], rcd[1])
            l_records.append((rcd, self._sf_bam, sf_tmp_out))
            #record=(rcd, self._sf_bam, sf_tmp_out)
            #self.run_molecule_statistic_for_BX_MI(record)

        pool = Pool(self._n_jobs)
        pool.map(unwrap_self_molecule_statistic_by_BX_MI, zip([self] * len(l_records), l_records), 1)
        pool.close()
        pool.join()

        # merge the tmp results:
        with open(sf_out, "w") as fout_rslt:
            fout_rslt.write("#MI\tBX\tChrm\tMolecule-Start\tMolecule-End\tMolecule-length\tNum-of-reads\n")
            for rcd in l_regions:
                sf_tmp_out = self._working_folder + "{0}_{1}.molecule_stat.txt".format(rcd[0], rcd[1])
                if os.path.isfile(sf_tmp_out) == True:
                    with open(sf_tmp_out) as fin_tmp:
                        for line in fin_tmp:
                            fout_rslt.write(line)


    def get_molecule_size_single_core(self, sf_out):
        samfile = pysam.AlignmentFile(self._sf_bam, "rb")
        m_chrs = {}
        l_algnmts = []
        pre_mi=None
        with open(sf_out,"w") as fout_rslt:
            fout_rslt.write("#MI\tMolecule-Start\tMolecule-End\tMolecule-length\n")
            for alignmt in samfile.fetch(until_eof=True):#this is retrieve all the reads, and no need to be sorted
                cur_mi=alignmt.reference_id
                if pre_mi!=None and pre_mi!=cur_mi:
                    max_chrm, i_left_pos, i_right_pos = self.statistic_one_molecule(m_chrs, l_algnmts)
                    if max_chrm!="*":
                        i_size = i_right_pos - i_left_pos
                        s_rcd = "{0}\t{1}\t{2}\t{3}\n".format(pre_mi, i_left_pos, i_right_pos, i_size)
                        fout_rslt.write(s_rcd)
                    m_chrs.clear()
                    del l_algnmts[:]
                pre_mi=cur_mi

                ch = "*"
                if alignmt.has_tag("ch"):
                    ch = alignmt.get_tag("ch")
                if ch not in m_chrs:
                    m_chrs[ch] = 1
                else:
                    m_chrs[ch] += 1

                cp = 0
                if alignmt.has_tag("cp"):
                    cp = alignmt.get_tag("cp")
                l_algnmts.append((ch, cp))
            #save the last one
            max_chrm, i_left_pos, i_right_pos = self.statistic_one_molecule(m_chrs, l_algnmts)
            if max_chrm != "*":
                i_size = i_right_pos - i_left_pos
                s_rcd = "{0}\t{1}\t{2}\t{3}\n".format(pre_mi, i_left_pos, i_right_pos, i_size)
                fout_rslt.write(s_rcd)
        samfile.close()

####
####