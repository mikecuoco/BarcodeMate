# BarcodeMate
Tool for pre-processing alignments with barcode


## Dependency

+ Python 2.7.10 or later

+ pysam 0.12 or later

## Usage

1. **Trim clipped reads**

	`python ./x_toolbox.py -T -b input.sorted.bam -r ${REF_FILE} -o output_trimmed.bam -p ./tmp -n 11`

2. **Convert reference indexed alignment to barcode (BX) indexed alignments**
	
	`python ./x_toolbox.py -C -b input.sorted.bam -o output_barcode_indexed.bam -p ./tmp -n 11`

3. **Convert reference indexed alignment to molecule barcode (MI) indexed alignments**
	
	`python ./x_toolbox.py -X -b input.sorted.bam -o output_barcode_indexed.bam -p ./tmp -n 11`

4. **Convert barcode indexed alignments to reference based alignment**
	
	`python ./x_toolbox.py -K -b input.sorted.bam -i header_template_bam.bam -o output_barcode_indexed.bam`

5. **Haplotype block statistic**

	`python ./x_toolbox.py -B -b input.sorted.bam  -o block_info.txt -n 11 -p ./tmp/`

6. **Molecule length statistic**

	`python ./x_toolbox.py -S -b barcode_indexed_sorted.bam  -o molecule_info.txt -n 11 -p ./tmp/`

7. **Generate barocde matrix**

	`python ./x_toolbox.py -G -b input.sorted.bam -i chrom_size.txt -o out_matrix.bed -n 11 -p ./tmp/ -w 10000 -k 250 -m 1`

	Matrix visualization:
	The exported matrix can be viewed through HiGlass: [NA12878 10X hg19](http://higlass.io/app/?config=ANPVPgtzT4S8fJxzH3SQPA)

- [ ] Generate matrix by haplotype information, for each block generate two matrices
