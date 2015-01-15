--------------------------------------------------------------------------------
Rover - ROVER-PCR Variant Caller: read-pair overlap considerate variant-calling
software for PCR-based massively parallel sequencing datasets
--------------------------------------------------------------------------------

Version: 2.0.0

Authors: Bernard J Pope (2,3), Tú Nguyen-Dumont (1), Fleur Hammet (1), 
         Daniel J Park (1) and Roger Li

         (1) Genetic Epidemiology Laboratory, Department of Pathology,
             The University of Melbourne.
         (2) Victorian Life Sciences Computation Initiative (VLSCI).
         (3) Department of Computing and Information Systems,
             The University of Melbourne.
         

Web:     https://github.com/bjpop/rover

License: BSD

Citation:

   Please cite Rover as follows:

   ROVER variant caller: read-pair overlap considerate variant-calling software
   applied to PCR-based massively parallel sequencing datasets.
   Bernard J Pope, Tú Nguyen-Dumont, Fleur Hammet and Daniel J Park.

   Source Code for Biology and Medicine 2014, 9:3
   DOI: 10.1186/1751-0473-9-3
   http://www.scfbm.org/content/9/1/3

Requirements: Python 2.7, and the PySam, PyVCF and biopython libraries
(http://code.google.com/p/pysam/)
(https://pypi.python.org/pypi/PyVCF)
(https://pypi.python.org/pypi/biopython)

--------------------------------------------------------------------------------
General description
--------------------------------------------------------------------------------

ROVER-PCR Variant Caller enables users to quickly and accurately identify
genetic variants from PCR-targeted, overlapping paired-end MPS datasets. The
open-source availability of the software and threshold tailorability enables
broad access for a range of PCR-MPS users.

ROVER-PCR Variant Caller is implemented in Python and runs on any Posix 
compliant system (Linux, OS X).

The software accepts a tab-delimited file (TSV) listing the coordinates of the
target-specific primers used for targeted enrichment based on a specified
genome-build. It also accepts aligned sequence files resulting from mapping
to the same genome-build. ROVER-PCR Variant Caller identifies the amplicon
a given read-pair represents and removes the primer sequences by using the
mapping co-ordinates and primer co-ordinates. It can also optionally check the 
primer sequences at the 5' ends to ensure that they are within a certain
threshold score away from the expected sequence, based on a gapped alignment.

It considers completely overlapping read-pairs with respect to primer-
intervening sequence. Only when a variant is observed in both reads of a 
read-pair does the signal contribute to a tally of read-pairs containing or not 
containing the variant. A user-defined threshold informs the minimum number of 
and proportion of read-pairs a variant must be observed in for a ‘call’ to be 
made. ROVER-PCR Variant Caller also reports the depth of coverage across 
amplicons to facilitate the identification of any regions that may require 
further screening.

--------------------------------------------------------------------------------
Command line usage:
--------------------------------------------------------------------------------

usage: rover [-h] [--version] --primers PRIMERS [--overlap OVERLAP]
             [--log FILE] --out FILE [--proportionthresh N] [--absthresh N]
             [--qualthresh N] [--primercheck SEQ] [--primerthresh N]
             [--gap_penalty N] [--id_info DBSNP] [--coverdir COVERDIR]
             [--datadir DATADIR] bams [bams ...]

Consider mapped reads to amplicon sites

positional arguments:
  bams                  BAM files containing mapped reads.

optional arguments:
  -h, --help            Show this help message and exit.
  --version             Show program's version number and exit.
  --primers PRIMERS     File name of primer coordinates in TSV format.
  --overlap OVERLAP     Minimum fraction overlap of read to block region.
                        Defaults to 0.9.
  --log FILE            Log progress in FILENAME, defaults to stdout.
  --out FILE            Name of output file containing called variants.
  --proportionthresh N  Keep variants which appear in this proportion of the
                        read pairs for a given target region, and bin
                        otherwise. Defaults to 0.05.
  --absthresh N         Only keep variants which appear in at least this many
                        read pairs. Defaults to 2.
  --qualthresh N        Minimum base quality score (phred).
  --primercheck SEQ     Rover will check the primer sequences at the 5' ends
                        by aligning them against the sequences in SEQ
  --primerthresh N      If the difference in score between the primer sequence 
                        in a read and the expected sequence is greater than N, 
                        the read will be discarded. Defaults to 5.0.
  --gap_penalty N       The deduction in score in the gapped alignment for gaps
                        gap_penalty * 2 for opening a gap and gap_penalty for
                        extending a gap. Defaults to 2.0.
  --id_info DBSNP       The dbsnp file in .vcf.gz format and accompanying .tbi
                        index file containing rs numbers for known variants. 
  --coverdir COVERDIR   Directory to write coverage files, defaults to current
                        working directory.
  --datadir DATADIR     Directory to write data files, defaults to current
                        working directory. 

Explanation of the arguments:

   -h

      Print a help message and exit.

   --version

      Print the version number of rover and exit.

   --primers PRIMERS

      Required.

      A tab-delimited list of the coordinates of the target-specific primers
      used for targeted enrichment based on a specified genome-build.
      The format is compatible with the BED file format:

      chromosome start end 

   --overlap OVERLAP

      Optional. Defaults to 0.9.

      Minimum fraction overlap of read to block region.

      0.5 means at least half of a block must be overlapped by a read
      for the read to be considered for that region. 1.0 would mean the entire
      block must be overlapped by the read.

   --log FILE

      Optional. Defaults to the standard output.

      Write a log of the program's progress to this file. 

   --out FILE

      Required.

      Name of the output file created by Rover. This file contains the
      variants called by the program. This file is in VCF format.

   --proportionthresh N

      Optional. Defaults to 0.05.

      Only keep variants which appear in this proportion of the read pairs for
      a given target region, and bin otherwise. A variant must appear in
      both reads of a pair to be counted. The proportion is calculated as
      follows:

          N = number of pairs containing this variant in both reads
          T = number of read pairs overlapping the target region

          proportion = N/T


      Note: variants must pass BOTH the proportionthresh and absthresh 
      thresholds to be kept. If they fail either test then they are binned.

      That is to say:

          if N/T  >= proportionthresh and N >= absthresh:
              keep the variant
          else:
              bin the variant

   --absthresh N

      Optional. Defaults to 2 read pairs.

      Only keep variants which appear in at least this many read pairs
      for a given target region.

      See comments above about proportionthresh.

   --qualthresh N

      Optional. Minimum phred quality score for bases appearing in SNVs and
      insertions. If this argument is set Rover will only consider SNVs and
      insertions where all DNA bases in those variants have a quality score
      greater than or equal to the argument. For example, setting this
      argument to 35 will cause Rover to discard any SNVs or insertions
      containing any bases with a score below 35. If the argument is not
      set then Rover will not consider quality scores in its decision
      to keep or discard a variant.

   --primercheck SEQ

      Optional. If provided, Rover will check the primer sequences of at the 5'
      ends of the reads, and assign a score based on a gapped alignment using 
      the pairwise2 module from Biopython. It will compare the primer sequence 
      in the read with the expected primer sequences located in SEQ, which is a 
      TSV in the following format:

      primer_name   expected_sequence

      Rover will also create .dat files containing statistics related to the
      primer checking process. These files will indicate the percentage of
      primers that are a certain score away from being an exact match. What the
      numbers actually mean will depend on the gap_penalty specified. A 
      directory called primer_data will need to be created, in which the .dat 
      files for each sample will be placed. Rover will not try to create this 
      directory. 

   --primerthresh N

      Optional. Defaults to 5.0.

      The maximum allowed difference in score for a primer from an exact match 
      before the read pair is discarded. 

   --gap_penalty N

      Optional. Defaults to 2.0. 

      The score deduction for the gapped alignment used in Rover (for checking 
      the primer sequences). A score deduction of gap_penalty * 2 is applied for
      opening a gap and a deduction of gap_penalty for extending a gap. 

   --id_info DBSNP

      Optional. Rover will find the rs numbers in DBSNP and show them in the 
      output VCF file if it can find the relevant number in DBSNP, which is a 
      vcf format file containing rs numbers for known variants in .vcf.gz format
      with an accompanying .tbi file (created by tabix).
 
   --coverdir COVERDIR

      Optional. Defaults to current working directory.

      Directory to write the coverage files. If the directory does not
      exist Rover will not try to create it. You must create the directory
      yourself.

      The format of the coverage files is a TSV with:

      chr     block_start     block_end       num_pairs

    --datadir DATADIR

      Optional. Defaults to current working directory. 

      Only relevant if primer checking is enabled. Data files will be created
      inside the directory DATADIR. "total.dat" contain statistics relating
      to the entire dataset, while "sample1.dat", "sample2.dat", etc. will
      contain statistics relating to the individual BAM files.

   bams [bams ...]

      One or more BAM files containing mapped reads.

--------------------------------------------------------------------------------
Example usage (should be all on one line)
--------------------------------------------------------------------------------

   rover --primers primer_coords.tsv --log rover_log --out variants.vcf
         --proportionthresh 0.15 --absthresh 2 --id_info dbsnp.vcf.gz 
         --coverdir coverage_files sample1.bam sample2.bam sample3.bam

   With primer checking:

   rover --primers primer_coord.tsv --log rover_log --out variants.vcf
         --proportionthresh 0.15 --absthresh 2 --id_info dbsnp.vcf.gz
         --primercheck primers.tsv --primerthresh 5.0 --gap_penalty 2.0
         --coverdir coverage_files --datadir primer_data
         sample1.bam sample2.bam sample3.bam

This assumes that the coordinates for your primer regions are in the file
primer_coords.tsv. The detected variants will appear in the output file
variants.vcf. The names of the samples will be taken from the prefix of
the bam file name, in this case "sample1" "sample2" and "sample3". 

If primer checking is enabled, the data file for the entire dataset will
be called will be output in primer_data/total.dat, and data files for
the individual BAM files will be output in primer_data/sample1.dat, 
primer_data/sample2.dat and primer_data/sample3.dat.

Coverage files containing the number of read pairs which mapped to
each region will be output in coverage_files/sample1.coverage
coverage_files/sample2.coverage and coverage_files/sample3.coverage.
A log file describing the actions taken by the program
will be stored in rover_log.

--------------------------------------------------------------------------------
