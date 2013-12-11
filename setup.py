#!/usr/bin/env python

from distutils.core import setup

setup(
    name='rover',
    version='1.0.1',
    author='Bernie Pope',
    author_email='bjpope@unimelb.edu.au',
    packages=['rover'],
    scripts=['rover/get_primer_blocks.py'],
    entry_points={
        'console_scripts': ['rover = rover.rover:main']
    },
    url='https://github.com/bjpop/rover',
    license='LICENSE.txt',
    description=(
        'ROVER-PCR Variant Caller: read-pair overlap '
        'considerate variant-calling software for PCR-based '
        'massively parallel sequencing datasets.'),
    long_description=(
        'The software accepts a tab-delimited file (TSV) '
        'listing the coordinates of the target-specific primers used '
        'for targeted enrichment based on a specified genome-build. It '
        'also accepts aligned sequence files resulting from mapping to '
        'the same genome-build. ROVER-PCR Variant Caller identifies the '
        'amplicon a given read-pair represents and removes the primer '
        'sequences by using the mapping co-ordinates and primer '
        'co-ordinates. It considers completely overlapping read-pairs '
        'with respect to primer-intervening sequence. Only when a variant '
        'is observed in both reads of a read-pair does the signal '
        'contribute to a tally of read-pairs containing or not containing '
        'the variant. A user-defined threshold informs the minimum number '
        'of and percentage of read-pairs a variant must be observed in for '
        'a "call" to be made. ROVER-PCR Variant Caller also reports the '
        'depth of coverage across amplicons to facilitate the identification '
        'of any regions that may require further screening.'),
    install_requires=[
        "pysam >= 0.7.5"
    ],
)
