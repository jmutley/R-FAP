R-FAP
=====

This repository contains the Perl source code for the R-FAP funtion prediction tool.

The tool is entirely contained within the Perl script 'R-FAP_annotation_for_multiple_genomes.pl' which, obviously, requires Perl to run, but also has several other dependencies. First, there are a set of file dependencies; in addition to the FASTA-formatted pan-genome (or perhaps more accuately, pan-proteome) of the taxon under study as well as the '.faa' files from each strain in question (grouped into a single directory), the annotations for the database are input as a tab-delimited file consisting of three columns: the database sequence name, the annotation tool used to determine function a priori, and the predicted funtion itself. 

Subsequently, there are the module depenencies. The most important is the Simrank k-mer matching module from DeSantis et al., which can be found at http://search.cpan.org/perldoc?String::Simrank

The Simrank module has a set of dependencies of its own, which will be spelled out here and included in this repository at a future date.

Once these modules are installed, the user may run the annotation tool in any directory (since all variables require complete filepaths to operate) by simply entering "perl R-FAP_annotation_for_multiple_genomes.pl" and following the on-screen directions.
