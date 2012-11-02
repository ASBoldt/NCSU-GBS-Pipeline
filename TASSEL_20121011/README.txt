================================================================================
NCSU Fork of TASSEL for Genotyping by Sequencing

Software version: alpha 0.0

Date or copyright date - original TASSEL code from standalone version 3, 11 Oct 2012

legal information - TASSEL is released under Lesser GPL v2 (see http://sourceforge.net/projects/tassel/)

Contact information:
Ross Whetten, ross_whetten@ncsu.edu

================================================================================
This program / pipeline was built as a tool for streamlining the processing of genotyping by sequencing data. It is a fork of the TASSEL software developed and made available by the USDA-ARS Lab with Cornell's Institute for Genomic Diversity (http://www.maizegenetics.net/).

This software takes two FASTQ files representing forward and reverse sequencing 
of clusters in the same lane of an Illumina flow cell and:
1) concatenates the forward and reverse sequences, preserving the forward sequence
   identifier
2) feeds the outputs from step 1 to a modified version FastqToTagCountPlugin that
   extracts 64-nt tags from both the forward and reverse sequences
3) Merges tag counts across lanes and creates TagsByTaxa files following the same
   strategies used by the standard GBS pipeline for single-end reads
4) Evaluates the similarity of read tag sequences to each other, identifies
   potential allelic variants, and tests for independence in a specified subset
   of the individuals represented in the TagsByTaxa table.
5) Exports a table of genotypes, encoded as codominant loci with heterozygous
   and homozygous individuals identified based on the presence of alleles 
   identified in step 4.

================================================================================
Hardware and software requirements

OS: 64-bit Linux (Ubuntu is the development environment)
RAM: 8 Gb minimum for one Hiseq lane of data; 16 Gb is better for multiple lanes
PROCESSOR: single core  
STORAGE: A single Hiseq lane of data can be almost 100 Gb; a terabyte of available 
space is not excessive
================================================================================
Installation instructions
--------------------------------------------------------------------------------
getting started tips
--------------------------------------------------------------------------------
documentation
================================================================================
Important known problems
================================================================================
Version history






