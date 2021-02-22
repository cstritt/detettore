*detettore* – a program to detect transposable element polymorphisms
================
February 2021

![](detettore_ad.png)

*Detettore* is a program to detect and characterize transposable element
(TE) polymorphisms using reference-aligned paired-end reads.

The program searches for:

  - **TE insertion polymorphisms (TIPs)**, i.e. TEs absent in the
    reference genome but present in a sequenced individual
  - **TE absence polymorphisms (TAPs)**, i.e. TEs present in the
    reference but absent in a sequenced individual


**Version 2, work in progress**.
Main changes:

  - all output in vcf format
  - easier parameters for filtering
  - tools for downstream analysis
  - compatibility with minimap2
