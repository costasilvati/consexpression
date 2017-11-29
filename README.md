# consexpression
Tool for RNA-Seq analysis.
# Install
Make sure that you the standard GNU build environment installed, as well as Python, Bowtie2 and TopHat. Users of Ubuntu Linux simply type:

sudo apt-get install python

sudo apt-get install bowtie2

sudo apt-get install tophat

sudo apt-get install r-base

Make download of consexpression. Unzip the file , change to the unzipped directory.

# Usage
Edit configuration file whit expermient data.

Configuration file path is: dao/CONFIG_tool.txt. A genome and GTF/GFF annotation file is necessary.

In consexpression directory type: python experiment.py dao/CONFIG_tool.txt
