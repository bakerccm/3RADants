This file provides a mapping from the original raw data files to simplified file names that work for demultiplexing with process_radtags.

This mapping is used by `snakemake` to generate soft links that point to the original data files.

**original** provides the path to the original raw data file *relative to the data folder*. (e.g. `data/Rawdata/folderA/fileA.fasta.gz` is shown as `Rawdata/folderA/fileA.fasta.gz`).

**link** provides the path for the soft link *relative to the links folder*. (e.g. out/data/plate1/file1.fasta.gz is shown as plate1/file1.fasta.gz).

The leading part of the path is supplied by the snakefile (which makes it easier to relocate the files if desired).
