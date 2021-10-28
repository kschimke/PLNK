# PLNK
PLNK runs alongside an Oxford Nanopore MinION sequencer, processing individual fast5 files using guppy for basecalling, C3POa for R2C2 consensus calling, and mappy for alignment before analyzing the library content.

Expects R2C2 reads using samples previously prepped for Illumina sequencing. 

## Dependencies

- [Python 3](https://www.python.org/downloads/)
- [NumPy](https://pypi.org/project/numpy/)
- [guppy]()
- [C3Poa](https://github.com/rvolden/C3POa)
- [mappy](https://pypi.org/project/mappy/)
