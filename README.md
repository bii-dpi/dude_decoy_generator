# Generate DUD-E style decoys in (only) Python

The Directory of Useful Decoys (Enhanced) provides [a useful webtool](http://dude.docking.org/generate) to generate DUD-E style decoys for custom ligand sets.

However, this tool can be quite slow. DUD-E's webtool claims it can create an average of 4.5 decoys per minute, after your job finishes waiting in the queue. Initial trials with our approach exhibit performance an order of magnitude faster.

**This repository thus offers a self-contained, purely Pythonic approach to generating decoys in the DUD-E style.**

Instructions:
1. Install Python 3.7, and the corresponding RDKit, tqdm, NumPy, Pandas, and SharedArray libraries.
2. Run 00_\* to 06_\* in sequence.

*This repository is still in alpha: feel free to raise any issues and questions as issues.*

Credit: [Dimorphite-DL](https://durrantlab.pitt.edu/dimorphite-dl/) is used to generate ligands' protonated states, an essential step of DUD-E style decoy creation.

TODO:
1. Add a setup.py to do 00-04, and rename 05-07.
2. Clean 06.
3. Enhance readme.

