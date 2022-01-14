# Generate DUD-E style decoys in (only) Python

The Directory of Useful Decoys (Enhanced) provides [a useful webtool](http://dude.docking.org/generate) to generate DUD-E style decoys for custom ligand sets.

However, this tool can be quite slow. DUD-E's webtool claims it can create an average of 4.5 decoys per minute, after your job finishes waiting in the queue. Initial trials with our approach exhibit performance an order of magnitude faster, depending on the number of CPUs available.

**This repository thus offers a self-contained, purely Pythonic approach to generating decoys in the DUD-E style.**

## Instructions

### A. Set up the data regarding compounds to be used as potential decoys

0. `00_download_ZINC_smiles.wget`: downloads SMILES for all ZINC15 compounds in stock to `data/raw_files/`.
1. `01_combine_smis.py`: combines the SMILES text files in `data/raw_files` into one file of ~15 million compounds, `all_zinc.smi`.
2. `02_get_fingerprints.py`: generates ECFP4 fingerprints for all compounds, and stores them in `data/fingerprints/ecfp4`. Support for other fingerprinting methods planned. Fingerprints are used in selecting the compounds most topologically-dissimilar to their corresponding actives.
3. `03_get_decoy_properties.py`: generates six physicochemical properties (molecular weight, logP, number of rotatable bonds, number of hydrogen bond acceptors and donors, and net charge) for all compounds, and stores them in `data/all_zinc_w_prop.csv". Decoys are selected to property-match their corresponding actives.
4. `04_sort_ligands.py`: stores arrays of these compunds sorted by these six properties. Pre-sorted arrays are stored to make property-matching during decoy selection faster through enabling binary search.

### B. Supply the SMILES of ligands to generate decoys for
A test file, `data/test_actives`, is available for trialling the following scripts' functionality.

5. `05_get_protonated.py`: generates SMILES corresponding to the protonated states of the supplied actives in pH 6-8. Decoys are found for both the original active and these forms before combination.
6. `06_get_candidates.py`: finds 3000-9000 property-matched decoy candidates for these actives and their protonated states.
7. `07_filter_candidates.py`: combines all decoy candidates across the ligands for which decoys are being found, and keeps only the 25% most topologically-dissimilar to any one ligand (according to the Tanimoto coefficient).
8. `08_write_decoys.py`: assigns candidates to actives from the filtered candidates in a manner that minimizes decoy duplication between actives.

To generate decoys for `data/test_actives`, for example:
python 05_get_protonated.py -i /data/test_actives
python 06_get_candidates.py test -i data/test_actives_protonated [num_cores]
python 07_filter_candidates.py $job_name
python 08_write_decoys.py $job_name [o]

## Requirements

Python 3.7

### pip recommended

1. progressbar2==3.55.0
2. numpy==1.21.5
3. pandas==1.3.5
4. SharedArray==3.2.1

### conda recommended

1. rdkit v2020.09.1.0. If only pip is preferred/available, the most recent compatible version of rdkit-pypi will also work.

Credit: [Dimorphite-DL](https://durrantlab.pitt.edu/dimorphite-dl/) is used to generate ligands' protonated states, an essential step of DUD-E style decoy creation.

Please raise bug reports and enhancement ideas in the Issues section.
