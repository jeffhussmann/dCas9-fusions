Installation
------------

First, clone up-to-date versions of this repo and of two dependencies from github:

```
git clone git@github.com:jeffhussmann/hits.git
git clone git@github.com:jeffhussmann/knock-knock.git
git clone git@github.com:jeffhussmann/dCas9-fusions.git
```

(Since dCas9-fusions is a private repo, you will probably need to set up an SSH key to authenticate with your github account - https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)

Then, in the order `hits`, then `knock-knock`, then `dCas9-fusions`, go into each cloned directory and install the corresponding python using pip by running this command:

```
pip install -e ./
```

You will also need to install some non-python dependencies. In a conda environment, run:
```
conda install blast
conda install samtools
```

Running
-------
Installation will have put an executable script `count_dCas9_fusion_domains` in your PATH.
Identify a directory `BASE_DIR` that you want to hold data and results.

First run

```count_dCas9_fusion_domains setup $BASE_DIR```

substituting your actual directory for `$BASE_DIR`. This will install references sequences for the domain libraries into `$BASE_DIR/target`:

```
$ tree $BASE_DIR/target
$BASE_DIR/target
├── pTwist-SFFV_Library1
│   ├── manifest.yaml
│   ├── ptwist-sffv.gb
│   ├── refs.fasta
│   ├── refs.fasta.fai
│   ├── refs.gff
│   ├── TwistLibrary1.fasta
│   └── xten16-2xnls-dcas9-xten80.gb
└── pTwist-SFFV_Library2
    ├── manifest.yaml
    ├── ptwist-sffv.gb
    ├── refs.fasta
    ├── refs.fasta.fai
    ├── refs.gff
    ├── TwistLibrary2.fasta
    └── xten16-2xnls-dcas9-xten80.gb
```

Then setup input data in `$BASE_DIR/data`, organized into subdirectories for each batch of samples. For example, here are the contents for a single batch `220610_pL1_pL2`:

```
$ tree $BASE_DIR/data
$BASE_DIR/data
└── 220610_pL1_pL2
    ├── pL1_trimmed.fastq.gz
    ├── pL2_trimmed.fastq.gz
    └── sample_sheet.yaml
 ```
 
The fastq.gz files should be the CCS reads from Pacbio output (i.e. corresponding to files originally named `demux/A02.demux.ccs_kinetics.pL1--pL1.fastq.gz` and `demux/A02.demux.ccs_kinetics.pL2--pL2.fastq.gz` from this batch) that have had barcodes trimmed from the beginning and ends, so that each read is expected to start and end with the relevant parts of the amplicon primers (`CATTGTCGATCCTACCATCCACTCGAC` and `CCATCCGCACGCATCTGGAATAAG`). Example code for how to do this trimming can be found in [notebooks/trimming_example.ipynb](notebooks/trimming_example.ipynb)

Each batch subdirectory should contain `sample_sheet.yaml` that associates sample names with fastq files and indicates which library was used for each sample:

```
$ cat sample_sheet.yaml
pL1:
    CCS_fastq_fn: pL1_trimmed.fastq.gz
    target_info: pTwist-SFFV_Library1

pL2:
    CCS_fastq_fn: pL2_trimmed.fastq.gz
    target_info: pTwist-SFFV_Library2
```

Then run `count_dCas9_fusion_domains parallel` to count occurences of domain pairs in each sample in parallel, with syntax:

```
$ count_dCas9_fusion_domains parallel -h
usage: count_domains parallel [-h] [--batch BATCH] [--stages STAGES]
                              base_dir max_procs

positional arguments:
  base_dir         the base directory to store input data, reference
                   annotations, and analysis output for a project
  max_procs        maximum number of samples to process at once

optional arguments:
  -h, --help       show this help message and exit
  --batch BATCH    if specified, the single batch name to process; if not
                   specified, all groups will be processed
  --stages STAGES
 ```
 
 After processing is finished, you can access a DataFrame with counts for each domain pair for each experiment as follows:
 
 ```
 $ ipython

In [1]: import dCas9_fusions.experiment

In [2]: base_dir = $BASE_DIR

In [3]: exps = dCas9_fusions.experiment.get_all_experiments(base_dir)

In [4]: exps
Out[4]:
{('220610_pL1_pL2',
  'pL1'): dCas9FusionExperiment: batch=220610_pL1_pL2, sample_name=pL1, base_dir=$BASE_DIR,
 ('220610_pL1_pL2',
  'pL2'): dCas9FusionExperiment: batch=220610_pL1_pL2, sample_name=pL2, base_dir=$BASE_DIR}

In [5]: exps['220610_pL1_pL2', 'pL1'].domain_counts
Out[5]:
         count  length  count_with_pseudocount  log10_count  fraction  fraction_with_floor  log10_fraction
N1  C1      13    4968                    13.1     1.117271  0.000890             0.000890       -3.050796
    C2      31    5946                    31.1     1.492760  0.002121             0.002121       -2.673378
    C3      24    6000                    24.1     1.382017  0.001642             0.001642       -2.784528
    C4       1    8001                     1.1     0.041393  0.000068             0.000068       -4.164739
    C5      12    5229                    12.1     1.082785  0.000821             0.000821       -3.085558
...        ...     ...                     ...          ...       ...                  ...             ...
N39 C35      1    7920                     1.1     0.041393  0.000068             0.000068       -4.164739
    C36      1    8526                     1.1     0.041393  0.000068             0.000068       -4.164739
    C37      8    7290                     8.1     0.908485  0.000547             0.000547       -3.261649
    C38      0    7299                     0.1    -1.000000  0.000000             0.000001       -6.000000
    C39      4    8274                     4.1     0.612784  0.000274             0.000274       -3.562679

[1369 rows x 7 columns]

```