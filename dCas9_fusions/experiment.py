import logging
import sys

import numpy as np
import pandas as pd

import hits.utilities

import knock_knock.experiment
import knock_knock.pacbio_experiment

from . import layout

memoized_property = hits.utilities.memoized_property

class dCas9FusionExperiment(knock_knock.pacbio_experiment.PacbioExperiment):
    @memoized_property
    def categorizer(self):
        return layout.Layout

    @memoized_property
    def domain_counts(self):
        refs = self.target_info.reference_sequences

        N_domains = [name for name in refs if name.startswith('N')]
        C_domains = [name for name in refs if name.startswith('C')]

        self.outcome_counts.sort_index(inplace=True)
        clean_counts = self.outcome_counts.loc['contains dCas9', 'clean domains']

        domain_counts = {}

        for N in N_domains:
            for C in C_domains:
                domain_counts[N, C] = clean_counts.get(f'{N},{C}', 0)

        domain_counts = pd.Series(domain_counts)

        dCas9_length = len(refs['XTEN16-2xNLS-dCas9-XTEN'])

        lengths = {(N, C): len(refs[N]) + len(refs[C]) + dCas9_length for N, C in domain_counts.index}

        df = pd.DataFrame({
            'count': domain_counts,
            'length': lengths,
        })

        df['count_with_pseudocount'] = df['count'] + 0.1
        df['log10_count'] = np.log10(df['count_with_pseudocount'])
        df['fraction'] = df['count'] / df['count'].sum()
        df['fraction_with_floor'] = df['fraction']
        df.loc[df['count'] == 0, 'fraction_with_floor'] = 1e-6
        df['log10_fraction'] = np.log10(df['fraction_with_floor'])

        return df

def get_all_experiments(base_dir):
    exps = {}
    batches = knock_knock.experiment.get_all_batches(base_dir)

    for batch in batches:
        sample_sheet = knock_knock.experiment.load_sample_sheet(base_dir, batch)

        if sample_sheet is None:
            print(f'Error: {batch} has no sample sheet')
            continue

        for name, description in sample_sheet.items():
            if isinstance(description, str):
                continue

            exp = dCas9FusionExperiment(base_dir, batch, name, description=description)
            exps[exp.batch, exp.sample_name] = exp

    return exps

def process_experiment_stage(base_dir,
                             batch_name,
                             sample_name,
                             stage,
                            ):
    sample_sheet = knock_knock.experiment.load_sample_sheet(base_dir, batch_name)

    if sample_sheet is None:
        print(f'Error: {batch_name} not found in {base_dir}')
        sys.exit(1)
    elif sample_name not in sample_sheet:
        print(f'Error: {sample_name} not found in {batch_name} sample sheet')
        sys.exit(1)
    else:
        description = sample_sheet[sample_name]

    exp = dCas9FusionExperiment(base_dir, batch_name, sample_name, description=description)

    logging.info(f'Starting {batch_name}: {sample_name} {stage}')

    exp.process(stage)

    logging.info(f'Finished {batch_name}: {sample_name} {stage}')
