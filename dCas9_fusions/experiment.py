import numpy as np
import pandas as pd

from . import layout
import hits.utilities
import knock_knock.pacbio_experiment

memoized_property = hits.utilities.memoized_property

class dCas9FusionExperiment(knock_knock.pacbio_experiment.PacbioExperiment):
    def __init__(self, base_dir, group, name, **kwargs):
        super().__init__(base_dir, group, name, **kwargs)

        label_offsets = {
            'SV40 NLS 1': 1,
            'SV40 NLS 2': 2,
            'reverse_primer': 1,
        }
        self.diagram_kwargs.update(label_offsets=label_offsets)

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
