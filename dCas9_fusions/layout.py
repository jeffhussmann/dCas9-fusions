import knock_knock.layout

from hits import interval, utilities, sam

memoized_property = utilities.memoized_property

class Layout(knock_knock.layout.Layout):
    category_order = [
        ('contains dCas9',
            ('clean domains',
             'messy domains',
            ),
        ),
        ('no dCas9',
            ('n/a',
            ),
        ),
        ('malformed layout',
            ('extra copy of primer',
             'missing a primer',
             'too short',
             'primer far from read edge',
             'primers not in same orientation',
             'no alignments detected',
            ),
        ),
    ]

    @memoized_property
    def full_length_dCas9_alignment(self):
        expected_start = self.target_info.features[('XTEN16-2xNLS-dCas9-XTEN', 'XTEN16')].start
        expected_end = self.target_info.features[('XTEN16-2xNLS-dCas9-XTEN', 'XTEN80')].end

        full_length_al = None

        for al in self.donor_alignments:
            covered = interval.get_covered_on_ref(al)
            if expected_start in covered and expected_end in covered:
                full_length_al = al
                break

        return full_length_al

    @memoized_property
    def domain_alignments(self):
        ti = self.target_info

        non_domain_als = [self.primer_alignments[5], self.primer_alignments[3], self.full_length_dCas9_alignment]

        uncovered = self.whole_read - interval.get_disjoint_covered(non_domain_als)

        gaps = uncovered.intervals
        if len(gaps) != 2:
            return None

        all_gap_covers = []

        for gap in gaps:
            gap_covers = []
            for al in self.extra_alignments:
                if sam.get_strand(al) != self.strand:
                    continue
                    
                uncovered = (gap - interval.get_covered(al)).total_length
                if uncovered < 10:
                    edit_distance = sam.total_edit_distance(al, ref_seq=ti.reference_sequences[al.reference_name])
                    gap_covers.append((edit_distance, uncovered, al))
            
            gap_covers = [al for edit_distance, uncoverd, al in sorted(gap_covers, key=lambda t: t[:2])]
            all_gap_covers.append(gap_covers)

        if self.strand == '-':
            all_gap_covers = all_gap_covers[::-1]

        if all(len(gap_covers) > 0 for gap_covers in all_gap_covers):
            return [gap_covers[0] for gap_covers in all_gap_covers]
        else:
            return None

    def categorize(self):
        self.details = 'n/a'
        self.relevant_alignments = self.original_alignments
        
        if self.seq is None or len(self.seq) <= self.target_info.combined_primer_length + 10:
            self.category = 'malformed layout'
            self.subcategory = 'too short'

        elif all(al.is_unmapped for al in self.alignments):
            self.category = 'malformed layout'
            self.subcategory = 'no alignments detected'

        elif self.extra_copy_of_primer:
            self.category = 'malformed layout'
            self.subcategory = 'extra copy of primer'

        elif self.missing_a_primer:
            self.category = 'malformed layout'
            self.subcategory = 'missing a primer'

        elif self.primer_strands[5] != self.primer_strands[3]:
            self.category = 'malformed layout'
            self.subcategory = 'primers not in same orientation'
        
        elif not self.primer_alignments_reach_edges:
            self.category = 'malformed layout'
            self.subcategory = 'primer far from read edge'

        elif self.full_length_dCas9_alignment:
            self.category = 'contains dCas9'
            if self.domain_alignments is not None:
                self.subcategory = 'clean domains'
                first_al, second_al = self.domain_alignments
                self.details = f'{first_al.reference_name},{second_al.reference_name}'

                self.relevant_alignments = self.domain_alignments + [self.primer_alignments[5],
                                                                     self.primer_alignments[3],
                                                                     self.full_length_dCas9_alignment,
                                                                    ] 

            else:
                self.subcategory = 'messy domains'
                self.details = 'n/a'

        else:
            self.category = 'no dCas9'
            self.subcategory = 'n/a'

        if self.strand == '-':
            self.relevant_alignments = [sam.flip_alignment(al) for al in self.relevant_alignments]

        self.categorized = True

        return self.category, self.subcategory, self.details

    def plot(self, relevant=True, **manual_diagram_kwargs):
        if not self.categorized:
            self.categorize()

        ti = self.target_info

        flip_target = ti.sequencing_direction == '-'

        label_offsets = {
            'SV40 NLS 1': 1,
            'SV40 NLS 2': 2,
            'reverse_primer': 1,
        }

        diagram_kwargs = dict(
            draw_sequence=False,
            flip_target=flip_target,
            split_at_indels=False,
            features_to_show=ti.features_to_show,
            center_on_primers=True,
            label_offsets=label_offsets,
        )

        for k, v in diagram_kwargs.items():
            manual_diagram_kwargs.setdefault(k, v)

        if relevant:
            als_to_plot = self.relevant_alignments
        else:
            als_to_plot = self.uncategorized_relevant_alignments

        diagram = knock_knock.visualize.ReadDiagram(als_to_plot,
                                                    ti,
                                                    **manual_diagram_kwargs,
                                                   )

        return diagram