import copy
import csv
import io
import json
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

import pandas as pd

import mag2ena
from App.utils import taxonomy_resolution as taxonomy


def node(taxid='1747', organism='Cutibacterium acnes', **changes):
    result = {
        'taxId': taxid, 'scientificName': organism, 'rank': 'species',
        'submittable': 'true', 'binomial': 'true', 'metagenome': 'false',
        'lineage': 'cellular organisms;Bacteria;', 'otherNames': [],
    }
    result.update(changes)
    return result


class TaxonomyRoutingRegressionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def frame(self, organism='Cutibacterium acnes', taxid='1747'):
        return pd.DataFrame([{
            'sample_name': 'mag_1', 'organism': organism,
            'tax_id': taxid, 'untouched': '001_A',
        }])

    def evidence(self, classification, **extra):
        return {
            'mag_1': {
                'row': {'user_genome': 'mag_1', 'classification': classification, **extra},
                'sources': [{'path': '/original/summary.tsv', 'sha256': 'fixture'}],
            },
        }

    def run_batch(self, classification='d__Bacteria;g__Cutibacterium;s__Cutibacterium acnes',
                  frame=None, evidence=None, payloads=None, **options):
        frame = self.frame() if frame is None else frame
        evidence = self.evidence(classification) if evidence is None else evidence
        payloads = {taxonomy.ENA + 'tax-id/1747': node()} if payloads is None else payloads

        def get(cache, url, kind='json'):
            result = payloads.get(url)
            if isinstance(result, Exception):
                raise result
            return copy.deepcopy(result)

        with patch.object(taxonomy.TaxonomyCache, 'get', new=get), \
             patch.object(taxonomy, 'index_metadata') as index:
            result = taxonomy.resolve_taxonomy_batch(
                frame, evidence=evidence, release='214.1',
                cache_path=self.root / 'cache.sqlite', **options,
            )
        return (*result, index)

    def test_existing_named_pair_works_with_two_column_summary(self):
        result, report, summary, index = self.run_batch()
        self.assertEqual(summary['resolved'], 1)
        self.assertEqual(summary['applied'], 0)
        self.assertEqual(report[0]['mapping_scope'], 'NAME_ID_CONSISTENCY_ONLY')
        self.assertEqual(report[0]['reference_accession'], '')
        self.assertEqual(result.loc[0, 'tax_id'], '1747')
        index.assert_not_called()

    def test_conventional_named_species_can_fill_missing_id(self):
        result, report, summary, _ = self.run_batch(
            frame=self.frame(taxid=''), payloads={
                taxonomy.ENA + 'any-name/Cutibacterium%20acnes': [node()],
            })
        self.assertEqual(summary['resolved'], 1)
        self.assertEqual(result.loc[0, 'tax_id'], '1747')
        self.assertEqual(report[0]['mapping_scope'], 'NAME_ID_CONSISTENCY_ONLY')

    def test_named_species_disagreement_cannot_be_silently_accepted(self):
        _, report, summary, _ = self.run_batch(
            classification='d__Bacteria;g__Othergenus;s__Othergenus differentii')
        self.assertEqual(len(summary['unresolved']), 1)
        self.assertEqual(report[0]['status'], 'METADATA_GTDB_SPECIES_CONFLICT')
        self.assertEqual(report[0]['candidate_verified'], 'no')

    def test_existing_generic_species_pair_survives_no_gtdb_species(self):
        organism = 'uncultured Mycoplasmoidaceae bacterium'
        result, report, summary, index = self.run_batch(
            classification='d__Bacteria;f__Mycoplasmoidaceae;g__;s__',
            frame=self.frame(organism, '3383102'), payloads={
                taxonomy.ENA + 'tax-id/3383102': node('3383102', organism, binomial='false'),
            })
        self.assertEqual(summary['resolved'], 1)
        self.assertEqual(report[0]['applied'], 'no')
        self.assertEqual(result.loc[0, 'tax_id'], '3383102')
        self.assertEqual(report[0]['proposed_rank'], 'species')
        index.assert_not_called()

    def test_above_species_placement_reference_does_not_invalidate_existing_generic_pair(self):
        organism = 'uncultured Mycoplasmoidaceae bacterium'
        evidence = self.evidence(
            'd__Bacteria;f__Mycoplasmoidaceae;g__;s__',
            closest_placement_reference='GCF_000005845.2',
        )
        _, _, summary, index = self.run_batch(
            evidence=evidence, frame=self.frame(organism, '3383102'), payloads={
                taxonomy.ENA + 'tax-id/3383102': node('3383102', organism, binomial='false'),
            })
        self.assertEqual(summary['resolved'], 1)
        index.assert_not_called()

    def test_new_generic_species_mapping_can_be_disabled(self):
        organism = 'uncultured Mycoplasmoidaceae bacterium'
        payloads = {
            taxonomy.ENA + 'any-name/uncultured%20Mycoplasmoidaceae%20bacterium': [
                node('3383102', organism, binomial='false')],
        }
        _, report, summary, _ = self.run_batch(
            classification='d__Bacteria;f__Mycoplasmoidaceae;g__;s__',
            frame=self.frame(organism, ''), payloads=payloads, apply_generic=False)
        self.assertEqual(len(summary['unresolved']), 1)
        self.assertIn('--apply-generic-taxonomy', report[0]['blocking_reason'])
        result, _, summary, _ = self.run_batch(
            classification='d__Bacteria;f__Mycoplasmoidaceae;g__;s__',
            frame=self.frame(organism, ''), payloads=payloads)
        self.assertEqual(summary['resolved'], 1)
        self.assertEqual(result.loc[0, 'tax_id'], '3383102')

    def test_generic_option_does_not_accept_genus_or_family_rank(self):
        organism = 'uncultured Mycoplasmoidaceae bacterium'
        for rank in ('genus', 'family'):
            with self.subTest(rank=rank):
                _, _, summary, _ = self.run_batch(
                    classification='d__Bacteria;f__Mycoplasmoidaceae;g__;s__',
                    frame=self.frame(organism, '3383102'), apply_generic=True, payloads={
                        taxonomy.ENA + 'tax-id/3383102':
                            node('3383102', organism, rank=rank, binomial='false'),
                    })
                self.assertEqual(len(summary['unresolved']), 1)

    def test_above_species_classification_cannot_create_named_descendant(self):
        _, report, summary, _ = self.run_batch(
            classification='d__Bacteria;g__Cutibacterium;s__',
            frame=self.frame(taxid=''), payloads={
                taxonomy.ENA + 'any-name/Cutibacterium%20acnes': [node()],
            })
        self.assertEqual(len(summary['unresolved']), 1)
        self.assertEqual(report[0]['status'], 'NO_SPECIES_ASSIGNMENT')

    def test_numbered_placeholder_reports_missing_columns_and_skips_index(self):
        _, report, summary, index = self.run_batch(
            classification='d__Bacteria;g__Agathobaculum;s__Agathobaculum sp003481705',
            frame=self.frame('Agathobaculum sp003481705', ''), payloads={})
        self.assertEqual(len(summary['unresolved']), 1)
        self.assertEqual(report[0]['status'], 'GTDB_REFERENCE_COLUMNS_MISSING')
        self.assertIn('original GTDB-Tk summary', report[0]['blocking_reason'])
        index.assert_not_called()

    def test_legacy_spaces_do_not_hide_raw_gtdb_species_suffix(self):
        original = self.frame('Ruminococcus G avistercoris', '')
        result, report, _, _ = self.run_batch(
            classification='d__Bacteria;g__Ruminococcus_G;s__Ruminococcus_G avistercoris',
            frame=original, payloads={})
        self.assertEqual(report[0]['status'], 'GTDB_REFERENCE_COLUMNS_MISSING')
        self.assertIn('Ruminococcus_G', report[0]['gtdb_classification'])
        pd.testing.assert_frame_equal(result, original)

    def test_legacy_spaces_do_not_hide_raw_gtdb_genus_suffix(self):
        _, report, summary, _ = self.run_batch(
            classification='d__Bacteria;g__Olsenella_F;s__',
            frame=self.frame('uncultured Olsenella F bacterium', ''), payloads={})
        self.assertEqual(len(summary['unresolved']), 1)
        self.assertEqual(report[0]['status'], 'NO_SPECIES_ASSIGNMENT')

    def test_recorded_reference_threshold_failure_has_no_name_fallback(self):
        raw = {
            'fastani_reference': 'GCF_000005845.2', 'fastani_reference_radius': '95',
            'fastani_taxonomy': 'd__Bacteria;g__Cutibacterium;s__Cutibacterium acnes',
            'fastani_ani': '90', 'fastani_af': '0.8',
        }
        with patch.object(taxonomy, 'resolve_named') as named:
            _, report, summary, _ = self.run_batch(
                evidence=self.evidence(raw['fastani_taxonomy'], **raw), payloads={})
        self.assertEqual(len(summary['unresolved']), 1)
        self.assertEqual(report[0]['status'], 'REFERENCE_BELOW_SPECIES_THRESHOLD')
        named.assert_not_called()

    def test_blank_reference_values_are_distinct_from_missing_columns(self):
        classification = 'd__Bacteria;g__Agathobaculum;s__Agathobaculum sp003481705'
        _, report, _, index = self.run_batch(
            evidence=self.evidence(classification, fastani_reference='N/A'),
            frame=self.frame('Agathobaculum sp003481705', ''), payloads={})
        self.assertEqual(report[0]['status'], 'NO_RECORDED_SPECIES_REFERENCE')
        index.assert_not_called()

    def test_named_lookup_service_failure_is_not_reported_as_a_match(self):
        _, report, summary, _ = self.run_batch(payloads={
            taxonomy.ENA + 'tax-id/1747': taxonomy.Unresolved('LOOKUP_FAILED', 'HTTP 503'),
        })
        self.assertEqual(report[0]['status'], 'LOOKUP_FAILED')
        self.assertEqual(len(summary['unresolved']), 1)

    def test_partial_success_still_stops_cli_before_credentials(self):
        metadata = self.root / 'input.tsv'
        summary = self.root / 'classification_only.tsv'
        report_path = self.root / 'report.tsv'
        prepared = self.root / 'prepared.tsv'
        frame = pd.concat([
            self.frame(),
            self.frame('Agathobaculum sp003481705', '').assign(sample_name='mag_2'),
        ], ignore_index=True)
        frame.to_csv(metadata, sep='\t', index=False)
        with summary.open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['user_genome', 'classification'])
            writer.writerow(['mag_1', 'd__Bacteria;g__Cutibacterium;s__Cutibacterium acnes'])
            writer.writerow(['mag_2', 'd__Bacteria;g__Agathobaculum;s__Agathobaculum sp003481705'])
        with patch.object(taxonomy.TaxonomyCache, 'get', return_value=node()), \
             patch.object(taxonomy, 'index_metadata') as index, \
             patch.object(mag2ena, 'resolve_webin_credentials') as credentials, \
             patch.object(mag2ena, 'build_and_submit') as submit, \
             redirect_stdout(io.StringIO()) as output:
            status = mag2ena.main([
                '--metadata', str(metadata), '--gtdbtk', str(summary),
                '--gtdb-release', '214.1', '--resolve-taxonomy',
                '--taxonomy-report-out', str(report_path),
                '--prepared-metadata-out', str(prepared),
                '--taxonomy-cache', str(self.root / 'cli_cache.sqlite'),
            ])
        self.assertEqual(status, 1)
        credentials.assert_not_called()
        submit.assert_not_called()
        index.assert_not_called()
        self.assertIn('GTDB_REFERENCE_COLUMNS_MISSING', output.getvalue())
        with report_path.open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual([row['resolved'] for row in rows], ['yes', 'no'])
        self.assertEqual(json.loads(rows[1]['gtdb_record_json'])['classification'],
                         'd__Bacteria;g__Agathobaculum;s__Agathobaculum sp003481705')


if __name__ == '__main__':
    unittest.main()
