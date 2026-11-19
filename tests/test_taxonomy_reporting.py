import copy
import csv
import io
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

import pandas as pd

import mag2ena
from App.utils import taxonomy_resolution as taxonomy


def node():
    return {
        'taxId': '1747', 'scientificName': 'Cutibacterium acnes', 'rank': 'species',
        'submittable': 'true', 'binomial': 'true', 'metagenome': 'false',
        'lineage': 'cellular organisms;Bacteria;', 'otherNames': [],
    }


class TaxonomyReportingTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def read_tsv(self, path):
        with Path(path).open(newline='') as handle:
            return list(csv.DictReader(handle, delimiter='\t'))

    def test_review_includes_populated_but_unresolved_ids_and_preserves_input(self):
        reports = [
            {'sample_name': 'accepted', 'resolved': 'yes', 'original_tax_id': '1747'},
            {'sample_name': 'blocked', 'resolved': 'no', 'original_tax_id': '562',
             'original_organism': 'Original name_A', 'status': 'ASSEMBLY_NOT_CURRENT',
             'blocking_reason': 'Recorded reference is suppressed', 'reference_accession': 'GCF_003481705.1'},
            {'sample_name': 'unknown', 'resolved': 'no', 'original_tax_id': '',
             'status': 'NO_SPECIES_ASSIGNMENT', 'blocking_reason': 'GTDB stops above species'},
        ]
        original = copy.deepcopy(reports)
        path = taxonomy.write_unresolved_taxonomy_report(self.root / 'evidence.tsv', reports)
        rows = self.read_tsv(path)
        self.assertEqual([row['sample_name'] for row in rows], ['blocked', 'unknown'])
        self.assertEqual(rows[0]['original_tax_id'], '562')
        self.assertEqual(rows[0]['original_organism'], 'Original name_A')
        self.assertIn('other genomic evidence', rows[0]['suggested_action'])
        self.assertIn('species-rank environmental taxon', rows[1]['suggested_action'])
        self.assertEqual(reports, original)

    def test_all_resolved_rerun_clears_stale_review_rows(self):
        report = self.root / 'evidence.tsv'
        path = taxonomy.write_unresolved_taxonomy_report(
            report, [{'sample_name': 'previously_blocked', 'resolved': 'no'}])
        self.assertEqual(len(self.read_tsv(path)), 1)
        taxonomy.write_unresolved_taxonomy_report(report, [{'sample_name': 'accepted', 'resolved': 'yes'}])
        self.assertEqual(self.read_tsv(path), [])
        self.assertIn('suggested_action', path.read_text().splitlines()[0])

    def test_writer_rejects_input_collision_without_overwriting_it(self):
        evidence = self.root / 'evidence.tsv'
        target = taxonomy.unresolved_taxonomy_path(evidence)
        target.write_text('original input\n')
        with self.assertRaisesRegex(ValueError, 'overwrite'):
            taxonomy.write_unresolved_taxonomy_report(
                evidence, [], protected_paths=[target])
        self.assertEqual(target.read_text(), 'original input\n')

    def test_service_failure_guidance_preserves_the_retry_path(self):
        action = taxonomy.taxonomy_suggested_action({'status': 'LOOKUP_FAILED'})
        self.assertIn('same cache', action)
        self.assertIn('online', action)
        self.assertNotIn('--apply-generic-taxonomy', action)

    def test_generic_candidate_and_scope_corrections_require_review(self):
        action = taxonomy.taxonomy_suggested_action({
            'status': 'GENERIC_REFERENCE_MAPPING',
            'blocking_reason': 'Generic replacement requires --apply-generic-taxonomy',
        })
        self.assertIn('Review', action)
        self.assertIn('--apply-generic-taxonomy', action)
        action = taxonomy.taxonomy_suggested_action({
            'status': 'EXACT_NAME_OR_SYNONYM_MATCH',
            'blocking_reason': 'Populated taxonomy needs a correction; use --taxonomy-scope all',
        })
        self.assertIn('Review', action)
        self.assertIn('--taxonomy-scope all', action)

    def test_cli_rejects_review_input_collision_before_lookup(self):
        evidence = self.root / 'evidence.tsv'
        metadata = taxonomy.unresolved_taxonomy_path(evidence)
        metadata.write_text('sample_name\torganism\ttax_id\none\tCutibacterium acnes\t1747\n')
        original = metadata.read_bytes()
        with patch.object(mag2ena, 'resolve_taxonomy_batch') as resolver, redirect_stdout(io.StringIO()):
            status = mag2ena.main([
                '--metadata', str(metadata), '--resolve-taxonomy',
                '--taxonomy-report-out', str(evidence), '--prepare-only',
                '--prepared-metadata-out', str(self.root / 'prepared.tsv'),
            ])
        self.assertEqual(status, 2)
        resolver.assert_not_called()
        self.assertEqual(metadata.read_bytes(), original)

    def test_cli_rejects_prepared_review_collision_before_lookup(self):
        metadata = self.root / 'input.tsv'
        metadata.write_text('sample_name\torganism\ttax_id\none\tCutibacterium acnes\t1747\n')
        evidence = self.root / 'evidence.tsv'
        with patch.object(mag2ena, 'resolve_taxonomy_batch') as resolver, redirect_stdout(io.StringIO()):
            status = mag2ena.main([
                '--metadata', str(metadata), '--resolve-taxonomy',
                '--taxonomy-report-out', str(evidence), '--prepare-only',
                '--prepared-metadata-out', str(taxonomy.unresolved_taxonomy_path(evidence)),
            ])
        self.assertEqual(status, 2)
        resolver.assert_not_called()

    def test_partial_batch_writes_all_outputs_and_still_blocks_submission(self):
        metadata = self.root / 'input.tsv'
        rows = [
            {'sample_name': 'accepted', 'organism': 'Cutibacterium acnes', 'tax_id': '1747',
             'untouched': '001_A'},
            {'sample_name': 'blocked', 'organism': 'Agathobaculum sp003481705', 'tax_id': '',
             'untouched': '002_B'},
        ]
        pd.DataFrame(rows).to_csv(metadata, sep='\t', index=False)
        evidence = self.root / 'evidence.tsv'
        prepared = self.root / 'prepared.tsv'
        with patch.object(taxonomy.TaxonomyCache, 'get', return_value=node()), \
             patch.object(mag2ena, 'resolve_webin_credentials') as credentials, \
             patch.object(mag2ena, 'build_and_submit') as submit, \
             redirect_stdout(io.StringIO()) as output:
            status = mag2ena.main([
                '--metadata', str(metadata), '--resolve-taxonomy',
                '--taxonomy-cache', str(self.root / 'cache.sqlite'),
                '--taxonomy-report-out', str(evidence),
                '--prepared-metadata-out', str(prepared),
            ])
        review_path = taxonomy.unresolved_taxonomy_path(evidence)
        self.assertEqual(status, 1)
        credentials.assert_not_called()
        submit.assert_not_called()
        self.assertEqual(self.read_tsv(prepared), rows)
        self.assertEqual(len(self.read_tsv(evidence)), 2)
        review = self.read_tsv(review_path)
        self.assertEqual([row['sample_name'] for row in review], ['blocked'])
        self.assertIn('original full GTDB-Tk summary', review[0]['suggested_action'])
        for path in (evidence, prepared, review_path):
            self.assertIn(str(path.resolve()), output.getvalue())
        self.assertIn('suggested_action', output.getvalue())

    def test_successful_preparation_writes_header_only_review(self):
        metadata = self.root / 'input.tsv'
        metadata.write_text('sample_name\torganism\ttax_id\none\tCutibacterium acnes\t1747\n')
        evidence = self.root / 'evidence.tsv'
        with patch.object(taxonomy.TaxonomyCache, 'get', return_value=node()), \
             redirect_stdout(io.StringIO()):
            status = mag2ena.main([
                '--metadata', str(metadata), '--resolve-taxonomy', '--prepare-only',
                '--taxonomy-report-out', str(evidence),
                '--taxonomy-cache', str(self.root / 'cache.sqlite'),
                '--prepared-metadata-out', str(self.root / 'prepared.tsv'),
            ])
        self.assertEqual(status, 0)
        self.assertEqual(self.read_tsv(taxonomy.unresolved_taxonomy_path(evidence)), [])


if __name__ == '__main__':
    unittest.main()
