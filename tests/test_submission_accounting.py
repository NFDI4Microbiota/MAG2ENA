import io
import json
import subprocess
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import Mock, patch

import polars as pl

import mag2ena


def receipt(aliases, errors=(), success=True):
    root = mag2ena.ET.Element('RECEIPT', {'success': str(success).lower()})
    for alias in aliases:
        mag2ena.ET.SubElement(root, 'SAMPLE', {'alias': alias, 'accession': 'ERS_' + alias})
    if errors:
        messages = mag2ena.ET.SubElement(root, 'MESSAGES')
        for message in errors:
            mag2ena.ET.SubElement(messages, 'ERROR').text = message
    return mag2ena.ET.tostring(root, encoding='unicode')


class SubmissionAccountingTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def frame(self, count=2):
        return pl.DataFrame([
            {'sample_name': f'mag_{index}', 'organism': 'Escherichia coli',
             'tax_id': '562', 'genome coverage': '20', 'assembly software': 'assembler',
             'platform': 'ILLUMINA'} for index in range(count)
        ])

    def run_submission(self, frame, receipts, java_results=None, missing_fastas=()):
        fasta_map = {
            alias: self.root / (alias + '.fasta.gz')
            for alias in frame.get_column('sample_name').to_list()
            if alias not in missing_fastas
        }
        inspections = {alias: {'sequence_count': 2} for alias in fasta_map}
        queued = Mock()
        queued.json.return_value = {'_links': {'poll': {'href': 'https://example.invalid/poll'}}}
        batches = []

        def build_xml(rows, submission, project):
            batches.append([row['sample_name'] for row in rows])
            return b'<WEBIN/>'

        with patch.object(mag2ena, 'validate_submission_taxonomy_records', return_value=[]), \
             patch.object(mag2ena, 'inspect_fasta_collection', return_value=inspections), \
             patch.object(mag2ena, '_build_metadata_xml', side_effect=build_xml), \
             patch.object(mag2ena.requests, 'post', return_value=queued) as post, \
             patch.object(mag2ena, '_poll_receipt', side_effect=receipts), \
             patch.object(mag2ena, '_deduplicate_fasta_headers', side_effect=lambda path: (path, False)), \
             patch.object(mag2ena.subprocess, 'run',
                          side_effect=java_results,
                          return_value=subprocess.CompletedProcess([], 0, 'Submitted successfully.', '')) as java, \
             redirect_stdout(io.StringIO()):
            result = mag2ena.build_and_submit(frame, {
                'study_accession': 'PRJEB1', 'ena_user': 'example_user',
                'ena_password': 'example_password', 'portal': 'test',
            }, fasta_map, self.root / 'results')
        saved = json.loads((self.root / 'results/result.json').read_text())
        self.assertEqual(saved, result)
        error_path = self.root / 'results/error.txt'
        return result, error_path.read_text() if error_path.exists() else '', java.call_count, post.call_count, batches

    def test_complete_receipt_counts_every_requested_mag(self):
        frame = self.frame()
        result, errors, java, _, _ = self.run_submission(frame, [receipt(frame['sample_name'])])
        self.assertEqual(result['expected'], 2)
        self.assertEqual(result['submitted'], 2)
        self.assertEqual(result['errors'], 0)
        self.assertEqual(result['metadata_errors'], 0)
        self.assertTrue(result['complete'])
        self.assertEqual(java, 2)
        self.assertEqual(errors, '')

    def test_missing_accession_is_a_failed_mag(self):
        result, errors, java, _, _ = self.run_submission(self.frame(), [receipt(['mag_0'])])
        self.assertEqual(result['submitted'], 1)
        self.assertEqual(result['errors'], 1)
        self.assertFalse(result['complete'])
        self.assertEqual(java, 1)
        self.assertIn('mag_1', errors)
        self.assertIn('webin_log.xml', errors)

    def test_all_metadata_rejected_never_reports_success(self):
        result, errors, java, _, _ = self.run_submission(
            self.frame(), [receipt([], ['Metadata rejected'], success=False)])
        self.assertEqual(result['submitted'], 0)
        self.assertEqual(result['errors'], 2)
        self.assertEqual(result['metadata_errors'], 1)
        self.assertFalse(result['complete'])
        self.assertEqual(java, 0)
        self.assertIn('Metadata rejected', errors)

    def test_failed_receipt_without_error_text_blocks_complete_status(self):
        frame = self.frame()
        result, errors, _, _, _ = self.run_submission(frame, [receipt(frame['sample_name'], success=False)])
        self.assertEqual(result['metadata_errors'], 1)
        self.assertFalse(result['complete'])
        self.assertIn('unsuccessful', errors)

    def test_existing_alias_accession_can_be_reused(self):
        message = 'The sample already exists, alias: "mag_1", accession: "ERS_previous"'
        result, _, java, _, _ = self.run_submission(
            self.frame(), [receipt(['mag_0'], [message], success=False)])
        self.assertEqual(result['submitted'], 2)
        self.assertEqual(result['errors'], 0)
        self.assertEqual(result['metadata_errors'], 0)
        self.assertTrue(result['complete'])
        self.assertEqual(java, 2)

    def test_receipt_alias_outside_the_batch_is_not_processed(self):
        result, _, java, _, _ = self.run_submission(self.frame(), [receipt(['foreign_mag'])])
        self.assertEqual(result['errors'], 2)
        self.assertFalse(result['complete'])
        self.assertEqual(java, 0)

    def test_duplicate_aliases_are_rejected_before_network_calls(self):
        frame = self.frame().with_columns(pl.lit('same_alias').alias('sample_name'))
        with patch.object(mag2ena, 'validate_submission_taxonomy_records', return_value=[]), \
             patch.object(mag2ena, '_build_metadata_xml', return_value=b'<WEBIN/>'), \
             patch.object(mag2ena.requests, 'post', side_effect=AssertionError('Unexpected network call')) as post, \
             self.assertRaisesRegex(ValueError, 'unique'):
            mag2ena.build_and_submit(frame, {
                'study_accession': 'PRJEB1', 'ena_user': 'example_user',
                'ena_password': 'example_password', 'portal': 'test',
            }, {}, self.root / 'results')
        post.assert_not_called()

    def test_empty_metadata_are_rejected_before_network_calls(self):
        with patch.object(mag2ena, 'validate_submission_taxonomy_records', return_value=[]), \
             patch.object(mag2ena.requests, 'post', side_effect=AssertionError('Unexpected network call')) as post, \
             self.assertRaisesRegex(ValueError, 'nonempty'):
            mag2ena.build_and_submit(self.frame().head(0), {
                'study_accession': 'PRJEB1', 'ena_user': 'example_user',
                'ena_password': 'example_password', 'portal': 'test',
            }, {}, self.root / 'results')
        post.assert_not_called()

    def test_assembly_failure_counts_once(self):
        frame = self.frame()
        result, _, java, _, _ = self.run_submission(frame, [receipt(frame['sample_name'])], java_results=[
            subprocess.CompletedProcess([], 0, 'Submitted successfully.', ''),
            subprocess.CompletedProcess([], 1, '', 'Assembly rejected'),
        ])
        self.assertEqual(result['submitted'], 1)
        self.assertEqual(result['errors'], 1)
        self.assertFalse(result['complete'])
        self.assertEqual(java, 2)

    def test_missing_fasta_is_reported_instead_of_silently_skipped(self):
        frame = self.frame()
        result, errors, java, _, _ = self.run_submission(
            frame, [receipt(frame['sample_name'])], missing_fastas={'mag_1'})
        self.assertEqual(result['submitted'], 1)
        self.assertEqual(result['errors'], 1)
        self.assertFalse(result['complete'])
        self.assertEqual(java, 1)
        self.assertIn('No matching FASTA', errors)

    def test_251_mags_use_two_batches_with_complete_accounting(self):
        frame = self.frame(251)
        aliases = frame.get_column('sample_name').to_list()
        result, _, java, posts, batches = self.run_submission(
            frame, [receipt(aliases[:250]), receipt(aliases[250:])])
        self.assertEqual([len(batch) for batch in batches], [250, 1])
        self.assertEqual(result['submitted'], 251)
        self.assertTrue(result['complete'])
        self.assertEqual(java, 251)
        self.assertEqual(posts, 2)
        self.assertTrue((self.root / 'results/submit_250.xml').is_file())
        self.assertTrue((self.root / 'results/webin_log_250.xml').is_file())

    def test_main_returns_failure_when_metadata_errors_remain(self):
        frame = self.frame()
        metadata = self.root / 'input.tsv'
        frame.write_csv(metadata, separator='\t')
        fasta_map = {alias: self.root / (alias + '.fasta.gz') for alias in frame['sample_name']}
        inspections = {alias: {'sequence_count': 2} for alias in fasta_map}
        with patch.object(mag2ena, 'validate_dataframe'), \
             patch.object(mag2ena, 'collect_fastas', return_value=fasta_map), \
             patch.object(mag2ena, '_validate_fasta_headers', return_value=[]), \
             patch.object(mag2ena, 'inspect_fasta_collection', return_value=inspections), \
             patch.object(mag2ena, 'resolve_webin_credentials', return_value=('example_user', 'example_password')), \
             patch.object(mag2ena, 'build_and_submit', return_value={
                 'submitted': 2, 'errors': 0, 'metadata_errors': 1, 'complete': False,
             }), redirect_stdout(io.StringIO()):
            status = mag2ena.main([
                '--metadata', str(metadata), '--fasta-dir', str(self.root),
                '--study-accession', 'PRJEB1',
            ])
        self.assertEqual(status, 1)


if __name__ == '__main__':
    unittest.main()