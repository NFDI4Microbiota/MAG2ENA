import copy
import csv
import hashlib
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


REFERENCE = 'GCF_000005845.2'
CLASSIFICATION = 'd__Bacteria;p__Pseudomonadota;g__Escherichia_A;s__Escherichia coli_A'


def node(taxid='562', name='Escherichia coli', **changes):
    record = {
        'taxId': taxid, 'scientificName': name, 'rank': 'species',
        'submittable': 'true', 'binomial': 'true', 'metagenome': 'false',
        'lineage': 'cellular organisms;Bacteria;Pseudomonadota;',
        'otherNames': [],
    }
    record.update(changes)
    return record


def saved_row(prefix='fastani_', classification=CLASSIFICATION):
    radius_field = prefix + ('radius' if prefix == 'closest_placement_' else 'reference_radius')
    return {
        'user_genome': 'mag_1', 'classification': classification,
        prefix + 'reference': REFERENCE, radius_field: '95.0',
        prefix + 'taxonomy': classification, prefix + 'ani': '98.5',
        prefix + 'af': '0.75', 'warnings': 'N/A', 'note': 'original_note_A',
        'extra_original_column': 'keep_this_value',
    }


class TaxonomyResolutionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.cache = self.root / 'cache.sqlite'

    def tearDown(self):
        self.temporary.cleanup()

    def frame(self, organism='Escherichia coli', taxid=''):
        return pd.DataFrame([{
            'sample_name': 'mag_1', 'organism': organism, 'tax_id': taxid,
            'untouched': '001_A', 'nullable': None,
        }])

    def write_summary(self, row=None, filename='summary.tsv'):
        row = row or saved_row()
        path = self.root / filename
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row), delimiter='\t')
            writer.writeheader()
            writer.writerow(row)
        return path

    def metadata(self, classification=CLASSIFICATION, reference=REFERENCE, **changes):
        row = {
            'accession': 'RS_' + reference, 'gtdb_taxonomy': classification,
            'gtdb_genome_representative': 'RS_' + reference,
            'ncbi_taxid': '562', 'ncbi_species_taxid': '562',
            'ncbi_organism_name': 'Escherichia coli',
            'ncbi_genbank_assembly_accession': '',
        }
        row.update(changes)
        path = self.root / 'bac120_metadata_r214.tsv'
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row), delimiter='\t')
            writer.writeheader()
            writer.writerow(row)
        return path

    def run_batch(self, frame=None, payloads=None, **options):
        frame = self.frame() if frame is None else frame
        if payloads is None:
            payloads = {
                taxonomy.ENA + 'any-name/Escherichia%20coli': [node()],
                taxonomy.ENA + 'tax-id/562': node(),
            }

        def lookup(cache, url, kind='json'):
            value = payloads.get(url)
            if isinstance(value, Exception):
                raise value
            return copy.deepcopy(value)

        with patch.object(taxonomy.TaxonomyCache, 'get', new=lookup):
            return taxonomy.resolve_taxonomy_batch(frame, cache_path=self.cache, **options)

    def reference_batch(self, *, raw=None, classification=CLASSIFICATION,
                        metadata_changes=None, payloads=None, frame=None, **options):
        raw = raw or saved_row(classification=classification)
        summary = self.write_summary(raw)
        metadata = self.metadata(classification=classification, **(metadata_changes or {}))
        if payloads is None:
            payloads = {
                taxonomy.NCBI + REFERENCE + '/dataset_report': {'reports': [{
                    'accession': REFERENCE, 'organism': {
                        'tax_id': 562, 'organism_name': 'Escherichia coli',
                    },
                }]},
                taxonomy.ENA + 'tax-id/562': node(),
            }
        return self.run_batch(
            frame=frame, payloads=payloads,
            evidence=taxonomy.read_gtdb_summaries([summary]), release='214.1',
            metadata_paths=[metadata], **options,
        )

    def test_exact_name_fills_consistent_pair_and_preserves_other_values(self):
        original = self.frame()
        result, report, summary = self.run_batch(original)
        self.assertEqual(result.loc[0, 'tax_id'], '562')
        self.assertEqual(summary['resolved'], 1)
        pd.testing.assert_frame_equal(result.drop(columns=['organism', 'tax_id']),
                                      original.drop(columns=['organism', 'tax_id']))
        self.assertEqual(original.loc[0, 'tax_id'], '')
        self.assertEqual(report[0]['applied'], 'yes')

    def test_nonexact_suggestion_is_rejected(self):
        result, report, _ = self.run_batch(payloads={
            taxonomy.ENA + 'any-name/Escherichia%20coli': [node(name='Escherichia fergusonii')],
        })
        self.assertEqual(result.loc[0, 'tax_id'], '')
        self.assertEqual(report[0]['status'], 'NAME_ID_MISMATCH_OR_UNRESOLVED')

    def test_exact_synonym_updates_to_current_scientific_name(self):
        record = node(otherNames=['Oldgenus coli:synonym'])
        result, _, _ = self.run_batch(self.frame('Oldgenus coli'), payloads={
            taxonomy.ENA + 'any-name/Oldgenus%20coli': [record],
        })
        self.assertEqual(tuple(result.loc[0, ['organism', 'tax_id']]), ('Escherichia coli', '562'))

    def test_ambiguous_exact_synonym_is_rejected(self):
        a = node(otherNames=['Oldgenus coli:synonym'])
        b = node('999', 'Othergenus coli', otherNames=['Oldgenus coli:synonym'])
        _, report, _ = self.run_batch(self.frame('Oldgenus coli'), payloads={
            taxonomy.ENA + 'any-name/Oldgenus%20coli': [a, b],
        })
        self.assertEqual(report[0]['status'], 'AMBIGUOUS_EXACT_NAME')

    def test_family_and_genus_cannot_supply_a_species(self):
        for rank in ('family', 'genus'):
            with self.subTest(rank=rank):
                _, report, summary = self.run_batch(payloads={
                    taxonomy.ENA + 'any-name/Escherichia%20coli': [node(rank=rank)],
                })
                self.assertEqual(len(summary['unresolved']), 1)
                self.assertIn('TAXON_NOT_SPECIES_RANK', report[0]['detail'])

    def test_subspecies_resolves_only_to_unique_species_ancestor(self):
        narrow = node('77', 'Escherichia coli subspecies alpha', rank='subspecies')
        result, _, _ = self.run_batch(
            self.frame(narrow['scientificName'], '77'), scope='all', payloads={
                taxonomy.ENA + 'tax-id/77': narrow,
                taxonomy.ENA_XML + '77': '<TAXON_SET><taxon><lineage><taxon taxId="562" rank="species"/></lineage></taxon></TAXON_SET>',
                taxonomy.ENA + 'tax-id/562': node(),
            })
        self.assertEqual(result.loc[0, 'tax_id'], '562')

    def test_generic_replacement_can_be_disabled(self):
        generic = node('999', 'Escherichia sp.', binomial='false')
        payloads = {taxonomy.ENA + 'any-name/Escherichia%20sp.': [generic]}
        result, report, summary = self.run_batch(self.frame('Escherichia sp.'), payloads=payloads, apply_generic=False)
        self.assertEqual(result.loc[0, 'tax_id'], '')
        self.assertIn('--apply-generic-taxonomy', report[0]['blocking_reason'])
        result, _, summary = self.run_batch(
            self.frame('Escherichia sp.'), payloads=payloads)
        self.assertEqual(summary['resolved'], 1)
        self.assertEqual(result.loc[0, 'tax_id'], '999')

    def test_existing_verified_generic_pair_needs_no_replacement_permission(self):
        record = node('999', 'Escherichia sp.', binomial='false')
        _, report, summary = self.run_batch(self.frame('Escherichia sp.', '999'), payloads={
            taxonomy.ENA + 'tax-id/999': record,
        })
        self.assertEqual(summary['resolved'], 1)
        self.assertEqual(report[0]['applied'], 'no')

    def test_metagenome_node_is_rejected_even_with_generic_option(self):
        _, _, summary = self.run_batch(apply_generic=True, payloads={
            taxonomy.ENA + 'any-name/Escherichia%20coli': [node(metagenome='true')],
        })
        self.assertEqual(len(summary['unresolved']), 1)

    def test_missing_submittable_flag_is_not_accepted(self):
        record = node()
        del record['submittable']
        _, report, _ = self.run_batch(payloads={
            taxonomy.ENA + 'any-name/Escherichia%20coli': [record],
        })
        self.assertIn('TAXON_NOT_SUBMITTABLE', report[0]['detail'])

    def test_populated_id_is_not_overwritten_in_empty_scope(self):
        original = self.frame('Escherichia coli', '999')
        result, report, summary = self.run_batch(original)
        pd.testing.assert_frame_equal(result, original)
        self.assertIn('--taxonomy-scope all', report[0]['blocking_reason'])
        self.assertEqual(len(summary['unresolved']), 1)

    def test_all_scope_permits_supported_correction(self):
        result, _, summary = self.run_batch(self.frame('Escherichia coli', '999'), scope='all')
        self.assertEqual(result.loc[0, 'tax_id'], '562')
        self.assertEqual(summary['applied'], 1)

    def test_gtdb_suffix_uses_recorded_reference_and_preserves_raw_evidence(self):
        raw = saved_row()
        result, report, summary = self.reference_batch(
            raw=raw, frame=self.frame('Escherichia coli_A', '999'), scope='all')
        self.assertEqual(summary['resolved'], 1)
        self.assertEqual(result.loc[0, 'organism'], 'Escherichia coli')
        self.assertEqual(report[0]['reference_accession'], REFERENCE)
        self.assertEqual(json.loads(report[0]['gtdb_record_json']), raw)
        source = json.loads(report[0]['gtdb_sources_json'])[0]
        self.assertEqual(source['sha256'], taxonomy.sha_file(Path(source['path'])))

    def test_every_saved_gtdb_species_uses_reference_evidence(self):
        classification = CLASSIFICATION.replace('coli_A', 'coli')
        _, report, _ = self.reference_batch(classification=classification)
        self.assertEqual(report[0]['reference_accession'], REFERENCE)

    def test_both_domains_and_reordered_metadata_preserve_identity(self):
        bacteria = self.write_summary()
        bacterial_metadata = self.metadata()
        archaeal_classification = 'd__Archaea;g__Methanobrevibacter;s__Methanobrevibacter smithii_A'
        archaeal = saved_row('closest_genome_', archaeal_classification)
        archaeal['user_genome'] = 'mag_archaea'
        archaeal['closest_genome_reference'] = 'GCF_000016525.1'
        archaeal_summary = self.write_summary(archaeal, 'archaea.tsv')
        archaeal_metadata = self.root / 'ar53_metadata_r214.tsv'
        with bacterial_metadata.open() as handle:
            metadata = next(csv.DictReader(handle, delimiter='\t'))
        metadata.update(
            accession='RS_GCF_000016525.1', gtdb_taxonomy=archaeal_classification,
            gtdb_genome_representative='RS_GCF_000016525.1',
            ncbi_taxid='2173', ncbi_species_taxid='2173',
            ncbi_organism_name='Methanobrevibacter smithii',
        )
        with archaeal_metadata.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(metadata), delimiter='\t')
            writer.writeheader()
            writer.writerow(metadata)
        original = pd.concat([
            self.frame('Methanobrevibacter smithii_A').assign(sample_name='mag_archaea'),
            self.frame('Escherichia coli_A'),
        ], ignore_index=True)
        result, reports, summary = self.run_batch(
            original,
            evidence=taxonomy.read_gtdb_summaries([bacteria, archaeal_summary]),
            release='214.1', metadata_paths=[bacterial_metadata, archaeal_metadata],
            payloads={
                taxonomy.NCBI + REFERENCE + '/dataset_report': {
                    'reports': [{'accession': REFERENCE, 'organism': {'tax_id': 562}}]},
                taxonomy.NCBI + 'GCF_000016525.1/dataset_report': {
                    'reports': [{'accession': 'GCF_000016525.1', 'organism': {'tax_id': 2173}}]},
                taxonomy.ENA + 'tax-id/562': node(),
                taxonomy.ENA + 'tax-id/2173': node(
                    '2173', 'Methanobrevibacter smithii', lineage='cellular organisms;Archaea;'),
            },
        )
        self.assertEqual(summary['resolved'], 2)
        self.assertEqual(list(result['sample_name']), ['mag_archaea', 'mag_1'])
        self.assertEqual(list(result['tax_id']), ['2173', '562'])
        self.assertEqual(reports[0]['gtdb_classification'], archaeal_classification)
        pd.testing.assert_frame_equal(result.drop(columns=['organism', 'tax_id']),
                                      original.drop(columns=['organism', 'tax_id']))

    def test_gtdb_markers_require_evidence_when_summary_is_missing(self):
        for organism in ('Escherichia coli_A', 'UBA123 sp000005845'):
            with self.subTest(organism=organism):
                _, report, _ = self.run_batch(self.frame(organism), payloads={})
                self.assertEqual(report[0]['status'], 'GTDB_EVIDENCE_REQUIRED')

    def test_import_preserves_underscores_and_omits_higher_rank_fallback(self):
        path = self.write_summary()
        evidence = taxonomy.read_gtdb_summaries([path])
        self.assertEqual(taxonomy.gtdb_organism_table(evidence).loc[0, 'organism'],
                         'Escherichia coli_A')
        row = saved_row(classification='d__Bacteria;g__Escherichia_A;s__')
        evidence = taxonomy.read_gtdb_summaries([self.write_summary(row, 'higher.tsv')])
        self.assertEqual(taxonomy.gtdb_organism_table(evidence).loc[0, 'organism'], '')
        _, report, _ = self.run_batch(evidence=evidence, payloads={})
        self.assertEqual(report[0]['status'], 'NO_SPECIES_ASSIGNMENT')

    def test_release_species_mismatch_is_blocked(self):
        raw = saved_row()
        metadata = self.metadata(classification=CLASSIFICATION.replace('coli_A', 'coli_B'))
        _, report, _ = self.run_batch(
            evidence=taxonomy.read_gtdb_summaries([self.write_summary(raw)]),
            release='214.1', metadata_paths=[metadata], payloads={})
        self.assertEqual(report[0]['status'], 'GTDB_RELEASE_SPECIES_MISMATCH')

    def test_missing_reference_is_blocked(self):
        _, report, _ = self.reference_batch(
            metadata_changes={'accession': 'RS_GCF_999999999.1'}, payloads={})
        self.assertEqual(report[0]['status'], 'REFERENCE_ABSENT_FROM_SELECTED_RELEASE')

    def test_below_threshold_and_nonfinite_metrics_are_blocked(self):
        for field, value in [('fastani_ani', '90'), ('fastani_af', '0.1'), ('fastani_ani', 'nan')]:
            with self.subTest(field=field, value=value):
                raw = saved_row()
                raw[field] = value
                _, _, summary = self.reference_batch(raw=raw, payloads={})
                self.assertEqual(len(summary['unresolved']), 1)

    def test_stricter_original_af_is_honoured(self):
        raw = saved_row()
        raw['fastani_af'] = '0.55'
        _, report, _ = self.reference_batch(raw=raw, min_af=0.65, payloads={})
        self.assertEqual(report[0]['status'], 'REFERENCE_BELOW_SPECIES_THRESHOLD')

    def test_modern_and_placement_reference_fields_are_supported(self):
        for prefix in ('closest_genome_', 'closest_placement_'):
            with self.subTest(prefix=prefix):
                _, report, summary = self.reference_batch(raw=saved_row(prefix))
                self.assertEqual(summary['resolved'], 1)
                self.assertEqual(report[0]['reference_evidence'], prefix.rstrip('_'))

    def test_snapshot_only_mapping_is_not_applied(self):
        _, report, summary = self.reference_batch(payloads={
            taxonomy.NCBI + REFERENCE + '/dataset_report': {'reports': []},
            taxonomy.ENA + 'tax-id/562': node(),
        })
        self.assertEqual(report[0]['status'], 'SNAPSHOT_ONLY_REFERENCE_MAPPING')
        self.assertEqual(report[0]['candidate_verified'], 'no')
        self.assertEqual(len(summary['unresolved']), 1)

    def test_wrong_domain_is_blocked(self):
        payloads = {
            taxonomy.NCBI + REFERENCE + '/dataset_report': {'reports': [{
                'accession': REFERENCE, 'organism': {'tax_id': 562},
            }]},
            taxonomy.ENA + 'tax-id/562': node(lineage='cellular organisms;Archaea;'),
        }
        _, report, _ = self.reference_batch(payloads=payloads)
        self.assertEqual(report[0]['status'], 'DOMAIN_MISMATCH')

    def test_duplicate_and_conflicting_sample_keys_are_rejected(self):
        original = self.frame()
        with self.assertRaisesRegex(ValueError, 'unique'):
            self.run_batch(pd.concat([original, original]))
        first = self.write_summary()
        raw = saved_row()
        raw['fastani_ani'] = '99'
        second = self.write_summary(raw, 'second.tsv')
        with self.assertRaisesRegex(ValueError, 'Conflicting'):
            taxonomy.read_gtdb_summaries([first, second])

    def test_service_failure_on_last_row_still_blocks(self):
        frame = pd.concat([self.frame().assign(sample_name=f'mag_{i}') for i in range(5)])
        _, report, summary = self.run_batch(frame, payloads={
            taxonomy.ENA + 'any-name/Escherichia%20coli':
                taxonomy.Unresolved('LOOKUP_FAILED', 'HTTP 503'),
        })
        self.assertEqual(len(summary['unresolved']), 5)
        self.assertTrue(all(row['status'] == 'LOOKUP_FAILED' for row in report))

    def test_report_serializes_full_evidence_and_protects_inputs(self):
        _, reports, _ = self.reference_batch()
        path = self.root / 'report.tsv'
        taxonomy.write_taxonomy_report(path, reports)
        with path.open() as handle:
            row = next(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(json.loads(row['gtdb_record_json'])['extra_original_column'],
                         'keep_this_value')
        with self.assertRaisesRegex(ValueError, 'overwrite'):
            taxonomy.write_taxonomy_report(path, reports, protected_paths=[path])

    def test_changed_local_metadata_cannot_silently_reuse_old_index(self):
        metadata = self.metadata()
        cache = taxonomy.TaxonomyCache(self.cache)
        try:
            taxonomy.index_metadata(cache, '214.1', metadata)
            metadata.write_text(metadata.read_text().replace('Escherichia coli_A', 'Other species'))
            with self.assertRaisesRegex(ValueError, 'changed'):
                taxonomy.index_metadata(cache, '214.1', metadata)
        finally:
            cache.db.close()

    def test_offline_http_cache_checks_its_digest(self):
        cache = taxonomy.TaxonomyCache(self.cache, offline=True)
        try:
            body = json.dumps(node())
            url = taxonomy.ENA + 'tax-id/562'
            cache.db.execute('INSERT INTO http VALUES (?,?,?,?,?)',
                             (url, 200, body, hashlib.sha256(body.encode()).hexdigest(), 0))
            cache.db.commit()
            self.assertEqual(cache.ena('562')['scientificName'], 'Escherichia coli')
            cache.db.execute('UPDATE http SET sha=?', ('bad',))
            cache.db.commit()
            with self.assertRaisesRegex(ValueError, 'checksum'):
                cache.ena('562')
        finally:
            cache.db.close()

    def test_refresh_retains_indexed_metadata(self):
        metadata = self.metadata()
        cache = taxonomy.TaxonomyCache(self.cache)
        taxonomy.index_metadata(cache, '214.1', metadata)
        cache.db.close()
        refreshed = taxonomy.TaxonomyCache(self.cache, refresh=True)
        try:
            self.assertEqual(refreshed.source_for('214.1', REFERENCE)['species'],
                             'Escherichia coli_A')
        finally:
            refreshed.db.close()

    def test_final_validator_checks_name_and_exact_species_rank(self):
        valid = self.frame(taxid='562').to_dict('records')
        self.assertEqual(taxonomy.validate_submission_taxonomy_records(valid, fetcher=lambda _: node()), [])
        for rank in ('genus', 'family', 'subspecies'):
            with self.subTest(rank=rank):
                errors = taxonomy.validate_submission_taxonomy_records(
                    valid, fetcher=lambda _: node(rank=rank))
                self.assertIn('TAXON_NOT_SPECIES_RANK', errors[0]['expected'])
        wrong = self.frame('Escherichia coli_A', '562').to_dict('records')
        self.assertTrue(taxonomy.validate_submission_taxonomy_records(wrong, fetcher=lambda _: node()))

    def test_cli_blocks_before_credentials_and_submission(self):
        metadata = self.root / 'input.tsv'
        self.frame('Escherichia coli_A').to_csv(metadata, sep='\t', index=False)
        report = self.root / 'report.tsv'
        prepared = self.root / 'prepared.tsv'
        with patch.object(taxonomy.TaxonomyCache, 'get', return_value=None), \
             patch.object(mag2ena, 'resolve_webin_credentials') as credentials, \
             patch.object(mag2ena, 'build_and_submit') as submit, redirect_stdout(io.StringIO()) as output:
            status = mag2ena.main([
                '--metadata', str(metadata), '--resolve-taxonomy',
                '--taxonomy-report-out', str(report),
                '--prepared-metadata-out', str(prepared),
                '--taxonomy-cache', str(self.cache),
            ])
        self.assertEqual(status, 1)
        self.assertIn('GTDB_EVIDENCE_REQUIRED', output.getvalue())
        credentials.assert_not_called()
        submit.assert_not_called()
        self.assertTrue(report.is_file())
        self.assertTrue(prepared.is_file())

    def test_cli_prepare_only_returns_nonzero_on_unresolved_taxonomy(self):
        metadata = self.root / 'input.tsv'
        self.frame('Escherichia coli_A').to_csv(metadata, sep='\t', index=False)
        with patch.object(taxonomy.TaxonomyCache, 'get', return_value=None), redirect_stdout(io.StringIO()):
            status = mag2ena.main([
                '--metadata', str(metadata), '--resolve-taxonomy', '--prepare-only',
                '--prepared-metadata-out', str(self.root / 'prepared.tsv'),
                '--taxonomy-cache', str(self.cache),
            ])
        self.assertEqual(status, 1)

    def test_cli_prevents_output_input_collision(self):
        metadata = self.root / 'input.tsv'
        self.frame().to_csv(metadata, sep='\t', index=False)
        before = metadata.read_bytes()
        with redirect_stdout(io.StringIO()):
            status = mag2ena.main([
                '--metadata', str(metadata), '--prepared-metadata-out', str(metadata),
                '--prepare-only',
            ])
        self.assertEqual(status, 2)
        self.assertEqual(metadata.read_bytes(), before)

    def test_cli_prevents_cache_input_collision(self):
        metadata = self.root / 'input.tsv'
        self.frame().to_csv(metadata, sep='\t', index=False)
        before = metadata.read_bytes()
        with redirect_stdout(io.StringIO()):
            status = mag2ena.main([
                '--metadata', str(metadata), '--resolve-taxonomy', '--prepare-only',
                '--taxonomy-cache', str(metadata),
                '--prepared-metadata-out', str(self.root / 'prepared.tsv'),
            ])
        self.assertEqual(status, 2)
        self.assertEqual(metadata.read_bytes(), before)

    def test_cli_without_resolver_blocks_mismatched_name_before_credentials(self):
        metadata = self.root / 'input.tsv'
        self.frame('Escherichia coli_A', '562').to_csv(metadata, sep='\t', index=False)
        with patch.object(taxonomy.TaxonomyCache, 'get', return_value=node()), \
             patch.object(mag2ena, 'load_fields_from_xml', return_value={}), \
             patch.object(mag2ena, 'resolve_webin_credentials') as credentials, \
             patch.object(mag2ena, 'build_and_submit') as submit, \
             redirect_stdout(io.StringIO()) as output:
            status = mag2ena.main([
                '--metadata', str(metadata), '--fasta-dir', str(self.root),
                '--study-accession', 'PRJEB_TEST',
            ])
        self.assertEqual(status, 2)
        self.assertIn('does not match', output.getvalue())
        credentials.assert_not_called()
        submit.assert_not_called()

    def test_cli_successful_preparation_writes_consistent_pair(self):
        metadata = self.root / 'input.tsv'
        self.frame().to_csv(metadata, sep='\t', index=False)
        prepared = self.root / 'prepared.tsv'
        with patch.object(taxonomy.TaxonomyCache, 'get', return_value=[node()]), \
             redirect_stdout(io.StringIO()):
            status = mag2ena.main([
                '--metadata', str(metadata), '--resolve-taxonomy', '--prepare-only',
                '--prepared-metadata-out', str(prepared),
                '--taxonomy-cache', str(self.cache),
            ])
        self.assertEqual(status, 0)
        output = pd.read_csv(prepared, sep='\t', dtype=str, keep_default_na=False)
        self.assertEqual(output.loc[0, 'tax_id'], '562')

    def test_cli_accepts_repeated_gtdb_summaries(self):
        args = mag2ena.build_parser().parse_args([
            '--metadata', 'input.tsv', '--gtdbtk', 'bacteria.tsv', '--gtdbtk', 'archaea.tsv',
        ])
        self.assertEqual(args.gtdbtk, [Path('bacteria.tsv'), Path('archaea.tsv')])


if __name__ == '__main__':
    unittest.main()
