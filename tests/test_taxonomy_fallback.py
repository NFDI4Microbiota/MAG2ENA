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
from urllib.parse import quote

import pandas as pd

import mag2ena
from App.utils import taxonomy_resolution as taxonomy


def node(taxid, name, **changes):
    record = {
        'taxId': str(taxid), 'scientificName': name, 'rank': 'species',
        'submittable': 'true', 'binomial': 'false', 'metagenome': 'false',
        'lineage': 'Bacteria; Testfamily; Testgenus;', 'otherNames': [],
    }
    record.update(changes)
    return record


def query_url(name, endpoint='any-name'):
    return taxonomy.ENA + endpoint + '/' + quote(name, safe='')


class TaxonomyFallbackTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.urls = []

    def tearDown(self):
        self.temporary.cleanup()

    def frame(self, organism='uncultured Testgenus bacterium', taxid=''):
        return pd.DataFrame([{'sample_name': 'mag_1', 'organism': organism,
                              'tax_id': taxid, 'keep': '001_A'}])

    def evidence(self, classification='d__Bacteria;f__Testfamily;g__Testgenus;s__', **extra):
        return {'mag_1': {'row': {'user_genome': 'mag_1', 'classification': classification,
                                 **extra}, 'sources': [{'path': '/original.tsv', 'sha256': 'original'}]}}

    def lookup(self, payloads):
        def get(cache, url, kind='json'):
            self.urls.append(url)
            value = payloads.get(url)
            if isinstance(value, Exception):
                raise value
            return copy.deepcopy(value)
        return get

    def payloads(self, name='Testgenus sp.', taxid='111', **changes):
        record = node(taxid, name, **changes)
        return {query_url(name): [record], taxonomy.ENA + 'tax-id/' + str(taxid): record}

    def run_batch(self, *, frame=None, evidence=None, payloads=None, **options):
        frame = self.frame() if frame is None else frame
        evidence = self.evidence() if evidence is None else evidence
        payloads = self.payloads() if payloads is None else payloads
        with patch.object(taxonomy.TaxonomyCache, 'get', new=self.lookup(payloads)):
            return taxonomy.resolve_taxonomy_batch(
                frame, evidence=evidence, release='214.1',
                cache_path=self.root / 'cache.sqlite', **options,
            )

    def test_genus_fallback_is_applied_by_default(self):
        result, report, summary = self.run_batch()
        self.assertEqual(summary['resolved'], 1)
        self.assertEqual(tuple(result.loc[0, ['organism', 'tax_id']]), ('Testgenus sp.', '111'))
        self.assertEqual(report[0]['primary_status'], 'NAME_ID_MISMATCH_OR_UNRESOLVED')
        self.assertEqual(report[0]['mapping_scope'], 'GENUS_LEVEL_ENVIRONMENTAL_TAXON')
        self.assertEqual(report[0]['fallback_source_rank'], 'genus')

    def test_genus_fallback_can_be_left_as_a_proposal(self):
        original = self.frame()
        result, report, summary = self.run_batch(frame=original, apply_generic=False)
        pd.testing.assert_frame_equal(result, original)
        self.assertEqual(len(summary['unresolved']), 1)
        self.assertEqual(report[0]['candidate_tax_id'], '111')
        self.assertIn('--no-apply-generic-taxonomy', report[0]['blocking_reason'])

    def test_genus_is_preferred_over_family(self):
        payloads = self.payloads()
        payloads.update(self.payloads('Testfamily bacterium', '222'))
        result, _, _ = self.run_batch(payloads=payloads)
        self.assertEqual(result.loc[0, 'tax_id'], '111')
        self.assertNotIn(query_url('Testfamily bacterium'), self.urls)

    def test_chamaesiphonaceae_family_fallback_accepts_cyanobacterium(self):
        evidence = self.evidence('d__Bacteria;p__Cyanobacteriota;f__Chamaesiphonaceae;g__JAAFHZ01;s__')
        result, report, summary = self.run_batch(
            frame=self.frame('uncultured JAAFHZ01 bacterium'), evidence=evidence,
            payloads=self.payloads('Chamaesiphonaceae cyanobacterium', '3107296', binomial='true'))
        self.assertEqual(summary['resolved'], 1)
        self.assertEqual(result.loc[0, 'tax_id'], '3107296')
        self.assertEqual(report[0]['fallback_source_rank'], 'family')
        self.assertEqual(report[0]['mapping_scope'], 'FAMILY_LEVEL_ENVIRONMENTAL_TAXON')

    def test_gemmatimonadaceae_family_fallback(self):
        result, report, _ = self.run_batch(
            frame=self.frame('uncultured AG11 bacterium'),
            evidence=self.evidence('d__Bacteria;f__Gemmatimonadaceae;g__AG11;s__'),
            payloads=self.payloads('Gemmatimonadaceae bacterium', '2026741', binomial='true'))
        self.assertEqual(result.loc[0, 'tax_id'], '2026741')
        self.assertEqual(report[0]['fallback_source_name'], 'Gemmatimonadaceae')

    def test_requested_cylindrospermum_mapping(self):
        result, _, summary = self.run_batch(
            frame=self.frame('uncultured Cylindrospermum bacterium'),
            evidence=self.evidence('d__Bacteria;p__Cyanobacteriota;g__Cylindrospermum;s__'),
            payloads=self.payloads('Cylindrospermum sp.', '2910429'))
        self.assertEqual(summary['resolved'], 1)
        self.assertEqual(result.loc[0, 'tax_id'], '2910429')

    def test_identifier_lookup_rejects_nearby_suggestions(self):
        exact = node('1948624', 'Gammaproteobacteria bacterium UBA4582', binomial='true')
        nearby = node('999', 'Gammaproteobacteria bacterium UBA4588', binomial='true')
        result, report, _ = self.run_batch(
            frame=self.frame('uncultured UBA4582 bacterium'),
            evidence=self.evidence('d__Bacteria;g__UBA4582;s__'), payloads={
                query_url('UBA4582', 'suggest-for-search'): [nearby, exact],
                taxonomy.ENA + 'tax-id/1948624': exact,
            })
        self.assertEqual(result.loc[0, 'tax_id'], '1948624')
        self.assertEqual(report[0]['mapping_scope'], 'IDENTIFIER_NAME_MATCH_ONLY')
        self.assertNotIn(taxonomy.ENA + 'tax-id/999', self.urls)

    def test_family_identifier_mapping_is_reported_as_family_evidence(self):
        record = node('1953107', 'Acidobacteria bacterium UBA2999', binomial='true')
        result, report, _ = self.run_batch(
            frame=self.frame('uncultured UBA2999 bacterium'),
            evidence=self.evidence('d__Bacteria;f__UBA2999;g__;s__'), payloads={
                query_url('UBA2999', 'suggest-for-search'): [record],
                taxonomy.ENA + 'tax-id/1953107': record,
            })
        self.assertEqual(result.loc[0, 'tax_id'], '1953107')
        self.assertEqual(report[0]['fallback_source_rank'], 'family')

    def test_ambiguous_identifier_is_not_selected(self):
        records = [node('101', 'Bacteria bacterium OLB17', binomial='true'),
                   node('102', 'Acidobacteria bacterium OLB17', binomial='true')]
        original = self.frame('uncultured OLB17 bacterium')
        result, report, summary = self.run_batch(
            frame=original, evidence=self.evidence('d__Bacteria;g__OLB17;s__'), payloads={
                query_url('OLB17', 'suggest-for-search'): records,
                taxonomy.ENA + 'tax-id/101': records[0], taxonomy.ENA + 'tax-id/102': records[1],
            })
        pd.testing.assert_frame_equal(result, original)
        self.assertEqual(len(summary['unresolved']), 1)
        self.assertEqual(report[0]['fallback_status'], 'NO_ENVIRONMENTAL_MAPPING')

    def test_unrelated_jake01_suggestion_is_not_selected(self):
        _, report, summary = self.run_batch(
            frame=self.frame('uncultured JAKEEG01 bacterium'),
            evidence=self.evidence('d__Bacteria;g__JAKEEG01;s__'), payloads={
                query_url('JAKEEG01', 'suggest-for-search'): [node('1526566', 'Planococcus sp. jake01')],
            })
        self.assertEqual(len(summary['unresolved']), 1)
        self.assertEqual(report[0]['candidate_tax_id'], '')

    def test_named_descendant_is_not_selected_as_a_generic_candidate(self):
        for record in (node('112', 'Testgenus namedii', binomial='true'),
                       node('112', 'Testgenus namedii', binomial='true', otherNames=['Testgenus sp.:synonym'])):
            with self.subTest(record=record):
                _, _, summary = self.run_batch(payloads={
                    query_url('Testgenus sp.'): [record], taxonomy.ENA + 'tax-id/112': record,
                })
                self.assertEqual(len(summary['unresolved']), 1)

    def test_invalid_environmental_nodes_remain_unresolved(self):
        for changes in ({'rank': 'genus'}, {'rank': 'family'}, {'rank': 'subspecies'},
                        {'submittable': 'false'}, {'metagenome': 'true'},
                        {'lineage': 'Archaea; Testgenus;'}):
            with self.subTest(changes=changes):
                _, _, summary = self.run_batch(payloads=self.payloads(**changes))
                self.assertEqual(len(summary['unresolved']), 1)

    def test_suffix_genus_fallback_preserves_raw_evidence(self):
        evidence = self.evidence('d__Bacteria;g__Olsenella_F;s__', extra='retain_A')
        result, report, _ = self.run_batch(
            frame=self.frame('uncultured Olsenella F bacterium'), evidence=evidence,
            payloads=self.payloads('Olsenella sp.', '2042683'))
        self.assertEqual(result.loc[0, 'organism'], 'Olsenella sp.')
        self.assertEqual(report[0]['fallback_source_name'], 'Olsenella_F')
        self.assertEqual(json.loads(report[0]['gtdb_record_json']), evidence['mag_1']['row'])
        self.assertEqual(json.loads(report[0]['gtdb_sources_json']), evidence['mag_1']['sources'])
        self.assertEqual(result.loc[0, 'keep'], '001_A')

    def test_suffix_name_is_canonicalized_for_unchanged_supplied_id(self):
        original = self.frame('Ruminococcus_E bromii_B', '40518')
        record = node('40518', 'Ruminococcus bromii', binomial='true')
        result, report, summary = self.run_batch(
            frame=original, evidence={}, payloads={taxonomy.ENA + 'tax-id/40518': record})
        self.assertEqual(summary['resolved'], 1)
        self.assertEqual(result.loc[0, 'organism'], 'Ruminococcus bromii')
        self.assertEqual(result.loc[0, 'tax_id'], '40518')
        self.assertEqual(report[0]['original_organism'], 'Ruminococcus_E bromii_B')
        self.assertEqual(taxonomy.validate_submission_taxonomy_records(
            result.to_dict('records'), fetcher=lambda _: record), [])

    def test_suffix_name_does_not_accept_an_unrelated_supplied_id(self):
        result, _, summary = self.run_batch(
            frame=self.frame('Ruminococcus_E bromii_B', '40518'), evidence={},
            payloads={taxonomy.ENA + 'tax-id/40518': node('40518', 'Othergenus differentii', binomial='true')})
        self.assertEqual(len(summary['unresolved']), 1)
        self.assertEqual(result.loc[0, 'organism'], 'Ruminococcus_E bromii_B')

    def test_generic_opt_out_still_allows_supplied_id_name_canonicalization(self):
        result, _, summary = self.run_batch(
            frame=self.frame('Lawsonella clevelandensis_A', '1528099'), evidence={},
            payloads={taxonomy.ENA + 'tax-id/1528099': node('1528099', 'Lawsonella clevelandensis', binomial='true')},
            apply_generic=False)
        self.assertEqual(summary['resolved'], 1)
        self.assertEqual(result.loc[0, 'organism'], 'Lawsonella clevelandensis')

    def test_populated_different_id_still_needs_scope_all(self):
        original = self.frame(taxid='999')
        result, report, summary = self.run_batch(frame=original)
        self.assertEqual(len(summary['unresolved']), 1)
        self.assertIn('--taxonomy-scope all', report[0]['blocking_reason'])
        pd.testing.assert_frame_equal(result, original)
        result, _, summary = self.run_batch(frame=original, scope='all')
        self.assertEqual(summary['resolved'], 1)
        self.assertEqual(result.loc[0, 'tax_id'], '111')

    def reference_fixture(self, state='current', ani='98'):
        classification = 'd__Bacteria;g__Testgenus;s__Testgenus sp000000001'
        evidence = self.evidence(classification,
            fastani_reference='GCF_000000001.1', fastani_reference_radius='95',
            fastani_taxonomy=classification, fastani_ani=ani, fastani_af='0.8')
        metadata = {'species': 'Testgenus sp000000001', 'representative': 'GCF_000000001.1',
                    'organism': 'Testgenus namedii', 'taxid': '222', 'species_taxid': '222',
                    'genbank': '', 'source': 'official_fixture'}
        payloads = self.payloads()
        payloads[taxonomy.NCBI + 'GCF_000000001.1/dataset_report'] = {'reports': [{
            'accession': 'GCF_000000001.1', 'assembly_info': {'assembly_status': state},
            'organism': {'tax_id': '222', 'organism_name': 'Testgenus namedii'},
        }]}
        payloads[taxonomy.ENA + 'tax-id/222'] = node('222', 'Testgenus namedii', binomial='true')
        return evidence, metadata, payloads

    def test_successful_reference_mapping_precedes_generic_fallback(self):
        evidence, metadata, payloads = self.reference_fixture()
        with patch.object(taxonomy, 'index_metadata'), \
             patch.object(taxonomy.TaxonomyCache, 'source_for', return_value=metadata):
            result, report, _ = self.run_batch(
                frame=self.frame('Testgenus sp000000001'), evidence=evidence, payloads=payloads)
        self.assertEqual(result.loc[0, 'tax_id'], '222')
        self.assertEqual(report[0]['status'], 'REFERENCE_SUPPORTED_SPECIES_MAPPING')
        self.assertEqual(report[0]['fallback_attempts_json'], '')
        self.assertNotIn(query_url('Testgenus sp.'), self.urls)

    def test_suppressed_reference_can_use_generic_fallback_and_keeps_reason(self):
        evidence, metadata, payloads = self.reference_fixture(state='suppressed')
        with patch.object(taxonomy, 'index_metadata'), \
             patch.object(taxonomy.TaxonomyCache, 'source_for', return_value=metadata):
            result, report, summary = self.run_batch(
                frame=self.frame('Testgenus sp000000001'), evidence=evidence, payloads=payloads)
        self.assertEqual(summary['resolved'], 1)
        self.assertEqual(result.loc[0, 'tax_id'], '111')
        self.assertEqual(report[0]['primary_status'], 'ASSEMBLY_NOT_CURRENT')
        self.assertIn('suppressed', report[0]['primary_detail'])
        self.assertEqual(report[0]['reference_accession'], 'GCF_000000001.1')

    def test_failed_ani_uses_only_the_generic_candidate(self):
        evidence, metadata, payloads = self.reference_fixture(ani='90')
        with patch.object(taxonomy, 'index_metadata'), \
             patch.object(taxonomy.TaxonomyCache, 'source_for', return_value=metadata):
            result, report, _ = self.run_batch(
                frame=self.frame('Testgenus sp000000001'), evidence=evidence, payloads=payloads)
        self.assertEqual(result.loc[0, 'tax_id'], '111')
        self.assertEqual(report[0]['primary_status'], 'REFERENCE_BELOW_SPECIES_THRESHOLD')
        self.assertNotIn(taxonomy.ENA + 'tax-id/222', self.urls)

    def test_fallback_queries_are_auditable(self):
        _, report, _ = self.run_batch()
        attempts = json.loads(report[0]['fallback_attempts_json'])
        self.assertEqual(attempts[0]['query'], 'Testgenus sp.')
        self.assertEqual(attempts[0]['returned_tax_ids'], ['111'])

    def test_reference_metadata_failure_still_allows_environmental_lookup(self):
        evidence, _, payloads = self.reference_fixture()
        with patch.object(taxonomy, 'index_metadata', side_effect=ValueError('metadata unavailable')):
            result, report, summary = self.run_batch(
                frame=self.frame('Testgenus sp000000001'), evidence=evidence, payloads=payloads)
        self.assertEqual(summary['resolved'], 1)
        self.assertEqual(result.loc[0, 'tax_id'], '111')
        self.assertEqual(report[0]['primary_status'], 'REFERENCE_METADATA_UNAVAILABLE')
        self.assertIn('metadata unavailable', report[0]['primary_detail'])

    def test_reference_metadata_failure_does_not_change_an_unrelated_valid_pair(self):
        evidence, _, _ = self.reference_fixture()
        frame = pd.concat([self.frame('Testgenus sp000000001'),
                           self.frame('Othergenus namedii', '333').assign(sample_name='mag_2')],
                          ignore_index=True)
        record = node('333', 'Othergenus namedii', binomial='true')
        with patch.object(taxonomy, 'index_metadata', side_effect=OSError('download unavailable')):
            result, report, summary = self.run_batch(frame=frame, evidence=evidence,
                payloads={taxonomy.ENA + 'tax-id/333': record})
        self.assertEqual(len(summary['unresolved']), 1)
        self.assertEqual(report[1]['resolved'], 'yes')
        self.assertEqual(result.loc[1, 'tax_id'], '333')

    def test_service_failure_is_reported_with_the_primary_failure(self):
        _, report, summary = self.run_batch(payloads={
            query_url('Testgenus sp.'): taxonomy.Unresolved('LOOKUP_FAILED', 'HTTP 503'),
        })
        self.assertEqual(len(summary['unresolved']), 1)
        self.assertEqual(report[0]['primary_status'], 'NAME_ID_MISMATCH_OR_UNRESOLVED')
        self.assertEqual(report[0]['fallback_status'], 'LOOKUP_FAILED')

    def test_offline_fallback_reuses_the_http_cache(self):
        cache_path = self.root / 'offline.sqlite'
        cache = taxonomy.TaxonomyCache(cache_path)
        payloads = self.payloads()
        payloads[query_url('uncultured Testgenus bacterium')] = []
        for url, value in payloads.items():
            body = json.dumps(value)
            cache.db.execute('INSERT INTO http VALUES (?,?,?,?,?)',
                             (url, 200, body, hashlib.sha256(body.encode()).hexdigest(), 0))
        cache.db.commit()
        cache.db.close()
        result, _, summary = taxonomy.resolve_taxonomy_batch(
            self.frame(), evidence=self.evidence(), cache_path=cache_path, offline=True)
        self.assertEqual(result.loc[0, 'tax_id'], '111')
        self.assertEqual(summary['requests'], 0)

    def test_plain_prepare_command_does_not_require_resolve_taxonomy(self):
        metadata = self.root / 'input.tsv'
        self.frame().to_csv(metadata, sep='\t', index=False)
        with redirect_stdout(io.StringIO()):
            status = mag2ena.main(['--metadata', str(metadata), '--prepare-only'])
        self.assertEqual(status, 0)

    def test_cli_generic_default_and_explicit_opt_out(self):
        for extra, expected_status in (([], 0), (['--no-apply-generic-taxonomy'], 1),
                                       (['--apply-generic-taxonomy'], 0)):
            with self.subTest(extra=extra):
                metadata = self.root / 'input.tsv'
                self.frame().to_csv(metadata, sep='\t', index=False)
                payloads = self.payloads()
                with patch.object(taxonomy.TaxonomyCache, 'get', new=self.lookup(payloads)), \
                     redirect_stdout(io.StringIO()):
                    status = mag2ena.main([
                        '--metadata', str(metadata), '--resolve-taxonomy', '--prepare-only',
                        '--taxonomy-cache', str(self.root / 'cli.sqlite'),
                        '--taxonomy-report-out', str(self.root / 'report.tsv'), *extra,
                    ])
                self.assertEqual(status, expected_status)

    def test_unresolved_fallback_stops_cli_before_credentials(self):
        metadata = self.root / 'input.tsv'
        self.frame().to_csv(metadata, sep='\t', index=False)
        report_path = self.root / 'report.tsv'
        with patch.object(taxonomy.TaxonomyCache, 'get', return_value=None), \
             patch.object(mag2ena, 'resolve_webin_credentials') as credentials, \
             patch.object(mag2ena, 'build_and_submit') as submit, redirect_stdout(io.StringIO()):
            status = mag2ena.main([
                '--metadata', str(metadata), '--resolve-taxonomy',
                '--taxonomy-report-out', str(report_path),
                '--taxonomy-cache', str(self.root / 'cli.sqlite'),
            ])
        self.assertEqual(status, 1)
        credentials.assert_not_called()
        submit.assert_not_called()
        with taxonomy.unresolved_taxonomy_path(report_path).open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(rows[0]['fallback_status'], 'NO_ENVIRONMENTAL_MAPPING')
        self.assertTrue(rows[0]['suggested_action'])


if __name__ == '__main__':
    unittest.main()