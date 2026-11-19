"""Resolve MAG taxonomy using ENA records and release-matched GTDB evidence."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import re
import sqlite3
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from contextlib import contextmanager
from pathlib import Path

import pandas as pd

ENA = 'https://www.ebi.ac.uk/ena/taxonomy/rest/'
ENA_XML = 'https://www.ebi.ac.uk/ena/browser/api/xml/'
NCBI = 'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/'
GTDB = 'https://data.gtdb.ecogenomic.org/releases/'
USER_AGENT = 'MAG2ENA-taxonomy'
FIELDS = ['metadata_file', 'sample_name', 'original_organism', 'original_tax_id',
          'status', 'detail', 'candidate_organism', 'candidate_tax_id',
          'proposed_organism', 'proposed_tax_id', 'apply_allowed',
          'current_ena_name', 'current_rank', 'current_submittable', 'current_binomial',
          'gtdb_release', 'gtdb_classification', 'assigned_gtdb_species',
          'reference_accession', 'reference_evidence', 'ANI', 'ANI_radius', 'AF',
          'reference_matches_release', 'release_representative',
          'release_ncbi_organism', 'release_ncbi_tax_id', 'release_ncbi_species_tax_id',
          'reference_taxonomy_source', 'reference_organism', 'reference_tax_id',
          'proposed_rank', 'proposed_submittable', 'proposed_binomial',
          'mapping_scope', 'release_species_clusters_sharing_ncbi_id',
          'lookup_utc', 'input_metadata_sha256', 'gtdb_sources_json', 'gtdb_record_json',
          'gtdb_metadata_sha256', 'ena_record_url', 'candidate_verified',
          'resolved', 'applied', 'blocking_reason',
          'primary_status', 'primary_detail', 'fallback_status', 'fallback_detail',
          'fallback_source_rank', 'fallback_source_name', 'fallback_query',
          'fallback_attempts_json', 'name_normalization']


class Unresolved(ValueError):
    def __init__(self, status, detail=''):
        super().__init__(detail or status)
        self.status, self.detail = status, detail


def utc():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def truth(value):
    return str(value).lower() in ('true', 't', '1')


def clean(value):
    value = str(value or '').strip()
    return '' if value.casefold() in ('n/a', 'na', 'none', 'nan', 'null') else value


def norm(value):
    return ' '.join(str(value or '').strip().casefold().split())


def accession(value):
    value = clean(value)
    return re.sub(r'^(RS_|GB_)', '', value)


def species(taxonomy):
    return next((t.strip()[3:] for t in str(taxonomy).split(';')
                 if t.strip().startswith('s__')), '')


def domain(taxonomy):
    return next((t.strip()[3:] for t in str(taxonomy).split(';')
                 if t.strip().startswith('d__')), '')



def sha_file(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


@contextmanager
def text_input(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as handle:
        yield handle



class TaxonomyCache:
    def __init__(self, path, offline=False, refresh=False):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.db = sqlite3.connect(str(path))
        self.offline, self.requests, self.hits, self.last = offline, 0, 0, 0
        self.db.executescript('''
        CREATE TABLE IF NOT EXISTS http(url TEXT PRIMARY KEY, code INTEGER,
          body TEXT, sha TEXT, fetched REAL);
        CREATE TABLE IF NOT EXISTS gtdb(release TEXT, accession TEXT, genbank TEXT,
          taxonomy TEXT, species TEXT, representative TEXT, taxid TEXT,
          species_taxid TEXT, organism TEXT, source TEXT,
          PRIMARY KEY(release, accession));
        CREATE INDEX IF NOT EXISTS gtdb_genbank ON gtdb(release, genbank);
        CREATE INDEX IF NOT EXISTS gtdb_species_id ON gtdb(release, species_taxid);
        CREATE TABLE IF NOT EXISTS sources(release TEXT, source TEXT, sha TEXT,
          rows INTEGER, fetched TEXT, PRIMARY KEY(release,source));
        ''')
        if refresh:
            self.db.execute('DELETE FROM http WHERE url LIKE ?', (ENA + '%',))
            self.db.execute('DELETE FROM http WHERE url LIKE ?', (ENA_XML + '%',))
            self.db.execute('DELETE FROM http WHERE url LIKE ?', (NCBI + '%',))
        self.db.commit()

    def get(self, url, kind='json'):
        cached = self.db.execute('SELECT code,body,sha,fetched FROM http WHERE url=?', (url,)).fetchone()
        if cached and (self.offline or time.time() - cached[3] < 7 * 86400):
            code, body, digest, _ = cached
            if hashlib.sha256(body.encode()).hexdigest() != digest:
                raise ValueError('HTTP cache checksum mismatch: {}'.format(url))
            self.hits += 1
        else:
            if self.offline:
                raise Unresolved('NOT_CACHED', url)
            for attempt in range(3):
                time.sleep(max(0, 0.4 - (time.monotonic() - self.last)))
                self.last = time.monotonic()
                self.requests += 1
                try:
                    req = urllib.request.Request(url, headers={'User-Agent': USER_AGENT,
                        'Accept': 'application/json' if kind == 'json' else 'application/xml'})
                    with urllib.request.urlopen(req, timeout=30) as response:
                        data = response.read(16 * 1024 * 1024 + 1)
                    if len(data) > 16 * 1024 * 1024:
                        raise ValueError('Unexpectedly large API response: {}'.format(url))
                    body, code = data.decode('utf-8'), 200
                    json.loads(body) if kind == 'json' else ET.fromstring(body)
                    break
                except urllib.error.HTTPError as exc:
                    if exc.code in (400, 404):
                        code, body = exc.code, ''
                        break
                    if attempt == 2 or exc.code not in (408, 429, 500, 502, 503, 504):
                        raise Unresolved('LOOKUP_FAILED', '{}: HTTP {}'.format(url, exc.code))
                    time.sleep(2 ** attempt)
                except (urllib.error.URLError, TimeoutError, OSError) as exc:
                    if attempt == 2:
                        raise Unresolved('LOOKUP_FAILED', '{}: {}'.format(url, exc))
                    time.sleep(2 ** attempt)
            self.db.execute('INSERT OR REPLACE INTO http VALUES (?,?,?,?,?)',
                (url, code, body, hashlib.sha256(body.encode()).hexdigest(), time.time()))
            self.db.commit()
        if code != 200:
            return None
        return json.loads(body) if kind == 'json' else body

    def ena(self, taxid):
        if not re.fullmatch(r'[1-9]\d*', str(taxid)):
            raise Unresolved('INVALID_TAX_ID', str(taxid))
        record = self.get(ENA + 'tax-id/' + str(taxid))
        if not isinstance(record, dict) or not record.get('taxId'):
            raise Unresolved('TAX_ID_NOT_RETURNED', str(taxid))
        return record

    def species_node(self, taxid):
        record = self.ena(taxid)
        if record.get('rank') == 'species':
            return record
        text = self.get(ENA_XML + str(record['taxId']), 'xml')
        if not text:
            raise Unresolved('NO_SPECIES_ANCESTOR', str(taxid))
        root = ET.fromstring(text)
        candidates = {n.get('taxId') for n in root.findall('.//lineage/taxon')
                      if n.get('rank') == 'species'}
        if len(candidates) != 1:
            raise Unresolved('NO_UNIQUE_SPECIES_ANCESTOR', str(taxid))
        return self.ena(candidates.pop())

    def source_for(self, release, reference):
        self.db.row_factory = sqlite3.Row
        rows = self.db.execute('SELECT * FROM gtdb WHERE release=? AND accession=?',
                              (release, reference)).fetchall()
        if not rows:
            rows = self.db.execute('SELECT * FROM gtdb WHERE release=? AND genbank=?',
                                  (release, reference)).fetchall()
        self.db.row_factory = None
        if not rows:
            raise Unresolved('REFERENCE_ABSENT_FROM_SELECTED_RELEASE', reference)
        rows = [dict(r) for r in rows]
        if len({r['taxonomy'] for r in rows}) != 1:
            raise Unresolved('REFERENCE_ALIAS_CONFLICT', reference)
        return rows[0]


class DigestReader:
    def __init__(self, raw):
        self.raw, self.digest = raw, hashlib.sha256()

    def read(self, size=-1):
        data = self.raw.read(size)
        self.digest.update(data)
        return data


def index_metadata(cache, release, source):
    source = str(source)
    marker = re.search(r'_r(\d{2,3})(?:\.|_)', source.rsplit('/', 1)[-1])
    if marker and marker[1] != release.split('.')[0]:
        raise ValueError('Metadata filename is R{}, but --release is {}'.format(marker[1], release))
    done = cache.db.execute('SELECT sha,rows FROM sources WHERE release=? AND source=?',
                            (release, source)).fetchone()
    if done:
        if not source.startswith('https://') and sha_file(Path(source)) != done[0]:
            raise ValueError('Cached GTDB metadata changed: use a new cache or restore ' + source)
        print('Using cached GTDB {} index: {} genomes ({})'.format(release, done[1], source), flush=True)
        return
    remote = source.startswith('https://')
    if remote and cache.offline:
        raise ValueError('GTDB metadata not indexed; run online once or supply --gtdb-metadata')
    print('Indexing GTDB {} metadata: {}'.format(release, source), flush=True)
    raw = urllib.request.urlopen(urllib.request.Request(source, headers={'User-Agent': USER_AGENT}),
                                 timeout=120) if remote else open(source, 'rb')
    digest = DigestReader(raw)
    stream = gzip.GzipFile(fileobj=digest) if source.endswith('.gz') else digest
    text = io.TextIOWrapper(stream, encoding='utf-8-sig', newline='') if source.endswith('.gz') else None
    if text is None:
        raw.close()
        # Remote metadata are streamed without retaining an extracted file.
        if remote:
            raise ValueError('Remote GTDB metadata must be gzip compressed')
        text = open(source, 'r', encoding='utf-8-sig', newline='')
    count, chunk = 0, []
    try:
        reader = csv.DictReader(text, delimiter='\t')
        required = {'accession', 'gtdb_taxonomy', 'gtdb_genome_representative',
                    'ncbi_taxid', 'ncbi_species_taxid', 'ncbi_organism_name'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError('Not an official GTDB metadata table: {}'.format(source))
        cache.db.execute('DELETE FROM gtdb WHERE release=? AND source=?', (release, source))
        cache.db.commit()
        for row in reader:
            if None in row or None in row.values():
                raise ValueError('Malformed GTDB metadata row in {}'.format(source))
            tax = row['gtdb_taxonomy']
            chunk.append((release, accession(row['accession']),
                          accession(row.get('ncbi_genbank_assembly_accession', '')), tax,
                          species(tax), accession(row['gtdb_genome_representative']),
                          clean(row['ncbi_taxid']), clean(row['ncbi_species_taxid']),
                          clean(row['ncbi_organism_name']), source))
            count += 1
            if len(chunk) == 2000:
                cache.db.executemany('INSERT OR REPLACE INTO gtdb VALUES (?,?,?,?,?,?,?,?,?,?)', chunk)
                cache.db.commit()
                chunk = []
            if count % 100000 == 0:
                print('  Indexed {} genomes'.format(count), flush=True)
        if chunk:
            cache.db.executemany('INSERT OR REPLACE INTO gtdb VALUES (?,?,?,?,?,?,?,?,?,?)', chunk)
        digest_value = digest.digest.hexdigest() if source.endswith('.gz') else sha_file(Path(source))
        cache.db.execute('INSERT OR REPLACE INTO sources VALUES (?,?,?,?,?)',
                         (release, source, digest_value, count, utc()))
        cache.db.commit()
        print('  Complete: {} genomes; SHA256 {}'.format(count, digest_value), flush=True)
    finally:
        text.close()
        raw.close()


def supported_reference(row, min_af):
    assigned = species(row.get('classification', ''))
    if not assigned:
        raise Unresolved('NO_SPECIES_ASSIGNMENT', 'Saved GTDB-Tk classification stops above species')
    reference_fields = ('fastani_reference', 'closest_genome_reference', 'closest_placement_reference')
    if not any(field in row for field in reference_fields):
        raise Unresolved(
            'GTDB_REFERENCE_COLUMNS_MISSING',
            'The supplied GTDB table has no reference columns. Supply the original '
            'GTDB-Tk summary with reference accession, taxonomy, radius, ANI and AF; '
            'a user_genome/classification import table cannot establish this mapping',
        )
    warnings = norm(row.get('warnings', '') + ' ' + row.get('note', ''))
    if 'incongruent' in warnings or 'not assigned to closest species' in warnings:
        raise Unresolved('SPECIES_REFERENCE_CONFLICT', 'Saved GTDB-Tk output reports conflicting/failed species evidence')
    for prefix in ('fastani_', 'closest_genome_', 'closest_placement_'):
        reference = accession(row.get(prefix + 'reference', ''))
        if not reference:
            continue
        if not re.fullmatch(r'GC[AF]_\d{9}\.\d+', reference):
            raise Unresolved('INVALID_REFERENCE_ACCESSION', reference)
        reference_species = species(row.get(prefix + 'taxonomy', ''))
        if reference_species != assigned:
            raise Unresolved('REFERENCE_SPECIES_DISAGREEMENT', '{} vs {}'.format(reference_species, assigned))
        try:
            ani = float(row[prefix + 'ani'])
            radius = float(row[prefix + 'reference_radius'] if prefix in ('fastani_', 'closest_genome_')
                           else row[prefix + 'radius'])
            af = float(row[prefix + 'af'])
        except (KeyError, ValueError):
            raise Unresolved('REFERENCE_EVIDENCE_INCOMPLETE', reference)
        if not all(math.isfinite(v) for v in (ani, radius, af)) or not 0 < radius <= 100 or not 0 <= ani <= 100 or not 0 <= af <= 1:
            raise Unresolved('INVALID_REFERENCE_METRICS', reference)
        if ani < radius or af < min_af:
            raise Unresolved('REFERENCE_BELOW_SPECIES_THRESHOLD', 'ANI {} / radius {}; AF {} / minimum {}'.format(ani, radius, af, min_af))
        return reference, prefix.rstrip('_'), ani, radius, af
    raise Unresolved('NO_RECORDED_SPECIES_REFERENCE', assigned)


def reference_taxonomy(cache, reference, metadata):
    report = cache.get(NCBI + urllib.parse.quote(reference, safe='') + '/dataset_report')
    reports = report.get('reports', []) if isinstance(report, dict) else []
    matches = [r.get('assembly', r) for r in reports if isinstance(r, dict)]
    matches = [r for r in matches if r.get('accession') == reference]
    if len(matches) > 1:
        raise Unresolved('AMBIGUOUS_ASSEMBLY_REPORT', reference)
    if matches:
        r = matches[0]
        if r.get('current_accession', reference) != reference:
            raise Unresolved('ASSEMBLY_VERSION_CHANGED', reference)
        info = r.get('assembly_info', {})
        state = str(info.get('assembly_status', '')).lower()
        if any(x in state for x in ('suppressed', 'withdrawn', 'replaced')):
            raise Unresolved('ASSEMBLY_NOT_CURRENT', reference + ': ' + state)
        organism = r.get('organism', {})
        taxid = str(organism.get('tax_id', ''))
        if re.fullmatch(r'[1-9]\d*', taxid):
            return taxid, organism.get('organism_name', ''), 'NCBI_ASSEMBLY', True
    # GCF/GCA prefixes are not interchangeable; use the recorded GenBank pair.
    ena_reference = reference if reference.startswith('GCA_') else metadata['genbank']
    if ena_reference:
        text = cache.get(ENA_XML + ena_reference, 'xml')
        if text:
            assemblies = [r for r in ET.fromstring(text).findall('.//ASSEMBLY')
                          if r.get('accession') == ena_reference]
            if len(assemblies) == 1:
                taxid = clean(assemblies[0].findtext('.//TAXON/TAXON_ID'))
                if re.fullmatch(r'[1-9]\d*', taxid):
                    return taxid, assemblies[0].findtext('.//TAXON/SCIENTIFIC_NAME', ''), 'ENA_ASSEMBLY', True
    # Snapshot-only mappings are reported but not applied.
    taxid = metadata['taxid']
    if re.fullmatch(r'[1-9]\d*', taxid):
        return taxid, metadata['organism'], 'GTDB_RELEASE_SNAPSHOT_ONLY', False
    raise Unresolved('REFERENCE_TAXON_UNAVAILABLE', reference)


def validate_node(node, expected_domain):
    if node.get('rank') != 'species':
        raise Unresolved('TAXON_NOT_SPECIES_RANK', str(node.get('rank')))
    if not truth(node.get('submittable')):
        raise Unresolved('TAXON_NOT_SUBMITTABLE', str(node.get('taxId')))
    if (truth(node.get('metagenome')) or 'metagenome' in norm(node.get('scientificName'))
            or 'metagenomes;' in norm(node.get('lineage'))):
        raise Unresolved('METAGENOME_TAXON_NOT_MAG_ORGANISM', node.get('scientificName', ''))
    if expected_domain and not re.search(r'(?<![A-Za-z])' + re.escape(expected_domain) + r'(?![A-Za-z])',
                                         str(node.get('lineage', '')), re.I):
        raise Unresolved('DOMAIN_MISMATCH', expected_domain + ' vs ' + str(node.get('lineage', '')))


def generic(node):
    name = norm(node.get('scientificName'))
    return not truth(node.get('binomial')) or bool(re.search(r'(?:\bsp\.|\bbacterium|\bcyanobacterium|\barchaeon|\bprokaryote)$', name))


def names(node):
    result = {norm(node.get('scientificName'))}
    classes = {'synonym', 'equivalent name', 'includes', 'misspelling', 'genbank synonym',
               'scientific name', 'homotypic synonym', 'heterotypic synonym', 'anamorph', 'teleomorph'}
    for value in node.get('otherNames', []):
        if isinstance(value, str):
            name, _, kind = value.rpartition(':')
            if kind in classes:
                result.add(norm(name))
        elif isinstance(value, dict):
            if value.get('nameClass', value.get('name_class')) in classes:
                result.add(norm(value.get('name')))
    return result - {''}


def set_proposal(out, node, status, scope, allowed):
    out.update(status=status, candidate_organism=node['scientificName'],
               candidate_tax_id=str(node['taxId']),
               proposed_organism=node['scientificName'] if allowed else '',
               proposed_tax_id=str(node['taxId']) if allowed else '', proposed_rank=node.get('rank', ''),
               proposed_submittable=str(node.get('submittable', '')).lower(),
               proposed_binomial=str(node.get('binomial', '')).lower(),
               mapping_scope=scope, apply_allowed='yes' if allowed else 'no',
               candidate_verified='yes',
               ena_record_url=ENA + 'tax-id/' + str(node['taxId']))


def resolve_named(row, saved, out, cache, apply_generic):
    current = None
    try:
        current = cache.ena(row['tax_id'].strip())
    except Unresolved as exc:
        if exc.status not in {'INVALID_TAX_ID', 'TAX_ID_NOT_RETURNED'}:
            raise
    if current:
        out.update(current_ena_name=current.get('scientificName', ''), current_rank=current.get('rank', ''),
                   current_submittable=str(current.get('submittable', '')).lower(),
                   current_binomial=str(current.get('binomial', '')).lower())
    expected = domain(saved.get('classification', ''))
    candidates = []
    if current and norm(row['organism']) in names(current):
        candidates = [current]
    else:
        result = cache.get(ENA + 'any-name/' + urllib.parse.quote(row['organism'].strip(), safe=''))
        for node in result if isinstance(result, list) else []:
            if norm(row['organism']) in names(node):
                candidates.append(node)
    usable = {}
    reasons = []
    for node in candidates:
        try:
            original_rank = node.get('rank')
            if original_rank in ('strain', 'subspecies', 'no rank'):
                node = cache.species_node(str(node['taxId']))
            validate_node(node, expected)
            usable[str(node['taxId'])] = (node, original_rank)
        except Unresolved as exc:
            reasons.append(exc.status)
    if len(usable) != 1:
        raise Unresolved('NAME_ID_MISMATCH_OR_UNRESOLVED' if not usable else 'AMBIGUOUS_EXACT_NAME',
                         '; '.join(reasons) or 'No unique exact scientific-name/synonym match; no fuzzy name replacement made')
    node, rank = next(iter(usable.values()))
    is_generic = generic(node)
    status = 'GENERIC_NAME_ID_MATCH' if is_generic else 'EXACT_NAME_OR_SYNONYM_MATCH'
    if rank in ('strain', 'subspecies', 'no rank'):
        status = 'MATCHED_NAME_SPECIES_ANCESTOR'
    set_proposal(out, node, status, 'GENERIC_ORGANISM' if is_generic else 'NAME_ID_CONSISTENCY_ONLY',
                 not is_generic or apply_generic)
    out['detail'] = 'Checks name/ID consistency, not independent genomic species identification'
    return node


def resolve_gtdb_reference(row, saved, out, cache, release, min_af, apply_generic):
    reference, evidence, ani, radius, af = supported_reference(saved, min_af)
    out.update(reference_accession=reference, reference_evidence=evidence, ANI=ani, ANI_radius=radius, AF=af)
    metadata = cache.source_for(release, reference)
    source = cache.db.execute('SELECT sha FROM sources WHERE release=? AND source=?',
                              (release, metadata['source'])).fetchone()
    out.update(release_representative=metadata['representative'], release_ncbi_organism=metadata['organism'],
               release_ncbi_tax_id=metadata['taxid'], release_ncbi_species_tax_id=metadata['species_taxid'],
               gtdb_metadata_sha256=source[0] if source else '')
    if metadata['species'] != species(saved.get('classification', '')):
        out['reference_matches_release'] = 'no'
        raise Unresolved('GTDB_RELEASE_SPECIES_MISMATCH',
                         'Saved species {} vs R{} {}'.format(species(saved.get('classification', '')), release, metadata['species']))
    out['reference_matches_release'] = 'yes'
    taxid, organism, source_kind, live = reference_taxonomy(cache, reference, metadata)
    out.update(reference_taxonomy_source=source_kind, reference_organism=organism, reference_tax_id=taxid)
    node = cache.species_node(taxid)
    validate_node(node, domain(saved.get('classification', '')))
    is_generic = generic(node)
    # Cluster overlap describes the selected release snapshot.
    count = 0
    if metadata['species_taxid']:
        count = cache.db.execute('SELECT COUNT(DISTINCT species) FROM gtdb WHERE release=? AND species_taxid=?',
                                 (release, metadata['species_taxid'])).fetchone()[0]
    out['release_species_clusters_sharing_ncbi_id'] = count
    status = 'GENERIC_REFERENCE_MAPPING' if is_generic else 'REFERENCE_SUPPORTED_SPECIES_MAPPING'
    if not live:
        status = 'SNAPSHOT_ONLY_REFERENCE_MAPPING'
    set_proposal(out, node, status, 'GENERIC_ORGANISM' if is_generic else
                 ('NCBI_CONCEPT_SPANS_MULTIPLE_GTDB_CLUSTERS' if count > 1 else 'MATCHED_REFERENCE_SPECIES'),
                 live and (not is_generic or apply_generic))
    if not live:
        out['candidate_verified'] = 'no'
        out['detail'] = 'Only the pinned GTDB snapshot supplied the assembly tax ID; not auto-applied'
    elif is_generic:
        out['detail'] = 'Reference has a generic/non-binomial taxon; no other cluster member was selected by name preference'
    elif count > 1:
        out['detail'] = 'Actual matched-reference taxon; NCBI species concept is broader than one GTDB cluster in the release snapshot'



def read_gtdb_summaries(paths):
    """Read original GTDB-Tk rows and retain their complete source evidence."""
    evidence = {}
    for value in paths or []:
        path = Path(value).resolve()
        source = {'path': str(path), 'sha256': sha_file(path)}
        with text_input(path) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            fields = reader.fieldnames or []
            if len(fields) != len(set(fields)) or not {'user_genome', 'classification'} <= set(fields):
                raise ValueError(f'GTDB-Tk summary needs unique user_genome/classification columns: {path}')
            for row in reader:
                if None in row or None in row.values():
                    raise ValueError(f'Malformed GTDB-Tk row in {path}')
                key = row['user_genome'].strip()
                if not key:
                    raise ValueError(f'Empty GTDB-Tk user_genome in {path}')
                if key in evidence and evidence[key]['row'] != row:
                    raise ValueError(f'Conflicting GTDB-Tk records for {key}')
                item = evidence.setdefault(key, {'row': dict(row), 'sources': []})
                if source not in item['sources']:
                    item['sources'].append(dict(source))
    return evidence


def gtdb_organism_table(evidence):
    """Import species assignments verbatim; higher ranks remain unassigned."""
    return pd.DataFrame(
        [{'sample_name': name, 'organism': species(item['row']['classification'])}
         for name, item in evidence.items()],
        columns=['sample_name', 'organism'],
    )


def _cell(value):
    return '' if pd.isna(value) else str(value)


def _gtdb_label(value):
    return bool(re.search(r'\bsp\d+\b|[A-Za-z0-9]_[A-Z]+(?:\b|_)', str(value)))


def _has_reference_value(saved):
    return any(clean(saved.get(field, '')) for field in (
        'fastani_reference', 'closest_genome_reference', 'closest_placement_reference',
    ))


def _deepest_gtdb_name(saved):
    parts = str(saved.get('classification', '')).split(';')
    for part in reversed(parts):
        _, separator, name = part.strip().partition('__')
        if separator and name.strip():
            return name.strip()
    return ''


def _resolve_primary_taxonomy_row(row, saved, out, cache, release, min_af, apply_generic):
    assigned = species(saved.get('classification', ''))
    needs_reference = (
        _gtdb_label(row['organism']) or _gtdb_label(_deepest_gtdb_name(saved))
    )
    if (assigned and _has_reference_value(saved)) or needs_reference:
        if not saved:
            raise Unresolved('GTDB_EVIDENCE_REQUIRED',
                             'Supply the original GTDB-Tk summary with --gtdbtk')
        resolve_gtdb_reference(row, saved, out, cache, release, min_af, apply_generic)
        return

    proposal = dict(out)
    node = resolve_named(row, saved, proposal, cache, apply_generic)
    if assigned and norm(assigned) not in names(node):
        raise Unresolved(
            'METADATA_GTDB_SPECIES_CONFLICT',
            f'Saved GTDB species {assigned!r} is not an exact name/synonym of '
            f"the proposed ENA taxon {node['scientificName']!r}; original reference evidence is needed",
        )
    if saved and not assigned and not generic(node):
        supplied_pair = False
        if clean(row['tax_id']):
            try:
                supplied = cache.ena(row['tax_id'].strip())
                supplied_pair = norm(row['organism']) in names(supplied)
            except Unresolved as exc:
                if exc.status not in {'INVALID_TAX_ID', 'TAX_ID_NOT_RETURNED'}:
                    raise
        if not supplied_pair:
            raise Unresolved(
                'NO_SPECIES_ASSIGNMENT',
                'Saved GTDB-Tk classification stops above species; an exact name '
                'lookup cannot establish a new named species assignment for this MAG',
            )
    out.update(proposal)


def _without_gtdb_suffixes(value):
    return re.sub(r'(?<=[A-Za-z0-9])_[A-Z]+(?=\s|$)', '', str(value)).strip()


def _fallback_levels(row, saved):
    ranks = {}
    for part in str(saved.get('classification', '')).split(';'):
        rank, separator, name = part.strip().partition('__')
        if separator and name.strip():
            ranks[rank] = name.strip()
    levels = [(rank, ranks[rank]) for rank in ('g', 'f') if rank in ranks]
    if not levels:
        name = re.sub(r'^uncultured\s+', '', row['organism'].strip(), flags=re.I)
        name = _without_gtdb_suffixes(name).split()
        if name and re.fullmatch(r'[A-Za-z][A-Za-z0-9_-]*', name[0]):
            label = name[0]
            if label.casefold() not in {'bacterium', 'archaeon', 'prokaryote',
                                         'unclassified', 'unidentified', 'unknown'}:
                levels.append(('f' if label.endswith('aceae') else 'g', label))
    return levels


def _fallback_query(cache, out, endpoint, query):
    attempts = json.loads(out.get('fallback_attempts_json') or '[]')
    attempt = {'endpoint': endpoint, 'query': query}
    attempts.append(attempt)
    try:
        records = cache.get(ENA + endpoint + '/' + urllib.parse.quote(query, safe=''))
        records = records if isinstance(records, list) else []
        attempt['returned_tax_ids'] = [str(node.get('taxId', '')) for node in records]
        return records
    except Unresolved as exc:
        attempt.update(status=exc.status, detail=exc.detail)
        if exc.status != 'NOT_CACHED':
            raise
        return []
    finally:
        out['fallback_attempts_json'] = json.dumps(attempts, ensure_ascii=False)


def _fallback_node(cache, candidate, expected_domain):
    node = cache.ena(str(candidate.get('taxId', '')))
    if str(node.get('taxId', '')) != str(candidate.get('taxId', '')):
        raise Unresolved('TAX_ID_CHANGED', str(candidate.get('taxId', '')))
    validate_node(node, expected_domain)
    return node


def _resolve_supplied_suffix_pair(row, saved, out, cache):
    original = row['organism'].strip()
    assigned = species(saved.get('classification', ''))
    source = assigned if assigned and original == assigned.replace('_', ' ') else original
    query = _without_gtdb_suffixes(source)
    if query == source or not re.fullmatch(r'[1-9]\d*', row['tax_id'].strip()):
        return False
    node = cache.ena(row['tax_id'].strip())
    validate_node(node, domain(saved.get('classification', '')))
    if norm(query) not in names(node):
        return False
    if assigned and norm(_without_gtdb_suffixes(assigned)) not in names(node):
        return False
    set_proposal(out, node, 'CANONICAL_NAME_FOR_SUPPLIED_ID',
                 'SUPPLIED_ID_NAME_CONSISTENCY', True)
    out.update(name_normalization=f'{source} -> {query}',
               fallback_status='CANONICAL_NAME_FOR_SUPPLIED_ID',
               fallback_query=query,
               detail='Verified the supplied tax ID against the suffix-normalized name; '
                      'the original GTDB name remains in the evidence report')
    return True


def resolve_environmental_fallback(row, saved, out, cache, apply_generic):
    """Find a submission taxon using a genus or family label."""
    expected_domain = domain(saved.get('classification', ''))
    descriptors = ('archaeon',) if expected_domain == 'Archaea' else ('bacterium', 'cyanobacterium')
    for rank, original in _fallback_levels(row, saved):
        label = _without_gtdb_suffixes(original)
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_-]*', label):
            continue
        queries = [f'{label} sp.'] if rank == 'g' else []
        queries += [f'{label} {descriptor}' for descriptor in descriptors]
        queries += [f'uncultured {label} {descriptor}' for descriptor in descriptors]
        queries.append(f'uncultured {label} sp.')
        for query in queries:
            candidates = _fallback_query(cache, out, 'any-name', query)
            usable = {}
            for candidate in candidates:
                if norm(query) not in names(candidate):
                    continue
                try:
                    node = _fallback_node(cache, candidate, expected_domain)
                    if norm(query) in names(node) and generic(node):
                        usable[str(node['taxId'])] = node
                except Unresolved as exc:
                    if exc.status in {'LOOKUP_FAILED', 'NOT_CACHED'}:
                        raise
            if len(usable) != 1:
                continue
            node = next(iter(usable.values()))
            status = 'GENUS_GENERIC_MAPPING' if rank == 'g' else 'FAMILY_GENERIC_MAPPING'
            scope = 'GENUS_LEVEL_ENVIRONMENTAL_TAXON' if rank == 'g' else 'FAMILY_LEVEL_ENVIRONMENTAL_TAXON'
            set_proposal(out, node, status, scope, apply_generic)
            out.update(fallback_status=status, fallback_source_rank='genus' if rank == 'g' else 'family',
                       fallback_source_name=original, fallback_query=query,
                       name_normalization=f'{original} -> {label}' if original != label else '',
                       detail='Environmental submission taxon selected from the supported lineage label; '
                              'no genomic species match is asserted')
            return True

        if not any(character.isdigit() for character in label):
            continue
        candidates = _fallback_query(cache, out, 'suggest-for-search', label)
        usable = {}
        pattern = re.compile(r'(?:\bsp\.|\bbacterium|\barchaeon|\bcyanobacterium)\s+'
                             + re.escape(label) + r'$', re.I)
        for candidate in candidates:
            if not pattern.search(candidate.get('scientificName', '')):
                continue
            try:
                node = _fallback_node(cache, candidate, expected_domain)
                if pattern.search(node.get('scientificName', '')):
                    usable[str(node['taxId'])] = node
            except Unresolved as exc:
                if exc.status in {'LOOKUP_FAILED', 'NOT_CACHED'}:
                    raise
        if len(usable) != 1:
            continue
        node = next(iter(usable.values()))
        set_proposal(out, node, 'IDENTIFIER_NAME_MAPPING', 'IDENTIFIER_NAME_MATCH_ONLY', apply_generic)
        out.update(fallback_status='IDENTIFIER_NAME_MAPPING',
                   fallback_source_rank='genus' if rank == 'g' else 'family',
                   fallback_source_name=original, fallback_query=label,
                   detail='Exact placeholder identifier found in the ENA organism name; '
                          'this is a name-based submission mapping, without ANI species confirmation')
        return True
    out.update(fallback_status='NO_ENVIRONMENTAL_MAPPING',
               fallback_detail='No unique submittable species-rank genus, family or exact-identifier candidate')
    return False


def _resolve_taxonomy_row(row, saved, out, cache, release, min_af, apply_generic,
                          reference_error=None):
    failure = None
    try:
        if reference_error and species(saved.get('classification', '')) and _has_reference_value(saved):
            raise reference_error
        _resolve_primary_taxonomy_row(row, saved, out, cache, release, min_af, apply_generic)
        if out['candidate_verified'] == 'yes':
            out['primary_status'] = out['status']
            out['primary_detail'] = out['detail']
            return
        failure = Unresolved(out['status'], out['detail'])
    except Unresolved as exc:
        failure = exc
    out.update(primary_status=failure.status, primary_detail=failure.detail)
    try:
        if _resolve_supplied_suffix_pair(row, saved, out, cache):
            return
        if resolve_environmental_fallback(row, saved, out, cache, apply_generic):
            return
    except Unresolved as exc:
        out.update(fallback_status=exc.status, fallback_detail=exc.detail)
        raise
    raise failure


def _validated_release(release):
    value = str(release or '').lstrip('rR')
    if not re.fullmatch(r'\d{2,3}\.\d+', value):
        raise ValueError('GTDB reference resolution requires --gtdb-release, e.g. 214.1')
    return value


def _prepare_reference_metadata(cache, evidence, names_to_resolve, release, metadata_paths):
    domains = {
        domain(evidence[name]['row']['classification'])
        for name in names_to_resolve if name in evidence
        and species(evidence[name]['row']['classification'])
        and _has_reference_value(evidence[name]['row'])
    }
    if not domains:
        return
    release = _validated_release(release)
    if metadata_paths:
        for path in metadata_paths:
            index_metadata(cache, release, Path(path).resolve())
    else:
        major = release.split('.')[0]
        archaeal_prefix = 'ar53' if int(major) >= 207 else 'ar122'
        for name, prefix in [('Bacteria', 'bac120'), ('Archaea', archaeal_prefix)]:
            if name in domains:
                url = f'{GTDB}release{major}/{release}/{prefix}_metadata_r{major}.tsv.gz'
                index_metadata(cache, release, url)


def resolve_taxonomy_batch(
    frame,
    *,
    evidence=None,
    release=None,
    cache_path=Path('taxonomy_lookup_cache.sqlite'),
    metadata_paths=(),
    min_af=None,
    scope='empty',
    apply_generic=True,
    refresh_http=False,
    offline=False,
    input_metadata_path='',
    input_metadata_sha256='',
):
    """Resolve each row and return prepared metadata, evidence and a summary."""
    if not {'sample_name', 'organism', 'tax_id'} <= set(frame.columns):
        raise ValueError('Metadata require sample_name, organism and tax_id columns')
    if scope not in {'empty', 'all'}:
        raise ValueError('Taxonomy scope must be empty or all')
    keys = [_cell(value).strip() for value in frame['sample_name']]
    if any(not key for key in keys) or len(keys) != len(set(keys)):
        raise ValueError('Metadata sample_name values must be nonempty and unique')
    evidence = evidence or {}
    if any(name in evidence and species(evidence[name]['row']['classification']) for name in keys):
        release = _validated_release(release)
    elif release:
        release = _validated_release(release)
    if min_af is None:
        min_af = 0.65 if release and int(release.split('.')[0]) < 207 else 0.5
    if not math.isfinite(min_af) or not 0 <= min_af <= 1:
        raise ValueError('GTDB minimum AF must be between 0 and 1')
    if offline and refresh_http:
        raise ValueError('Offline taxonomy lookup cannot refresh HTTP responses')

    result = frame.copy(deep=True)
    reports = []
    cache = TaxonomyCache(Path(cache_path), offline=offline, refresh=refresh_http)
    try:
        reference_error = None
        try:
            _prepare_reference_metadata(cache, evidence, keys, release, metadata_paths)
        except (Unresolved, ValueError, OSError, ET.ParseError) as exc:
            reference_error = Unresolved(
                'REFERENCE_METADATA_UNAVAILABLE', f'{type(exc).__name__}: {exc}'
            )
        consecutive_failures = 0
        for position, (_, original) in enumerate(frame.iterrows()):
            row = {column: _cell(original[column]) for column in frame.columns}
            key = keys[position]
            item = evidence.get(key, {})
            saved = item.get('row', {})
            out = dict.fromkeys(FIELDS, '')
            out.update(
                metadata_file=str(input_metadata_path),
                sample_name=row['sample_name'], original_organism=row['organism'],
                original_tax_id=row['tax_id'], gtdb_release=release or '',
                gtdb_classification=saved.get('classification', ''),
                assigned_gtdb_species=species(saved.get('classification', '')),
                gtdb_sources_json=json.dumps(item.get('sources', []), ensure_ascii=False),
                gtdb_record_json=json.dumps(saved, ensure_ascii=False),
                input_metadata_sha256=input_metadata_sha256, lookup_utc=utc(),
                apply_allowed='no', candidate_verified='no', resolved='no', applied='no',
            )
            try:
                if consecutive_failures >= 5:
                    raise Unresolved('SERVICE_UNAVAILABLE',
                                     'Five consecutive lookup failures; rerun when the services recover')
                _resolve_taxonomy_row(
                    row, saved, out, cache, release, min_af, apply_generic,
                    reference_error=reference_error,
                )

                current = (row['organism'], row['tax_id'])
                candidate = (out['candidate_organism'], out['candidate_tax_id'])
                if out['candidate_verified'] != 'yes':
                    out['blocking_reason'] = out['detail'] or 'Current reference taxonomy could not be verified'
                elif current == candidate:
                    out['resolved'] = 'yes'
                elif out['apply_allowed'] != 'yes':
                    out['blocking_reason'] = ('Generic mappings are disabled by --no-apply-generic-taxonomy; '
                                              'omit that option or use --apply-generic-taxonomy')
                elif (scope == 'empty' and clean(row['tax_id']) and clean(row['tax_id']) != candidate[1]):
                    out['blocking_reason'] = 'Populated taxonomy needs a correction; use --taxonomy-scope all'
                else:
                    result.iat[position, result.columns.get_loc('organism')] = candidate[0]
                    result.iat[position, result.columns.get_loc('tax_id')] = candidate[1]
                    out.update(resolved='yes', applied='yes')
            except Unresolved as exc:
                out.update(status=exc.status, detail=exc.detail, blocking_reason=exc.detail)
            except (ValueError, ET.ParseError, OSError, KeyError) as exc:
                out.update(status='LOOKUP_FAILED', detail=f'{type(exc).__name__}: {exc}')
                out['blocking_reason'] = out['detail']
            consecutive_failures = (
                consecutive_failures + 1 if out['status'] == 'LOOKUP_FAILED'
                else consecutive_failures if out['status'] == 'SERVICE_UNAVAILABLE' else 0
            )
            reports.append(out)
        summary = {
            'resolved': sum(row['resolved'] == 'yes' for row in reports),
            'applied': sum(row['applied'] == 'yes' for row in reports),
            'unresolved': [
                {'sample_name': row['sample_name'], 'status': row['status'],
                 'reason': row['blocking_reason'] or row['detail']}
                for row in reports if row['resolved'] != 'yes'
            ],
            'status_counts': dict(Counter(row['status'] for row in reports)),
            'requests': cache.requests,
            'cache_hits': cache.hits,
        }
        return result, reports, summary
    finally:
        cache.db.close()


def write_taxonomy_report(path, reports, *, protected_paths=()):
    path = Path(path)
    if path.resolve() in {Path(value).resolve() for value in protected_paths}:
        raise ValueError(f'Taxonomy report would overwrite an input or cache: {path}')
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(reports)


def unresolved_taxonomy_path(report_path):
    """Place the unresolved review beside the complete evidence report."""
    path = Path(report_path)
    return path.with_name(path.stem + '.unresolved.tsv')


def taxonomy_suggested_action(row):
    """Explain the next step without assigning or searching for a taxon."""
    status = str(row.get('status', ''))
    reason = str(row.get('blocking_reason', '') or row.get('detail', ''))
    if '--apply-generic-taxonomy' in reason:
        return (
            'Generic mappings are enabled by default. Review the candidate and remove '
            '--no-apply-generic-taxonomy, or use --apply-generic-taxonomy explicitly.'
        )
    if '--taxonomy-scope all' in reason:
        return 'Review the proposed correction; use --taxonomy-scope all only if that correction is intended.'
    if status in {'ASSEMBLY_NOT_CURRENT', 'ASSEMBLY_VERSION_CHANGED'}:
        return (
            'Review the recorded reference assembly status and suppression/replacement reason. '
            'Reassess the MAG taxonomy using other genomic evidence before assigning a replacement.'
        )
    if status in {'GTDB_REFERENCE_COLUMNS_MISSING', 'GTDB_EVIDENCE_REQUIRED',
                  'REFERENCE_EVIDENCE_INCOMPLETE', 'MAG_ABSENT_FROM_SAVED_GTDB_RESULTS'}:
        return (
            'Supply the original full GTDB-Tk summary with matching sample IDs and '
            'reference accession, taxonomy, radius, ANI and AF. Retain the original evidence.'
        )
    if status in {'GTDB_RELEASE_SPECIES_MISMATCH', 'REFERENCE_ABSENT_FROM_SELECTED_RELEASE',
                  'REFERENCE_SPECIES_DISAGREEMENT', 'REFERENCE_ALIAS_CONFLICT',
                  'METADATA_GTDB_SPECIES_CONFLICT'}:
        return (
            'Verify the original GTDB release, sample match, reference accession/version '
            'and conflicting names. Resolve the discrepancy before applying a mapping.'
        )
    if status in {'REFERENCE_BELOW_SPECIES_THRESHOLD', 'SPECIES_REFERENCE_CONFLICT',
                  'INVALID_REFERENCE_METRICS', 'NO_RECORDED_SPECIES_REFERENCE'}:
        return (
            'Review the original reference and ANI/AF evidence and thresholds. If species '
            'support is absent, obtain further identification evidence or curate an appropriate '
            'species-rank environmental taxon; do not infer a named species from proximity alone.'
        )
    if status == 'SNAPSHOT_ONLY_REFERENCE_MAPPING':
        return 'Verify the current taxonomy for the recorded reference/version; the historical snapshot alone was not accepted.'
    if status in {'LOOKUP_FAILED', 'SERVICE_UNAVAILABLE', 'NOT_CACHED',
                  'NOT_LOOKED_UP_AFTER_SERVICE_FAILURE', 'REFERENCE_METADATA_UNAVAILABLE'}:
        return 'Retry online with the same cache when the taxonomy services are available; offline mode needs previously cached responses.'
    return (
        'Review a submittable species-rank environmental taxon consistent with the supported '
        'lineage and verify its ENA name/ID. If no existing taxon covers the organism, request '
        'one through ENA Webin Register taxonomy. Record the evidence for the curated choice.'
    )


def write_unresolved_taxonomy_report(report_path, reports, *, protected_paths=()):
    """Write every unresolved row, including rows with nonempty tax IDs."""
    path = unresolved_taxonomy_path(report_path)
    protected = {Path(value).resolve() for value in (*protected_paths, report_path)}
    if path.resolve() in protected:
        raise ValueError(f'Unresolved taxonomy review would overwrite an input/output/cache: {path}')
    fields = [
        'sample_name', 'original_organism', 'original_tax_id', 'gtdb_classification',
        'status', 'blocking_reason', 'reference_accession', 'candidate_organism',
        'candidate_tax_id', 'primary_status', 'primary_detail', 'fallback_status',
        'fallback_detail', 'fallback_source_rank', 'fallback_source_name',
        'fallback_query', 'name_normalization', 'suggested_action',
    ]
    rows = []
    for report in reports:
        if report.get('resolved') == 'yes':
            continue
        row = {field: report.get(field, '') for field in fields}
        row['blocking_reason'] = report.get('blocking_reason') or report.get('detail', '')
        row['suggested_action'] = taxonomy_suggested_action(report)
        rows.append(row)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    return path


def taxonomy_failure_message(summary, report_path, *, unresolved_path=None, prepared_path=None):
    failures = summary['unresolved']
    lines = [f'Taxonomy unresolved for {len(failures)} MAG(s); submission stopped.']
    review = unresolved_path or unresolved_taxonomy_path(report_path)
    lines.append(f'Unresolved taxonomy review: {Path(review).resolve()}')
    lines.append(f'Complete taxonomy evidence: {Path(report_path).resolve()}')
    if prepared_path:
        lines.append(f'Prepared metadata for correction: {Path(prepared_path).resolve()}')
    else:
        lines.append('Use --prepared-metadata-out to save the prepared table for correction.')
    for item in failures[:10]:
        lines.append(f"  {item['sample_name']}: {item['status']}: {item['reason']}")
    if len(failures) > 10:
        lines.append(f'  {len(failures) - 10} additional failures are listed in the unresolved review.')
    lines.append(
        'Review suggested_action for every unresolved row, record supported corrections '
        'in a copy of the metadata and validate that curated file before submission.'
    )
    return '\n'.join(lines)


def validate_submission_taxonomy_records(records, *, fetcher=None):
    """Check final species IDs and their corresponding names before submission."""
    errors, fetched = [], {}
    if fetcher is None:
        cache = TaxonomyCache(Path(':memory:'))
        fetcher = cache.ena
    else:
        cache = None
    try:
        for number, row in enumerate(records, start=1):
            taxid = _cell(row.get('tax_id', '')).strip()
            organism = _cell(row.get('organism', '')).strip()
            try:
                if not re.fullmatch(r'[1-9]\d*', taxid):
                    raise Unresolved('INVALID_TAX_ID', 'A numeric species-level tax ID is required')
                if taxid not in fetched:
                    fetched[taxid] = fetcher(taxid)
                node = fetched[taxid]
                if not isinstance(node, dict) or node.get('_error') or node.get('_not_found'):
                    raise Unresolved('TAXON_LOOKUP_FAILED', str(node))
                if str(node.get('taxId', '')) != taxid:
                    raise Unresolved('TAX_ID_CHANGED', 'Resolve the current tax ID before submission')
                validate_node(node, '')
                if norm(organism) not in names(node):
                    raise Unresolved('ORGANISM_TAX_ID_MISMATCH',
                                     f"{organism!r} does not match {node.get('scientificName', '')!r}")
            except (Unresolved, ValueError, OSError, ET.ParseError) as exc:
                reason = f'{exc.status}: {exc.detail}' if isinstance(exc, Unresolved) else str(exc)
                errors.append({
                    'row': number, 'sample_name': _cell(row.get('sample_name', '')),
                    'tax_id': taxid, 'expected': reason,
                })
        return errors
    finally:
        if cache:
            cache.db.close()