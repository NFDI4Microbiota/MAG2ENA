<p align="center">
  <img src="App/imgs/mag2ena_logo.png" alt="MAG2ENA logo" width="300"/>
</p>

# MAG2ENA

**Version 1.0.7** — environment-aware curation, validation, submission, and
post-submission metadata maintenance for
**metagenome-assembled genomes (MAGs)** to the **European Nucleotide Archive
(ENA)** using the [GSC MIMAG / ERC000047](https://www.ebi.ac.uk/ena/browser/view/ERC000047)
checklist.

MAG2ENA connects parent-sample metadata, CheckM/CheckM2 quality estimates,
GTDB-Tk classifications, species-level taxonomy checks, curator-reviewed
environmental ontology resolution, FASTA preparation, ENA Webin submission,
and metadata-only updates to existing MAG samples. It is available as a hosted
[web application](https://web.app.ufz.de/mag2ena) and as a command-line tool
for reproducible submission pipelines.

Developed within [NFDI4Microbiota](https://nfdi4microbiota.de/), a German National Research Data Infrastructure (NFDI) consortium that serves as the central hub for the microbiology community, providing access to data, analysis services, metadata standards, and training — with a mission to make microbial research data FAIR (Findable, Accessible, Interoperable, and Reusable).

The tool provides two complementary interfaces:

| Interface | Best for |
|-----------|----------|
| **Hosted web app** | Guided submission and metadata updates, unified provenance, optional result encryption, and asynchronous job tracking; no installation is required |
| **CLI** (`mag2ena.py`) | Scriptable preparation and serial submission of new MAG collections without provenance recording |
| **Web source** (`App/`) | Included for groups that want to inspect, adapt, or self-host the web application |

Most users should use either the hosted web app or the CLI. The code under
`App/` does not need to be configured or run unless you are deploying your own
MAG2ENA web service. Self-hosting additionally requires Redis, an RQ worker,
persistent job storage, and Java; completion emails require an optional Mailgun
API key. Users of the hosted service never need to provide a Mailgun key.

---

# Quick Start

## Hosted web app

Open [web.app.ufz.de/mag2ena](https://web.app.ufz.de/mag2ena). No local
installation, Redis instance, worker, or Mailgun configuration is needed.

## CLI or self-hosted deployment

```bash
git clone https://github.com/NFDI4Microbiota/MAG2ENA.git
cd mag2ena
```

---

# Self-hosting the Web Application

The following options run the source in `App/` as your own web deployment. They
are not required for using the public hosted service. If you only need the CLI,
install the project with `uv sync` or `pip install -e .` and continue to the
CLI Tool section below; Redis and Mailgun are not required for CLI use.

## Option 1 — Docker (recommended for self-hosting)

Runs the web app with Redis and Java pre-configured, while persisting job history across container restarts:

```bash
docker build -t mag2ena .

docker run -p 8501:8501 \
  -e API_KEY="your-mailgun-api-key" \
  -v "$(pwd)/App/jobs:/app/App/jobs" \
  -v "$(pwd)/App/task-results:/app/App/task-results" \
  mag2ena
```

`API_KEY` is the optional Mailgun API key used by the web deployment to send
completion emails. Omit that line if email notification is not required; all
submission, update, job-tracking, and download features remain available.

Open `http://localhost:8501`.

## Option 2 — uv (local self-hosting)

```bash
# Install uv
curl -LsSf https://astral.sh/uv/install.sh | sh

# Install dependencies
uv sync

# Start Redis
redis-server --daemonize yes

# Start RQ worker (the bootstrap sets a safe working directory/import path)
uv run python App/worker.py &

# Start the web app
uv run streamlit run App/app.py
```

## Option 3 — pip (local self-hosting)

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .

redis-server --daemonize yes
python App/worker.py &
streamlit run App/app.py
```

## External requirements

Java is required for Webin-CLI. Install on Debian/Ubuntu:

```bash
sudo apt install default-jre
```

MAG2ENA uses Webin-CLI 11.0.0. Its JAR is bundled at:

```
App/webin-cli-11.0.0.jar
```

---

# Web Application

The hosted service is ready to use without local setup. The same five-page
Streamlit application is provided in `App/` for optional self-hosting.
Submission jobs run through Redis Queue, so processing continues after the
browser is closed. Assemblies within one job are submitted serially; separate
jobs run concurrently only when the deployment has multiple workers.

## Pages

### Home
Platform overview, evaluation highlights, and the complete three-step guided
submission workflow.

### Jobs
Track queued, running, successful, or failed submissions by Job ID; unlock
password-protected jobs and download their results.

### Update
Discover existing MAG samples from a study or a list of sample accessions, edit
sample metadata in batches of at most 1,000, validate the exact diff, and submit
an ENA metadata-only update without replacing assemblies. The same parent autofill,
taxonomy, environmental-resolution, fill-column, and quality-import assistance used
on Home is available during updates with unified provenance. This workflow is
available in the web application only.

### Help & Tutorials
Workflow guidance, safeguards, input specifications, FAQs, and tutorials.

### About
Project summary, authors, acknowledgements, data availability, and contact details.

---

## Submit Workflow

### Step 1 — Metadata

Enter metadata manually in the editable table, upload a TSV, or download the
provided template. Before FASTA upload, MAG2ENA reports checklist problems by row
and field.

**Metadata assistance tools:**

| Tool | Description |
|------|-------------|
| **ENA Autofill** | Enter a source sample accession (e.g. `ERS123456`) to auto-populate environmental fields. MAG2ENA uses the sample XML first and, if unavailable, recovers only metadata values that agree across assemblies exactly linked to that sample |
| **Taxonomy Resolver** | Look up NCBI Taxonomy IDs in bulk and reject non-species or metagenome taxonomies before submission |
| **Environmental Context Resolver** | Group shared parent evidence, resolve the three environmental fields with role-aware ontology constraints, review confidence and alternatives, and apply only curator-approved terms |
| **Fill Column** | Apply one value across an entire column, with enum-aware dropdowns for controlled vocabulary fields |
| **Import Quality Files** | Parse and merge CheckM, CheckM2, and GTDB-Tk output files directly into the metadata table |

**Supported quality file formats:**

| Tool | Expected filename |
|------|-------------------|
| CheckM v1 | `storage.tsv` |
| CheckM2 | `quality_report.tsv` |
| GTDB-Tk | `gtdbtk.bac120.summary.tsv` or `gtdbtk.ar53.summary.tsv` |

GTDB-Tk organism names are parsed with rank-aware logic: species names are preferred, falling back to genus/family with an `uncultured` prefix when only higher-rank assignments are available.

### Step 2 — FASTA Files

Upload one `.fasta.gz` file per MAG. The filename (without the `.fasta.gz` suffix) must exactly match the `sample_name` in the metadata table. One-sequence MAGs are detected automatically. When the submitter confirms that they represent complete chromosomes, MAG2ENA generates the ENA chromosome list using the exact FASTA sequence identifier and an ENA-safe chromosome name. An incomplete one-sequence MAG requires a single-contig exception from the ENA helpdesk and must not be declared as a complete chromosome.

### Step 3 — ENA Submission

| Field | Details |
|-------|---------|
| Study | Create a new study (name, title ≥ 20 chars, description ≥ 20 chars) or provide an existing accession (`PRJEBXXXX`) |
| Credentials | ENA Webin username and password |
| Portal | **Testing** (default — always validate here first) or **Production** |
| Email | Optional — receive a completion notification when the deployment has the Mailgun `API_KEY` environment variable configured |

An optional results password encrypts the completed job artifacts at rest. MAG2ENA
does not retain this password and cannot recover it. Clicking **Submit to ENA**
enqueues the job; save the resulting Job ID, and the results password when used,
to retrieve the outcome from **Jobs**.

---

## Update Workflow

Use **Update** to correct or enrich metadata for MAG samples already registered
in ENA. This changes the existing sample records; it does not create replacement
accessions or upload assemblies again.

1. Select **Testing** or **Production** and enter the Webin account that owns the records.
2. Discover MAG samples from a study accession or provide sample accessions directly.
3. Select a review batch of at most 1,000 samples and retrieve their authoritative ENA XML.
4. Edit the table manually or use parent autofill, taxonomy and environmental resolution, fill-column operations, or quality-file imports.
5. Review the exact field-level diff. MAG2ENA runs local ERC000047/taxonomy checks and a non-mutating ENA validation request.
6. Confirm the update. The worker retrieves the records again to detect concurrent changes before sending the metadata-only update.

Sample accessions, original aliases, study relationships, and assembly files are
not editable. Each completed update includes the original XML, submitted XML,
rollback XML, field-level diff, ENA receipts, and TSV/JSON provenance.

> Production updates alter existing ENA sample records. Use the correct portal
> and owning Webin account, review the diff carefully, and retain the rollback XML.

---

## Jobs Page

Look up any job by its ID to see:

- **Queued / Running**: Current position in the queue
- **Success**: Number of submitted MAGs, errors, assigned accessions, and downloadable artifacts
- **Failed**: Error details

**Downloadable outputs per job:**

| File | Contents |
|------|----------|
| `submit.xml` | ENA metadata XML sent to WEBIN v2 |
| `webin_log.xml` | ENA response containing sample accessions |
| `manifests.zip` | Per-MAG Webin-CLI manifest files |
| `provenance.tsv` / `provenance.json` | Unified human- and machine-readable web-workflow history |
| `success.txt` | Assemblies submitted successfully |
| `error.txt` | Assemblies that failed with error details |
| `result.json` | Final job summary |

Metadata-update jobs instead provide `before.xml`, `modify.xml`, `rollback.xml`,
`metadata_diff.tsv`, the ENA validation and update receipts, and the same unified
TSV/JSON provenance. Only samples with changed fields are sent to ENA.

---

# CLI Tool

`mag2ena.py` provides the web workflow's preparation safeguards and the same ENA
submission pathway without provenance recording or browser job tracking. Metadata
registration is automatically split into batches of up to 1,000 rows, and assembly
uploads remain serial to avoid the submission errors observed with parallel Webin-CLI
uploads.

The CLI creates new submissions but does not modify existing ENA sample records.
Post-submission metadata changes must use the web application's **Update** page,
which provides interactive diff review, ENA dry-run validation, concurrency checks,
rollback XML, and provenance.
## Installation with uv (Linux/macOS)

The CLI requires Git, Python 3.11 or newer, and Java 17 or newer. The repository pins Python 3.11 in `.python-version`; `uv` manages the Python environment and dependencies, while Java must be installed separately.

Install `uv` if it is not already available:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
source "$HOME/.local/bin/env"
uv --version
```
Clone the repository if you have not already done so, then enter its directory:

```bash
git clone https://github.com/NFDI4Microbiota/MAG2ENA.git
cd MAG2ENA
```

From the repository root, create the isolated Python environment using the committed dependency lockfile:

```bash
uv python install 3.11
uv sync --locked --python 3.11
```

This creates `.venv/` in the repository and installs the locked dependencies.

Check Java:

```bash
java -version
```

If Java is missing or older than 17, you will have to install it. On Debian/Ubuntu run:

```bash
sudo apt update
sudo apt install default-jre-headless
java -version
```

Verify MAG2ENA and the exact Webin-CLI JAR configured for its submission command:

```bash
uv run --locked python mag2ena.py --version
uv run --locked python mag2ena.py --help

webin_jar="$(uv run --locked python -c 'from mag2ena import WEBIN_JAR; print(WEBIN_JAR)')"
printf 'Configured Webin-CLI: %s\n' "$webin_jar"
java -jar "$webin_jar" -version
```

These checks do not submit data, but insted provide versions and confirmations on used packages. Downloading a different Webin-CLI JAR does not automatically switch MAG2ENA to it: the `WEBIN_JAR` path in `mag2ena.py` must refer to the intended file, so if you download a new version the call needs to be updated there. 

To verfy the MAG2ENA version, activate the environment first:

```bash
source .venv/bin/activate
python mag2ena.py --version
```

Which should return the current version of MAG2ENA and confirm the installation. You can now proceed with the CLI submission.

## Usage

```bash
python mag2ena.py --version

# Submit with a new study
python mag2ena.py \
  --metadata examples/metadata.tsv \
  --fasta-dir examples/fasta \
  --study-name 'My MAG study' \
  --study-title 'Descriptive title longer than 20 characters' \
  --study-description 'Study description longer than 20 characters'

# Submit to an existing study
python mag2ena.py \
  --metadata examples/metadata.tsv \
  --fasta-dir examples/fasta \
  --study-accession PRJEBXXXX

# Submit to production (default is test)
python mag2ena.py ... --portal prod

# Confirm complete chromosome-level one-sequence MAGs
python mag2ena.py ... --single-contig-as-chromosome
```

The CLI prompts for the Webin username and password after local validation. The
password is entered without terminal echo, and neither credential is accepted as a
command-line argument.

**Always validate using `--portal test` (the default) before submitting to production.**

## Optional metadata preparation

The CLI can reuse parent metadata, import CheckM/CheckM2 and GTDB-Tk output, resolve taxonomy from exact ENA names or original GTDB reference evidence, and create a four-column environmental review table. GTDB reference resolution requires the database release used for classification:

~~~bash
python mag2ena.py \
  --metadata metadata.tsv \
  --checkm quality_report.tsv \
  --gtdbtk gtdbtk.bac120.summary.tsv \
  --gtdb-release 214.1 \
  --autofill-parent-metadata \
  --resolve-taxonomy \
  --taxonomy-report-out taxonomy_report.tsv \
  --environment-review-out environment_review.tsv
  --prepared-metadata-out prepared_metadata.tsv \
  --prepare-only
~~~

Edit the decision column in the review TSV and mark approved rows `accept`. Apply
those decisions during a later preparation or submission:

~~~bash
python mag2ena.py \
  --metadata metadata.tsv \
  --autofill-parent-metadata \
  --environment-decisions environment_review.tsv \
  --prepared-metadata-out approved_metadata.tsv \
  --prepare-only
~~~

Autofill and imports preserve existing cells, and environmental proposals are never applied without explicit acceptance. Use `--environment-scope all` when reviewing existing environmental values as well as empty cells.

Repeat `--gtdbtk` to supply both bacterial and archaeal summaries. Original GTDB species names retain underscores during import; higher-rank classifications do not supply a species name. The taxonomy report preserves the complete original GTDB row, its classification, source paths and SHA256 hashes, and the reference accession, ANI, AF and selected metadata release.

Repeat `--gtdbtk` to supply both bacterial and archaeal summaries. Original GTDB species names retain underscores during import; higher-rank classifications do not supply a species name. The taxonomy report preserves every supplied GTDB column, its classification, source paths and SHA256 hashes, and any recorded reference accession, ANI, AF and selected metadata release.

A two-column `user_genome/classification` table supports classification import and exact ENA name/ID checks; it cannot provide reference-based mappings. Use the original full GTDB-Tk summary for placeholders and suffixed taxa, including cases where old metadata have already replaced underscores with spaces. When a species assignment has a recorded reference accession, the reference checks remain mandatory and failures do not fall back to a name lookup. When no reference is recorded, conventional species names can be checked against exact ENA names and synonyms; this verifies naming consistency rather than an ANI-based assignment. Existing verified name/ID pairs, including generic species-rank taxa, are not rejected merely because the GTDB classification stops above species. A new generic mapping still requires `--apply-generic-taxonomy`; a higher-rank GTDB classification cannot establish a new named descendant species.
Unresolved or conflicting rows remain in the report and cause exit code 1 before credentials are requested or submission starts, including in `--prepare-only` mode. Prepared metadata are still written for inspection. Without `--resolve-taxonomy`, submission preflight checks every supplied name/ID pair directly against ENA. Offline lookup is for preparation; submission preflight always checks current ENA records.

## Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `--metadata` | Yes | Path to metadata TSV |
| `--fasta-dir` | For submission | Directory containing `.fasta.gz` files |
| `--checkm` / `--gtdbtk` | No | Import quality estimates or original GTDB-Tk classifications; repeat `--gtdbtk` for both domains |
| `--autofill-parent-metadata` | No | Fill empty cells from ENA parent records with strict assembly fallback |
| `--resolve-taxonomy` | No | Resolve exact ENA names or release-matched GTDB references |
| `--gtdb-release` | For GTDB resolution | Original database point release, e.g. `214.1` |
| `--gtdb-metadata` | No | Local official metadata TSV/TSV.GZ; repeat for both domains |
| `--gtdb-min-af` | No | Original minimum AF; default 0.5 from R207, 0.65 earlier |
| `--taxonomy-scope` | No | Fill `empty` tax IDs (default), or permit supported corrections with `all` |
| `--taxonomy-report-out` | No | Detailed evidence TSV; defaults to `<metadata-stem>.taxonomy.tsv` |
| `--taxonomy-cache` | No | Reusable SQLite cache path |
| `--taxonomy-refresh-http` | No | Refresh cached live responses, retaining GTDB indexes |
| `--taxonomy-offline` | No | Prepare from cached evidence; submission still requires live checks |
| `--apply-generic-taxonomy` | No | Permit verified generic species-rank replacements |
| `--environment-review-out` | No | Write four-column environmental proposals for review |
| `--environment-decisions` | No | Apply proposals explicitly accepted in a review TSV |
| `--environment-scope` | No | Resolve `empty` fields (default) or normalize `all` environmental fields |
| `--prepared-metadata-out` | No | Save metadata after assistance and review decisions |
| `--prepare-only` | No | Prepare metadata without requesting credentials or submitting |
| `--non-interactive` | No | Read both credentials from environment variables and never prompt |
| `--portal` | No | `test` (default) or `prod` |
| `--single-contig-as-chromosome` | Conditional | Required when any FASTA has one sequence; confirms all such MAGs are complete chromosomes |
| `--study-accession` | No* | Existing study accession |
| `--study-name` | No* | New study short name |
| `--study-title` | No* | New study title (≥ 20 chars) |
| `--study-description` | No* | New study description (≥ 20 chars) |
| `--release-date` | No | Hold-until date in `YYYY-MM-DD`; must be within two years |
| `--output-dir` | No | New or empty results directory |
| `--version` | No | Print the MAG2ENA version and exit |

\* Either `--study-accession` or all three `--study-*` fields are required.

## Outputs

Each submission creates a timestamped results directory:

```
mag2ena_cli_results_<timestamp>/
├── metadata.submitted.tsv
├── submit.xml
├── webin_log.xml
├── manifests.zip
├── success.txt
├── error.txt             # Present only when errors occur
└── result.json
```

---

# Input Format

## Metadata TSV

Tab-separated file. One row per MAG.

**Mandatory fields:**

| Field | Notes |
|-------|-------|
| `sample_name` | Unique identifier — must match FASTA filename |
| `organism` | Scientific name (e.g. `uncultured Firmicutes bacterium`) |
| `tax_id` | NCBI Taxonomy ID |
| `genome coverage` | Positive number |
| `platform` | Sequencing platform (controlled vocabulary) |
| `assembly software` | Software name and version |
| `assembly quality` | Controlled vocabulary |
| `completeness score` | 0–100 |
| `contamination score` | 0–100 |
| `completeness software` | e.g. `CheckM2` |
| `binning software` | Software name |
| `binning parameters` | Free text |
| `metagenomic source` | Source metagenome accession |
| `isolation_source` | Free text |
| `collection date` | ISO 8601 format |
| `geographic location (country and/or sea)` | Controlled vocabulary |
| `geographic location (latitude)` | Decimal degrees |
| `geographic location (longitude)` | Decimal degrees |
| `broad-scale environmental context` | Biome term, preferably `label [ENVO:identifier]` |
| `local environmental context` | Smaller-grain ENVO or interoperable OBO term, preferably `label [identifier]` |
| `environmental medium` | Environmental-material term, preferably `label [ENVO:identifier]` |
| `taxonomic identity marker` | e.g. `multi-marker approach` |
| `sample derived from` | Source sample accession (e.g. `ERS123456`) |
| `project name` | Free text |

Any additional columns are submitted as ENA `SAMPLE_ATTRIBUTE` entries without modification.

A template TSV is available from the submission workflow on Home or at `template.tsv`.

The environmental resolver groups rows with the same parent sample and evidence, so one reviewed decision can be applied to many MAGs. It also reads explicit ENVO identifiers from parent-sample titles and descriptions, but reuses a term only when its label or surrounding language identifies the appropriate MIxS role (biome, local context, or environmental material). Semantic-role checks lower the confidence of material, spatial-feature, and process terms proposed for the wrong field. It never changes metadata during candidate generation. Users explicitly accept proposals before application and can inspect alternatives. In the web app, applications can be undone during the session and proposals, evidence, alternatives, confidence, ontology version, and curator decisions are recorded in the unified provenance ledger. The CLI instead exposes proposals and decisions through its four-column review TSV.

## Unified provenance ledger

The Home and Update workflows keep append-only histories of manual edits,
bioinformatics-file imports, parent autofill, taxonomy and environmental assistance,
curator decisions, validation, generated artifacts, and ENA outcomes. Submission
provenance additionally records FASTA checksums, transformations, and manifests;
Update provenance records original values, proposed changes, validation, and the
rollback artifact. Ledgers can be downloaded as `provenance.tsv` or
`provenance.json` and remain inside the encrypted archive when job encryption is
enabled. Credentials are never included.

## FASTA Files

- Format: gzip-compressed (`.fasta.gz`)
- One file per MAG
- Filename without extension must match `sample_name` exactly
- Maximum 1,000 MAGs per web job; the CLI accepts larger metadata tables in batches

---

# Credentials & Secrets

In the web app, Webin credentials are held only as long as needed to authenticate
the queued submission or metadata update. They are excluded from the job database,
provenance, downloadable artifacts, and logs, and are removed from the persisted
queue job when the worker starts processing it.

The CLI requests credentials only after local preflight. Password input uses a
non-echoing prompt, and neither credential has a command-line option, preventing
them from being stored in shell history or exposed in the Python process command line.
The Webin-CLI subprocess receives sensitive arguments through a temporary
mode-`0600` Java response file that is deleted immediately after use.

For automation, provide credentials through the environment or a secret manager:

~~~bash
export MAG2ENA_WEBIN_USER="Webin-XXXXX"
read -rsp "Webin password: " MAG2ENA_WEBIN_PASSWORD
export MAG2ENA_WEBIN_PASSWORD
echo
python mag2ena.py ... --non-interactive
unset MAG2ENA_WEBIN_USER MAG2ENA_WEBIN_PASSWORD
~~~

For operators self-hosting the web app, the optional Mailgun completion-email
feature reads its API key from the `API_KEY` environment variable. This is a
server-side deployment secret: do not ask MAG2ENA users to provide it, do not
commit it to the repository, and do not bake it into a container image. The
hosted MAG2ENA service and CLI users do not configure this value.

When using Docker, pass the environment variable at runtime to avoid baking credentials into the image:

```bash
docker run -p 8501:8501 \
  -e API_KEY="your-mailgun-api-key" \
  mag2ena
```

---

# Submission Workflow (Internal)

```
Metadata TSV
      │
      ▼
ERC000047 validation
      │   ├─ format and controlled-vocabulary checks
      │   └─ species-rank, non-metagenome taxonomy preflight
      │
      ▼
FASTA integrity, structure, header, and filename checks
      │
      ▼
ENA XML generation (SAMPLE_SET + optional PROJECT_SET)
      │
      ▼
WEBIN v2 metadata submission → sample accessions returned
      │
      ▼
Per-MAG manifest generation
      │   └─ confirmed complete one-sequence chromosome? → chromosome list generated
      ▼
Serial Webin-CLI 11.0.0 assembly submission
      │
      ▼
Results: accessions, receipts, manifests, and success/error logs
      │
      └─ Web only: provenance TSV/JSON and optional encrypted archive
```

The web metadata-update path is separate:

```
Study or sample accession → authenticated ENA sample XML → assisted editing
      → exact diff → local validation → ENA dry run → explicit approval
      → concurrent-change check → metadata-only MODIFY → receipt and rollback XML
```

---

# Common Errors

**Missing FASTA file**
```
Missing FASTA files for: sample_X
```
The FASTA filename (without `.fasta.gz`) does not match `sample_name` in the metadata table.

**Validation errors**
Check for: wrong date format (use ISO 8601), coordinates with wrong notation, enum values not matching the checklist, or empty mandatory fields.

**Taxonomy rejected**
MAG2ENA requires every `tax_id` to resolve through ENA as exactly species rank
and not represent a metagenome. Use the species-level parent for subspecies or
strain assignments; genus and higher-rank IDs remain blocked.

**Existing sample cannot be loaded or updated**
Select the portal where the sample was originally registered and use the Webin
account that owns it. A public record may be readable without being editable by
the current account.

**Existing sample changed after review**
MAG2ENA cancels the update if ENA metadata changed after the batch was loaded.
Reload the record, review the new diff, and validate it again.

**Incomplete single-contig MAG**
Do not select the complete-chromosome confirmation. ENA requires a single-contig
exception for an incomplete one-sequence assembly; contact `datasubs@ebi.ac.uk`
before submitting it.

**Webin-CLI errors**
See `error.txt` in the job output. Common causes: invalid manifest field values, temporary ENA API issues, or missing metadata values that passed local validation but failed server-side.

---

# References

- ENA Submission Portal: https://www.ebi.ac.uk/ena/browser/submit
- ENA Checklist ERC000047: https://www.ebi.ac.uk/ena/browser/view/ERC000047
- GSC MIxS Standard: https://www.nature.com/articles/nbt1366
- MIxS Term Browser: https://w3id.org/mixs/
- ENVO Ontology: https://www.ebi.ac.uk/ols4/ontologies/envo

---

# Acknowledgements

This work was supported by the Deutsche Forschungsgemeinschaft (DFG, German
Research Foundation) through NFDI4Microbiota — National Research Data
Infrastructure for Microbiota Research (project 460129525); the Canadian
International Development Research Centre (IDRC, grant 109981); the UK Foreign,
Commonwealth and Development Office; the São Paulo Research Foundation (FAPESP,
grant 2024/10958-1); the Google PhD Fellowship; and the Initiative and Networking
Fund of the Helmholtz Association through the Helmholtz Metadata Collaboration
project call (project ZT-I-PF-3-108).

# License and Contributions

MAG2ENA is released under the [MIT License](LICENSE). Issues and contributions
are welcome on [GitHub](https://github.com/NFDI4Microbiota/MAG2ENA).

# Authors

- **Breno L. S. de Almeida** — corresponding author ([brenoslivio@usp.br](mailto:brenoslivio@usp.br))
- **Anderson P. Avila-Santos** — corresponding author ([anderson.santos@ufz.de](mailto:anderson.santos@ufz.de))
- Martin Bole
- Sanchita Kamath
- Robson P. Bonidia
- Peter F. Stadler
- André C. P. L. F. de Carvalho
- **Ulisses N. da Rocha** — corresponding author ([ulisses.rocha@ufz.de](mailto:ulisses.rocha@ufz.de))
