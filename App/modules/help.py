# modules/help.py
from pathlib import Path

import streamlit as st
from utils.css_injection import inject_css


_VIDEO_DIR = Path(__file__).resolve().parent.parent / "videos"


def tutorials():
    st.markdown("### Video Tutorials")
    st.markdown(
        "Here you will find practical use-case video tutorials that guide you through "
        "MAG2ENA and demonstrate its submission and metadata-update workflows:"
    )

    with st.expander("**Use case 1: Submitting MAGs through MAG2ENA**"):
        _l, col, _r = st.columns([1, 3, 1])
        with col:
            video = _VIDEO_DIR / "Tutorial1.mp4"
            if video.exists():
                st.video(str(video))
            else:
                st.warning("Tutorial1.mp4 is not available in this deployment.")

    with st.expander("**Use case 2: Updating metadata for MAGs in ENA**"):
        _l, col, _r = st.columns([1, 3, 1])
        with col:
            video = _VIDEO_DIR / "Tutorial2.mp4"
            if video.exists():
                st.video(str(video))
            else:
                st.warning("Tutorial2.mp4 is not available in this deployment.")


def usage_guide():
    st.markdown("### From assembled MAGs to ENA accessions")
    st.markdown(
        "MAG2ENA is designed for the point where genome reconstruction is finished but "
        "submission-ready metadata is not. It brings parent records, quality estimates, "
        "taxonomy, environmental ontology curation, validation, deposition, and subsequent "
        "metadata updates into one traceable workflow."
    )

    col1, col2 = st.columns(2, gap="large")
    with col1:
        st.markdown("#### Who is it for?")
        st.markdown(
            """
            - Researchers depositing one study without learning ENA XML or invoking Webin-CLI directly.
            - Bioinformaticians connecting CheckM/CheckM2 and GTDB-Tk outputs to a submission pipeline.
            - Data stewards reviewing consistent metadata across MAGs derived from shared samples.
            - Submitters correcting or enriching metadata for MAGs already registered in ENA.
            - Consortia processing large MAG collections through a web interface or CLI.

            You should already have gzip-compressed MAG assemblies and the information needed
            to describe how they were assembled and binned.
            """
        )
    with col2:
        st.markdown("#### What do you get?")
        st.markdown(
            """
            - Preflighted MIMAG metadata with row- and field-level error reporting.
            - Registered ENA study and sample accessions, plus one assembly submission per MAG.
            - Generated submission XML and per-MAG Webin manifests.
            - ENA receipts, Webin-CLI success/error logs, and an overall result summary.
            - Metadata-only updates with an exact diff and rollback XML for existing MAG samples.
            - A unified provenance ledger in TSV and JSON for curation review and reproducibility.
            """
        )

    st.markdown("### Guided submission workflow")
    st.markdown(
        """
        1. **Load metadata.** Try the bundled two-MAG example, download the ERC000047-shaped template, or upload an existing TSV.
        2. **Reuse parent evidence.** ENA autofill retrieves collection date, coordinates, country, isolation source, metagenomic source, and environmental context without overwriting populated cells.
        3. **Connect bioinformatics outputs.** Import CheckM/CheckM2 quality reports and GTDB-Tk summaries instead of copying results manually.
        4. **Resolve identifiers.** Search ENA taxonomy and review field-aware ENVO proposals for the three environmental context fields.
        5. **Run preflight.** Fix checklist, format, vocabulary, taxonomy, and FASTA problems before transmission.
        6. **Submit asynchronously.** Register metadata through ENA Webin REST and deposit assemblies through the bundled Webin-CLI 9.0.3.
        7. **Retrieve and audit.** Open the Job ID to download accessions, receipts, XML, manifests, logs, and final provenance.
        """
    )

    st.markdown("### Update metadata already in ENA")
    st.markdown(
        "The **Update** page changes sample metadata only; it never uploads or replaces "
        "assemblies. Load MAG samples from a study accession or paste sample accessions "
        "directly. Study results are divided into review batches of at most 1,000. "
        "MAG2ENA retrieves each authoritative sample record with your Webin account, "
        "preserves fields outside the editor, and submits only samples with an explicit "
        "field-level change. Parent autofill, taxonomy resolution, the environmental "
        "resolver, fill-column operations, and CheckM/CheckM2 or GTDB-Tk imports remain "
        "available for curator-approved edits."
    )
    st.markdown(
        "Before an update can be queued, the current diff must pass local checklist and "
        "taxonomy checks and an ENA dry run. The worker then retrieves the records again "
        "and cancels the whole batch if any changed since review. Results include the "
        "before XML, MODIFY XML, rollback XML, diff, receipt, and unified TSV/JSON provenance."
    )

    st.markdown("### Metadata assistance with curator control")
    st.markdown(
        """
        | Capability | What MAG2ENA does | What remains under user control |
        |---|---|---|
        | **ENA parent autofill** | Reads the parent sample XML; if unavailable, keeps only values that agree across exactly linked assembly records | Existing cells are preserved and conflicting fallback values are omitted |
        | **Quality import** | Maps CheckM/CheckM2 completeness, contamination, and software fields by MAG name | The table remains editable before validation |
        | **GTDB-Tk import** | Extracts the lowest available taxonomic assignment | ENA/NCBI taxon selection is separately checked |
        | **Taxonomy resolver** | Suggests ENA taxonomy records and fills selected IDs | Preflight accepts only species-rank, non-metagenome taxa |
        | **Environmental resolver** | Groups shared evidence, searches field-appropriate ontology space, ranks candidates, and explains confidence | No proposal changes metadata until the curator approves it |
        """
    )

    st.markdown("### Input files")
    st.markdown(
        """
        - **Metadata:** tab-separated values (`.tsv`). Mandatory fields are marked with `*` in the editor.
        - **Assemblies:** one gzip-compressed FASTA (`.fasta.gz`) per MAG. Its filename without the suffix must exactly match `sample_name`.
        - **CheckM v1:** `storage.tsv` with `Bin Id`, `Completeness`, and `Contamination`.
        - **CheckM2:** `quality_report.tsv` with `Name`, `Completeness`, and `Contamination`.
        - **GTDB-Tk:** bacterial or archaeal summary TSV with `user_genome` and `classification`.
        """
    )

    st.markdown("### Submission safeguards")
    st.markdown(
        """
        - Required fields, regular expressions, and enumerated vocabularies are derived from ERC000047 and reported by row and field.
        - Every taxon ID must resolve through the ENA taxonomy service at exactly species rank. Higher ranks, subspecies/strain records, and metagenome taxonomies are blocked.
        - FASTA gzip integrity, record structure, and sample-name matching are checked; duplicate sequence identifiers are detected and normalized before assembly submission.
        - Complete one-sequence chromosomes require explicit confirmation and receive an automatically generated chromosome list.
        - Incomplete one-sequence MAGs are blocked because ENA requires a single-contig exception from its helpdesk.
        - Testing and Production are separate choices; use Testing first for a disposable end-to-end check.
        """
    )
    st.caption(
        "Checklist validation uses the ERC000047 XML bundled with the application. ENA remains "
        "the final authority and may apply additional server-side rules."
    )

    st.markdown("### Provenance and results protection")
    st.markdown(
        "The provenance ledger records sources, before/after values, imports, resolver evidence "
        "and decisions, validation outcomes, file hashes, generated artifacts, and ENA responses. "
        "An optional results password encrypts completed job artifacts at rest. MAG2ENA does not "
        "store that password and cannot recover it."
    )

    st.markdown("### Command-line workflow")
    st.markdown(
        "MAG2ENA includes a CLI for reproducible pipelines and collections larger than "
        "the 1,000-MAG web limit. It supports parent-record autofill, CheckM/CheckM2 and "
        "GTDB-Tk imports, species-level taxonomy resolution, and a four-column environmental "
        "review TSV. The CLI does not record provenance or create browser Jobs. Metadata "
        "registration is batched automatically, while assembly uploads remain serial."
    )
    st.markdown("**Prepare metadata without submitting**")
    st.code(
        """python mag2ena.py \\
  --metadata metadata.tsv \\
  --checkm quality_report.tsv \\
  --gtdbtk gtdbtk.bac120.summary.tsv \\
  --autofill-parent-metadata \\
  --resolve-taxonomy \\
  --environment-review-out environment_review.tsv \\
  --prepared-metadata-out prepared_metadata.tsv \\
  --prepare-only""",
        language="bash",
    )
    st.markdown("**Submit to an existing study**")
    st.code(
        """python mag2ena.py \\
  --metadata prepared_metadata.tsv \\
  --fasta-dir fasta/ \\
  --study-accession PRJEBXXXX""",
        language="bash",
    )
    st.caption(
        "Submission defaults to Webin Testing. After preflight, interactive runs prompt for "
        "the Webin username and a non-echoed password; neither credential has a command-line "
        "option. For automation, use MAG2ENA_WEBIN_USER and MAG2ENA_WEBIN_PASSWORD with "
        "--non-interactive and obtain the secret from a secure environment or secret manager."
    )


def faq():
    st.markdown("### Frequently Asked Questions")

    with st.expander("Which ENA checklist does this tool use?"):
        st.write(
            "ERC000047 — the MIMAG package from the Genomic Standards Consortium. It defines the "
            "minimum information required for a Metagenome-Assembled Genome submission."
        )
    with st.expander("Why does MAG2ENA require a species-level taxonomy ID?"):
        st.write(
            "MAG2ENA applies a conservative MAG submission policy: each tax_id must resolve to "
            "exactly species rank in the ENA taxonomy service and must not represent a metagenome. "
            "Genus or higher-rank IDs and subspecies/strain IDs are blocked. For an assignment below "
            "species rank, use its species-level parent taxon. This is checked again by the background "
            "worker immediately before submission."
        )
    with st.expander("What does a 'group' mean in the Environmental resolver?"):
        st.write(
            "A review group is a set of MAG rows that share the same parent sample and contextual "
            "evidence for one target field. Broad-scale context, local context, and environmental "
            "medium are separate review units. One curator decision can therefore be applied "
            "consistently to related MAGs without copying it row by row."
        )
    with st.expander("Does autofill or the Environmental resolver overwrite my metadata?"):
        st.write(
            "ENA autofill fills only empty cells. Environmental candidates are shown as proposals "
            "with evidence, alternatives, and a confidence tier; they are applied only after you "
            "select and approve them. Resolver applications are reversible during the session and "
            "recorded in provenance."
        )
    with st.expander("Do I need a pre-existing ENA study (project)?"):
        st.write(
            "No. You can create a new study as part of the submission by providing a study name, title, "
            "and description. If you already have a study accession (e.g. PRJEB12345), select "
            "'Existing study accession' and paste it in."
        )
    with st.expander("What is the Testing portal and when should I use it?"):
        st.write(
            "ENA's testing environment exercises metadata registration and assembly upload without creating "
            "lasting public records; test objects are discarded by ENA. Run a complete test first to catch "
            "account, metadata, manifest, and FASTA issues before selecting Production."
        )
    with st.expander("My MAG is a single contig — do I need to create a chromosome list?"):
        st.write(
            "MAG2ENA detects one-sequence FASTA files automatically. If you confirm that the "
            "sequence represents a complete chromosome, MAG2ENA generates the required chromosome "
            "list before passing it to Webin-CLI. Do not confirm an incomplete one-sequence MAG: "
            "ENA requires these submissions to request a single-contig exception through "
            "datasubs@ebi.ac.uk."
        )
    with st.expander("What happens if some MAGs succeed and others fail?"):
        st.write(
            "Each MAG is submitted independently via Webin-CLI. Successes and errors are logged separately "
            "in the job directory. Check the Jobs tab for a breakdown and the full Webin output."
        )
    with st.expander("Are asynchronous jobs also run in parallel?"):
        st.write(
            "Asynchronous means the submission continues after you leave the browser and can be reopened "
            "with its Job ID. Each worker handles one job at a time, and assemblies within a job are submitted "
            "serially. Separate jobs run simultaneously only when the deployed service has multiple workers; "
            "otherwise they wait safely in the queue."
        )
    with st.expander("Can I resubmit a failed MAG without re-running everything?"):
        st.write(
            "Currently each submission session is independent. For partial re-submissions, filter your "
            "metadata TSV to only the failed samples and start a new submission."
        )
    with st.expander("Can I update metadata for MAGs that are already in ENA?"):
        st.write(
            "Yes. Open Update and discover MAG samples from a study accession or enter sample "
            "accessions directly. MAG2ENA edits the existing sample objects through ENA MODIFY; "
            "it does not resubmit assemblies or create replacement sample accessions. You must "
            "use a Webin account authorised to edit those records. Study results are reviewed "
            "in batches of at most 1,000."
        )
    with st.expander("Is there a limit on how many MAGs I can submit at once?"):
        st.write(
            "The guided web workflow accepts up to 1,000 MAGs per job. The CLI is intended for larger "
            "collections: it automatically divides metadata registration into batches of up to 1,000 "
            "and keeps Webin-CLI assembly uploads serial for reliability."
        )
    with st.expander("What is included in the unified provenance ledger?"):
        st.write(
            "The TSV and JSON ledgers cover parent-record autofill, file imports, manual cell changes, "
            "taxonomy assistance, environmental proposals and curator decisions, validation, FASTA "
            "preparation, generated artifacts, and ENA outcomes. Sensitive credentials are excluded."
        )
    with st.expander("Can I password-protect my submission results?"):
        st.write(
            "Yes. Create an optional results password in the final submission form. After processing, "
            "MAG2ENA encrypts the complete job output—including XML, manifests, accession logs, and the "
            "result summary—and removes the plaintext copies. The Jobs page requires the same password "
            "to decrypt the results temporarily for viewing. MAG2ENA does not retain the password, so it "
            "cannot be recovered or reset."
        )
    with st.expander("How is my data handled?"):
        st.write(
            "Uploaded FASTA files are written to a temporary directory and deleted automatically after the "
            "submission job completes. Webin credentials are never written to the MAG2ENA database, job "
            "directory, or downloadable logs. They are held temporarily while the job waits in the background "
            "queue and removed from the persisted job as soon as a worker starts it. If you create a results "
            "password, completed job artifacts are stored only as an authenticated encrypted archive; plaintext "
            "is extracted temporarily after a successful unlock and removed when you leave the Jobs page. Do "
            "not submit data you are not authorised to make public."
        )


def runUI():
    inject_css()

    with st.expander("Using the platform"):
        st.info("""
            This page explains what MAG2ENA automates, where curator decisions remain essential,
            how submission safeguards work, and what artifacts are produced.

            For questions not covered here, contact the corresponding authors:
            brenoslivio@usp.br, anderson.santos@ufz.de, ulisses.rocha@ufz.de.
        """)

    usage_guide()
    faq()
    tutorials()

if __name__ == "__main__":
    runUI()
