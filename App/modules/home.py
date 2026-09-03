# modules/home.py
import base64
import streamlit as st
from utils.css_injection import inject_css


def runUI():
    """Render the value proposition and the MAG submission workflow."""
    inject_css()

    with open("imgs/mag2ena_logo.png", "rb") as f:
        logo_b64 = base64.b64encode(f.read()).decode("utf-8")

    st.markdown(
        f"""
        <div class="container-max">
          <section class="hero">
            <img src="data:image/png;base64,{logo_b64}" alt="MAG2ENA logo" width="420">
            <p class="hero-sub">From fragmented bioinformatics outputs to reviewed, traceable ENA records</p>
            <span class="badge">Genomic Standards Consortium package for MAGs (MIMAGS) supported</span>
            <span class="badge">ENA Webin-CLI 9.0.3</span>
          </section>
        </div>
        """,
        unsafe_allow_html=True,
    )

    _cards = [
        (
            "🧬",
            "Reuse, do not retype",
            "Recover parent-sample metadata from ENA and import CheckM/CheckM2 and GTDB-Tk outputs.",
        ),
        (
            "🌍",
            "Resolve with evidence",
            "Generate field-aware ENVO proposals with confidence, alternatives, provenance, and curator approval.",
        ),
        (
            "🛡️",
            "Catch failures early",
            "Preflight ERC000047 fields, species-rank taxonomy, controlled vocabularies, and FASTA structure.",
        ),
        (
            "🚀",
            "Submit and maintain",
            "Submit new MAGs or update existing ENA metadata with validation, provenance, and rollback support.",
        ),
    ]
    _columns = st.columns(4)
    for _col, (_icon, _title, _text) in zip(_columns, _cards):
        with _col:
            st.markdown(
                f'<div class="feature-card">'
                f'<span class="feature-icon">{_icon}</span>'
                f'<div class="feature-title">{_title}</div>'
                f'<div class="feature-text">{_text}</div>'
                f'</div>',
                unsafe_allow_html=True,
            )

    st.divider()

    st.info(
        "**Already submitted your MAGs?** Open **Update** to retrieve them by study or "
        "sample accession, improve their metadata, and submit metadata-only changes without "
        "uploading the assemblies again."
    )

    st.markdown("### Demonstrated on real ENA submissions")
    metric_cols = st.columns(4)
    metric_cols[0].metric("MAGs accepted", "100 / 100")
    metric_cols[0].caption("Marine and human-gut evaluation set")
    metric_cols[1].metric("Environmental values recovered", "204 / 300")
    metric_cols[1].caption("68.0% filled from parent records")
    metric_cols[2].metric("Review units", "126 vs 300")
    metric_cols[2].caption("58.0% fewer through parent-aware grouping")
    metric_cols[3].metric("Resolver proposals", "60 / 126")
    metric_cols[3].caption("34 in the high-confidence tier")
    st.caption(
        "Figures summarize the MAG2ENA evaluation described in the accompanying manuscript; "
        "proposal confidence supports review and is not a substitute for curator judgement."
    )

    st.divider()

    with st.expander("Why MAG2ENA?", expanded=False):
        st.markdown(
            "MAG submission usually means reconciling a checklist, parent-sample records, "
            "quality reports, taxonomy output, ontology terms, XML, manifests, and Webin-CLI. "
            "MAG2ENA connects those steps while keeping the researcher in control of semantic decisions."
        )

        st.markdown("#### What makes the workflow different")
        st.markdown(
            """
            - **Parent-aware curation:** ENA sample XML is queried first. If it is unavailable, MAG2ENA can conservatively recover values that agree across assemblies exactly linked to that sample; conflicting values are omitted.
            - **Field-aware environmental resolution:** broad-scale context, local context, and environmental medium are resolved as distinct MIxS roles rather than receiving the same free-text term.
            - **Curator control:** ontology candidates are proposals, not silent replacements. Users see the evidence and confidence tier before applying a term.
            - **Submission-grade preflight:** checklist constraints and taxonomy policy are checked before credentials or files are sent to ENA.
            - **Metadata maintenance:** existing MAG samples can be retrieved, edited, validated, and updated without replacing their assemblies.
            - **End-to-end provenance:** autofill, imports, manual edits, resolver decisions, validation, file preparation, and ENA outcomes are captured in downloadable TSV and JSON ledgers.
            """
        )

        st.markdown("#### How it works")
        st.markdown(
            """
            1. Load the two-MAG example, download the template, or upload an existing metadata TSV.
            2. Retrieve reusable collection, location, source, and environmental metadata from each parent ENA sample.
            3. Import CheckM/CheckM2 quality estimates and GTDB-Tk classifications, then resolve species-level NCBI taxonomy IDs.
            4. Review grouped, field-specific environmental ontology proposals and apply only the terms you approve.
            5. Validate every row against ERC000047, ENA taxonomy, and MAG2ENA's pre-submission safeguards.
            6. Upload matching `.fasta.gz` assemblies. Complete one-sequence chromosomes are handled explicitly; incomplete single-contig MAGs are blocked pending an ENA exception.
            7. Submit to Testing or Production and use the Job ID to retrieve accessions, logs, manifests, XML, and provenance.
            """
        )
        st.info(
            "Start with the **Testing** portal. Test records are temporary, allowing the complete "
            "workflow to be checked before a production deposition. Newly registered Webin accounts "
            "may also require time to activate."
        )

        st.markdown("#### Web for guided curation, CLI for pipelines")
        st.info(
            "The web app supports guided jobs of up to **1,000 MAGs** with provenance and "
            "asynchronous Jobs tracking. The MAG2ENA CLI provides the same metadata "
            "assistance and ENA pathway without provenance, automatically batches metadata "
            "registration, and submits assemblies serially for reliability."
        )
        st.code(
            """python mag2ena.py \\
  --metadata metadata.tsv \\
  --fasta-dir fasta/ \\
  --study-accession PRJEBXXXX""",
            language="bash",
        )
        st.caption(
            "After local preflight, the CLI securely prompts for the Webin username and a "
            "non-echoed password; credentials are not accepted as command-line arguments."
        )
        st.markdown(
            "See the [README on GitHub](https://github.com/NFDI4Microbiota/MAG2ENA) "
            "for metadata-preparation commands, non-interactive secret handling, all CLI "
            "options, and deployment guidance."
        )

        st.markdown(
            """
            **References**
            - [GSC MIMAG checklist: ERC000047](https://www.ebi.ac.uk/ena/browser/view/ERC000047)
            - [MIxS Term Browser](https://w3id.org/mixs/)
            - [ENA Webin Submission Portal](https://www.ebi.ac.uk/ena/browser/submit)
            - [ENVO Ontology](https://www.ebi.ac.uk/ols4/ontologies/envo)
            """
        )

    st.markdown("## Start a submission")
    st.caption(
        "Try the bundled two-MAG example for a complete, previously test-validated walkthrough, "
        "or begin with your own metadata and assemblies."
    )

    from modules import submit

    submit.runUI()

if __name__ == "__main__":
    runUI()
