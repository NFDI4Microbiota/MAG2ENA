# modules/about.py
import streamlit as st
from utils.css_injection import inject_css


def runUI():
    inject_css()

    with st.container(border=True):
        st.markdown("**About MAG2ENA**")
        st.markdown(
            """
            MAG2ENA supports submission and post-submission metadata maintenance for
            metagenome-assembled genomes (MAGs) in the European Nucleotide Archive (ENA).
            It combines parent-record metadata autofill, CheckM/CheckM2 and GTDB-Tk imports,
            curator-reviewed environmental ontology resolution, ERC000047 validation, FASTA
            preparation, provenance, and Webin submission without requiring users to write XML
            or run Webin-CLI manually.

            The web application supports guided submissions of up to 1,000 MAGs per job and
            metadata-only updates to existing MAG samples without replacing their assemblies. The
            MAG2ENA command-line interface provides the same metadata assistance and ENA
            submission pathway for larger, scriptable workflows, without web provenance or Jobs tracking.
            """
        )

        col_github, col_ena = st.columns(2)
        with col_github:
            st.link_button(
                "GitHub",
                "https://github.com/NFDI4Microbiota/MAG2ENA",
                width="stretch",
            )
        with col_ena:
            st.link_button(
                "ENA Webin Portal",
                "https://www.ebi.ac.uk/ena/browser/submit",
                width="stretch",
            )

    with st.container(border=True):
        st.markdown("**Authors**")
        st.markdown(
            """
            This platform was developed by:

            - **Breno L. S. de Almeida** (corresponding author; [brenoslivio@usp.br](mailto:brenoslivio@usp.br))
            - **Anderson P. Avila-Santos** (corresponding author; [anderson.santos@ufz.de](mailto:anderson.santos@ufz.de))
            - Martin Bole
            - Sanchita Kamath
            - Robson P. Bonidia
            - Peter F. Stadler
            - André C. P. L. F. de Carvalho
            - **Ulisses N. da Rocha** (corresponding author; [ulisses.rocha@ufz.de](mailto:ulisses.rocha@ufz.de))
            """
        )

    with st.container(border=True):
        st.markdown("**Acknowledgements**")
        st.markdown(
            """
            This work was supported by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation)
            through NFDI4Microbiota – National Research Data Infrastructure for Microbiota Research, project
            number 460129525. This work has been funded by the Canadian International Development Research
            Centre (IDRC) under the Grant Agreement 109981, and the UK government's Foreign, Commonwealth and
            Development Office. The views expressed here do not necessarily reflect those of the UK government's
            Foreign, Commonwealth and Development Office, IDRC, or IDRC's Board of Governors. Breno L. S. de
            Almeida has been funded by the São Paulo Research Foundation (FAPESP), grant #2024/10958-1, and the
            Google PhD Fellowship. This project (ZT-I-PF-3-108) was funded by the Initiative and Networking Fund
            of the Helmholtz Association in the framework of the Helmholtz Metadata Collaboration project call.

            We also acknowledge the open-source libraries and tools that made this work possible.
            """
        )

    with st.container(border=True):
        st.markdown("**Data Availability**")
        st.markdown(
            """
            The web application is available at
            [web.app.ufz.de/mag2ena](https://web.app.ufz.de/mag2ena). Source code and the CLI are
            available under the MIT License at
            [github.com/NFDI4Microbiota/MAG2ENA](https://github.com/NFDI4Microbiota/MAG2ENA).
            """
        )

    with st.container(border=True):
        st.markdown("**Contact**")
        st.markdown(
            """
            Questions or issues? Reach out to the corresponding authors:

            - [brenoslivio@usp.br](mailto:brenoslivio@usp.br)
            - [anderson.santos@ufz.de](mailto:anderson.santos@ufz.de)
            - [ulisses.rocha@ufz.de](mailto:ulisses.rocha@ufz.de)
            """
        )


if __name__ == "__main__":
    runUI()
