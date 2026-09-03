from datetime import datetime

import modules
import streamlit as st
import utils
from streamlit_option_menu import option_menu


def clear_job_state():
    """Clear the active Jobs lookup without affecting submission progress."""
    modules.jobs.clear_decrypted_job()
    st.session_state.pop("job_id", None)
    st.session_state.pop("job_input", None)


@st.dialog("ℹ️ Use notice", width="medium")
def cookie_dialog():
    st.markdown(
        """
        This web server uses **session cookies solely** to ensure proper functionality.

        No personal tracking or persistent cookies are employed.

        This platform is released under the [MIT license](https://opensource.org/licenses/MIT). **Free for academic and commercial use**.
        """
    )

    st.session_state["cookie"] = True

    if st.button("Got it", width="stretch", type="primary"):
        st.rerun()


def runUI():
    st.set_page_config(page_title = "MAG2ENA", page_icon = "imgs/mag2ena_logo.png", initial_sidebar_state = "expanded", layout="wide")

    utils.inject_css()

    manual_idx = st.session_state.pop("_nav_index", None)
    page = option_menu(
        None,
        ["Home", "Jobs", "Update", "Help & Tutorials", "About"],
        icons=["house", "gear-wide", "pencil-square", "question-circle", "info-circle"],
        menu_icon="cast",
        default_index=0,
        orientation="horizontal",
        manual_select=manual_idx,
    )

    if page != "Jobs":
        clear_job_state()
    if page != "Update":
        modules.modify.clear_session_state()

    if page == "Home":
        modules.home.runUI()

    elif page == "Update":
        modules.modify.runUI()

    elif page == "Jobs":
        modules.jobs.runUI()

    elif page == "About":
        modules.about.runUI()

    elif page == "Help & Tutorials":
        modules.help.runUI()

    # Show dialog once per session
    if "cookie" not in st.session_state:
        cookie_dialog()

    st.markdown(
        f"""
        <hr>
        <div style="text-align:center; font-size: 0.9em;">
        © {datetime.now().year} MAG2ENA — Released under the 
        <a href="https://opensource.org/licenses/MIT" target="_blank">MIT License</a>. Free for academic and commercial use.
        </div>
        <br>
        """,
        unsafe_allow_html=True
    )

if __name__ == "__main__":
    runUI()
