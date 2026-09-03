def inject_css(*args, **kwargs):
    """Load Streamlit-only styling lazily so core helpers remain CLI-safe."""
    from .css_injection import inject_css as _inject_css

    return _inject_css(*args, **kwargs)
