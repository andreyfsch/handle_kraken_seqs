import pytest
import pandas as pd
from download import webscrap


@pytest.mark.slow
def test_fetch_kraken2_collections_live():
    """
    Integration test: fetch the live Kraken2 collections page,
    parse the table, and validate core structure.
    """
    # Skip test if no internet
    import socket
    try:
        socket.create_connection(("benlangmead.github.io", 443), timeout=5)
    except OSError:
        pytest.skip("No internet connection")

    df = webscrap.fetch_kraken2_collections()
    # Check it is a DataFrame and has expected columns
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert "Collection" in df.columns
    # Check at least one row has the expected Kraken2 archive URLs
    url_cols = ["tarball_url", "md5_url", "library_report_url"]
    for col in url_cols:
        assert col in df.columns
        # Should be non-null for at least one collection
        assert df[col].notnull().any()
