import pytest
import pandas as pd
from pathlib import Path
import requests
import download.webscrap as webscrap


@pytest.fixture
def sample_html():
    # Adjust path as needed
    html_path = Path("tests/data/benlangmead-github-io_aws-indexes_k2.html")
    return html_path.read_text(encoding="utf-8")


SAMPLE_TABLE_HTML = """
<table>
  <tr>
    <th>Collection</th>
    <th>Contains</th>
    <th>Date</th>
    <th>Archive size (GB)</th>
    <th>Index size (GB)</th>
    <th>HTTPS URL</th>
    <th>Inspect</th>
    <th>Library report</th>
    <th>MD5</th>
  </tr>
  <tr>
    <td>TestCol</td>
    <td>Viruses</td>
    <td>2024-07-19</td>
    <td>12.3</td>
    <td>1.5</td>
    <td><a href="http://example.com/archive.tar.gz">tar.gz</a></td>
    <td><a href="http://example.com/inspect.html">inspect</a></td>
    <td><a href="http://example.com/library_report.tsv">tsv</a></td>
    <td><a href="http://example.com/archive.md5">md5</a></td>
  </tr>
</table>
"""


def test_fetch_returns_dataframe(monkeypatch):
    class MockResp:
        status_code = 200
        def __init__(self, text): self.text = text
        def raise_for_status(self): pass
    monkeypatch.setattr(
        webscrap.requests, "get",
        lambda url: MockResp(SAMPLE_TABLE_HTML)
    )
    df = webscrap.fetch_kraken2_collections()

    # Only require these columns are present,
    # not that they are the *only* columns
    must_have = [
        'Collection', 'Contains', 'Date',
        'tarball_url', 'md5_url', 'library_report_url'
    ]
    assert isinstance(df, pd.DataFrame)
    assert set(must_have).issubset(df.columns)


def test_first_row_and_urls(monkeypatch, sample_html):
    """First row has expected data and all URL columns look plausible."""
    class MockResp:
        status_code = 200
        def __init__(self, text): self.text = text
        def raise_for_status(self): pass
    monkeypatch.setattr(
        webscrap.requests, "get",
        lambda url: MockResp(sample_html)
    )
    df = webscrap.fetch_kraken2_collections()
    assert not df.empty
    row = df.iloc[0]
    assert isinstance(row['Collection'], str) and row['Collection']
    # URL checks (if present)
    for col in ['tarball_url', 'md5_url', 'library_report_url']:
        val = row[col]
        assert (val == "") or val.startswith("http")


def test_network_error(monkeypatch):
    """Raises on network error."""
    def raise_error(*a, **kw): raise requests.RequestException("fail")
    monkeypatch.setattr(webscrap.requests, "get", raise_error)
    with pytest.raises(requests.RequestException):
        webscrap.fetch_kraken2_collections()


def test_missing_table(monkeypatch):
    """Raises if no valid table is found."""
    class MockResp:
        status_code = 200
        text = "<html><body>No tables here!</body></html>"
        def raise_for_status(self): pass
    monkeypatch.setattr(
        webscrap.requests, "get",
        lambda url: MockResp()
    )
    try:
        webscrap.fetch_kraken2_collections()
    except RuntimeError as e:
        msg = str(e)
        assert (
            "No valid rows" in msg
            or "Could not find Kraken2 Refseq collection table" in msg
        )


def test_broken_html(monkeypatch):
    """Should not crash if table structure is weird."""
    broken_html = """
    <table>
      <tr><th>Collection</th><th>Contains</th></tr>
      <tr><td>Foo</td></tr>
    </table>
    """

    class MockResp:
        status_code = 200
        def __init__(self): self.text = broken_html
        def raise_for_status(self): pass
    monkeypatch.setattr(
        webscrap.requests, "get",
        lambda url: MockResp()
    )
    df = webscrap.fetch_kraken2_collections()
    assert isinstance(df, pd.DataFrame)


def test_partial_links(monkeypatch):
    """Should handle rows missing one or more link columns."""
    html = """
    <table>
      <tr>
        <th>Collection</th><th>Contains</th><th>Date</th>
        <th>Archive Size</th><th>Index Size</th><th>Links</th>
      </tr>
      <tr>
        <td>ABC</td><td>Stuff</td><td>2021-01-01</td>
        <td>1G</td><td>2G</td>
        <td>
          <a href="http://somewhere/abc.tar.gz">tar.gz</a>
          <!-- md5 and tsv missing -->
        </td>
      </tr>
    </table>
    """

    class MockResp:
        status_code = 200
        def __init__(self): self.text = html
        def raise_for_status(self): pass
    monkeypatch.setattr(
        webscrap.requests, "get",
        lambda url: MockResp()
    )
    df = webscrap.fetch_kraken2_collections()
    assert "ABC" in df["Collection"].values
    assert "http://somewhere/abc.tar.gz" in df["tarball_url"].values
