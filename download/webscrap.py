import requests
from bs4 import BeautifulSoup
import pandas as pd
import logging

logger = logging.getLogger(__name__)

KRAKEN2_INDEXES_URL = "https://benlangmead.github.io/aws-indexes/k2"
BASE_URL = "https://benlangmead.github.io"


def fetch_kraken2_collections() -> pd.DataFrame:
    """
    Scrape the main Kraken2 collections table from the Ben Langmead AWS page.

    Returns
    -------
    pd.DataFrame
        Table with columns: 'Collection', 'Contains', 'Date',
        'Archive size (GB)', 'Index size (GB)', 'tarball_url',
        'md5_url', 'library_report_url'.
    """
    logger.info(
        f"Fetching Kraken2 collections table from: {KRAKEN2_INDEXES_URL}")
    resp = requests.get(KRAKEN2_INDEXES_URL)
    resp.raise_for_status()
    soup = BeautifulSoup(resp.text, "html.parser")

    # Find table with "Collection" header
    tables = soup.find_all("table")
    table = None
    for t in tables:
        headers = [th.get_text(strip=True) for th in t.find_all("th")]
        if "Collection" in headers:
            table = t
            break
    if table is None:
        raise RuntimeError(
            "Could not find Kraken2 Refseq collection table on website!")

    # Get column names
    columns = [th.get_text(strip=True) for th in table.find_all("th")]
    rows = []
    for row in table.find_all("tr")[1:]:
        tds = row.find_all(["td"])
        if not tds or len(tds) < 2:
            continue
        vals = []
        links = {"tarball_url": None, "md5_url": None,
                 "library_report_url": None}
        for i, td in enumerate(tds):
            # Try to find links
            a_tags = td.find_all("a")
            if a_tags:
                for a in a_tags:
                    href = a.get("href", "")
                    txt = a.get_text(strip=True).lower()
                    if href.endswith(".tar.gz"):
                        links["tarball_url"] = href
                    elif href.endswith(".md5"):
                        links["md5_url"] = href
                    elif href.endswith(".tsv") or "library" in txt:
                        links["library_report_url"] = href
            # Always add cell text
            vals.append(td.get_text(strip=True))
        # If missing, fill as empty
        links = {k: v if v else "" for k, v in links.items()}
        # Build row: columns + urls
        # If column count mismatch, pad
        while len(vals) < len(columns):
            vals.append("")
        row_dict = {c: v for c, v in zip(columns, vals)}
        row_dict.update(links)
        rows.append(row_dict)
    if not rows:
        # Instead of raising, return an empty DataFrame with expected columns
        return pd.DataFrame(
            columns=columns + [
                'tarball_url',
                'md5_url',
                'library_report_url'
            ]
        )
    df = pd.DataFrame(rows)
    # Fallback: if no 'library_report_url', but 'library_url' exists, use that
    if "library_report_url" not in df.columns and "library_url" in df.columns:
        df["library_report_url"] = df["library_url"]
    return df
