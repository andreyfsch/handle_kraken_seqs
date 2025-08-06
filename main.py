"""
main.py

Demonstrates programmatic invocation of the Kraken2 CLI pipeline using
the run_cli function from cli.py. This allows scripting of all dataset
management and generation steps from Python, rather than the shell.

Each CLI command can be called by passing a list of arguments just as
you would on the command line, but as Python lists.

Examples
--------
# List available collections
>>> from cli import run_cli
>>> run_cli(['list-collections'])

# Download a Kraken2 database
>>> run_cli([
...     'download-db',
...     '--collection', 'PlusPF',
...     '--outdir', './kraken_downloads'
... ])

# Download all genomes listed in the library_report.tsv
>>> run_cli([
...     'download-genomes',
...     '--library-report', './kraken_downloads/library_report.tsv',
...     '--output-dir', './kraken_downloads/genomes',
...     '--threads', '8'
... ])

# Generate CSV datasets from genomes
>>> run_cli([
...     'generate-csv',
...     '--min-subseq-len', '100',
...     '--max-subseq-len', '400',
...     '--compression',
...     '--parallel'
... ])

Notes
-----
All CLI commands and options are available programmatically.
"""

from cli import run_cli


def main():
    # Example: List available collections
    run_cli(['list-collections'])

    # Example: Download a Kraken2 database (uncomment to use)
    # run_cli([
    #     'download-db',
    #     '--collection', 'PlusPF',
    #     '--outdir', './kraken_downloads'
    # ])

    # Example: Download genomes (uncomment to use)
    # run_cli([
    #     'download-genomes',
    #     '--library-report',
    #     './kraken_downloads/library_report.tsv',
    #     '--output-dir', './kraken_downloads/genomes',
    #     '--threads', '8'
    # ])

    # Example: Generate CSV datasets (uncomment to use)
    # run_cli([
    #     'generate-csv',
    #     '--min-subseq-len', '100',
    #     '--max-subseq-len', '400',
    #     '--compression',
    #     '--parallel'
    # ])


if __name__ == '__main__':
    main()
