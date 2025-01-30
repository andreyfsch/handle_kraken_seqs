import argparse
import os

import webscrap

home = os.path.expanduser('~')

parser = argparse.ArgumentParser(prog="handle_kraken_seqs", usage="python -m %(prog)s [options]")
parser.add_argument("-p", "--path", type=str, default=home+"/kraken2", help="Path of the kraken2 folder")
parser.add_argument("-c", "--collection", type=str, default="viral", help="Name of the kraken2 collection", choices=webscrap.collections.keys())
parser.add_argument("-d", "--download-genomes", action="store_true", help="Wether to download genome files")

parser.parse_args()
