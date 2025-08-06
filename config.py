import logging
import random

KRAKEN_PATH = "/home/andrey/generate_kraken_dataset/kraken2"
KRAKEN_DATABASE = "viral"
CSV_OUTPUT_PATH = "."
MIN_SUBSEQ_LEN = 100
MAX_SUBSEQ_LEN = 512
COMPRESSED_CSVS = True
MIN_NUM_SEQS = 5000
PERCENTAGE = 98.0
RNG_SEED = random.Random(42)
LOGGING_LEVEL = logging.INFO
