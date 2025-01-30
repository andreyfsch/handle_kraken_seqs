import pandas as pd
import csv
from datasets import Dataset, DatasetDict

virus_df = pd.read_csv('/home/andrey/kraken2/viral/virus_seqs_100.csv')
virus_dataset = Dataset.from_pandas(virus_df)

virus_train_testvalid = virus_dataset.train_test_split(test_size=0.1)
virus_test_valid = virus_train_testvalid['test'].train_test_split(test_size=0.5)

virus_train_test_valid_dataset = DatasetDict({
    'train': virus_train_testvalid['train'],
    'test': virus_test_valid['test'],
    'valid': virus_test_valid['train']})

virus_train_df = virus_train_test_valid_dataset['train'].to_pandas()
virus_test_df = virus_train_test_valid_dataset['test'].to_pandas()
virus_valid_df = virus_train_test_valid_dataset['valid'].to_pandas()

virus_train_df.to_csv('/home/andrey/kraken2/viral/virus_train_100.csv', index=False)
virus_test_df.to_csv('/home/andrey/kraken2/viral/virus_test_100.csv', index=False)
virus_valid_df.to_csv('/home/andrey/kraken2/viral/virus_valid_100.csv', index=False)

print('Generating tsv with 6-mer files for DNABERT-6...')
for filename in ['virus_train_100', 'virus_test_100', 'virus_valid_100']:
    with open('/home/andrey/kraken2/viral/'+filename+'.csv', 'r') as csv_file, open('/home/andrey/kraken2/viral/'+filename+'.tsv', 'w', newline='') as tsv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        tsv_writer = csv.writer(tsv_file, delimiter='\t')
        first = True
        for row in csv_reader:
            if not first:
                row[0] = ' '.join([row[0][i:i+6] for i in range(0, len(row[0]), 6)])
            else:
                first = False
            tsv_writer.writerow(row)
print('Done!')
