import shutil
import urllib.request as request
from contextlib import closing
from urllib.error import URLError
import gzip
import os
import sys

print('Downloading genomes...')
with open("/home/andrey/kraken2/viral/library_report.tsv", "r") as f:
    f.readline() # skip header line
    line = f.readline()
    while(line):
        # viral   >NC_073604.1 Pseudomonas phage SPA05, complete genome   ftp://ftp.ncbi.nlm.nih.gov/...
        seqid = line.split()[1][1:]
        url = line.split()[-1]
        newpath = '/home/andrey/kraken2/viral/genomes/'+seqid
        if not os.path.exists(newpath):
            os.makedirs(newpath)
        try:
            with closing(request.urlopen(url)) as r:
                with open('/home/andrey/kraken2/viral/genomes/'+seqid+'/genome.fna.gz', 'wb') as g:
                    shutil.copyfileobj(r, g)
                g.close()
        except URLError as e:
            if e.reason.find('No such file or directory') >= 0:
                raise Exception('FileNotFound')
            else:
                raise Exception(f'Something else happened. "{e.reason}"')
        with gzip.open('/home/andrey/kraken2/viral/genomes/'+seqid+'/genome.fna.gz', 'rb') as g:
            file_content = g.read()
            with open('/home/andrey/kraken2/viral/genomes/'+seqid+'/genome.fna', 'wb') as h:
                h.write(file_content)
            h.close()
        g.close()
        os.remove('/home/andrey/kraken2/viral/genomes/'+seqid+'/genome.fna.gz')
        sys.stdout.write(f'\033[K{ i } of { n_lines - 1 } genomes\r')
        line = f.readline()
f.close()
sys.stdout.write('\n')
print('Done!')
