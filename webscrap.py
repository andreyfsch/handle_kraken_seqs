try:
    from bs4 import BeautifulSoup
except ImportError:
    import pip
    pip.main(['install', '--user', 'beautifulsoup4'])
    from bs4 import BeautifulSoup
import requests

collections = {}

url = "https://benlangmead.github.io/aws-indexes/k2"

response = requests.get(url)

if response.status_code != 200:
    raise ValueError('Kraken2 website unavailable')

soup = BeautifulSoup(response.text, 'html.parser')

section = soup.find('section')
table = section.find('tbody')

trs = table.find_all('tr')
for tr in trs:
    tds = tr.find_all('td')
    collection_name = tds[0].contents[0].lower().split()[0]
    collections[collection_name] = {}
    
    library_a = tds[7].find('a')
    if library_a:
        library_url = library_a.attrs['href']
        collections[collection_name]['library'] = library_url
    else:
        del collections[collection_name]

    collections[collection_name]['date'] = tds[2].contents[0]

