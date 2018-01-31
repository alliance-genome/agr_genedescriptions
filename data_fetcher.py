import gzip
import urllib


class DataFetcher:
    def __init__(self, raw_files_source: str, release_version: str, species: str, project_id: str):
        self.raw_files_source = raw_files_source
        self.release_version = release_version
        self.species = species
        self.project_id = project_id

    def get_gene_data(self):
        address = self.raw_files_source + '/' + self.release_version + '/species/' + self.species + '/' + \
                  self.project_id + '/annotation/' + self.species + '.' + self.project_id + '.' + self.release_version \
                  + '.geneIDs.txt.gz'
        with urllib.request.urlopen(address) as url:
            gzipFile = gzip.GzipFile(fileobj=url)
            for line in gzipFile:
                fields = line.decode("utf-8").split(',')
                if fields[4] != "Dead":
                    yield fields[1]

        pass
