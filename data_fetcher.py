import gzip
import urllib
from collections import namedtuple

Gene = namedtuple('Gene', ['id', 'name'])


class WBDataFetcher:
    """data fetcher for WormBase raw files"""
    def __init__(self, raw_files_source: str, release_version: str, species: str, project_id: str):
        """create a new data fetcher

        :param raw_files_source: base url where to fetch the raw files
        :type raw_files_source: str
        :param release_version: WormBase release version for the input files
        :type release_version: str
        :param species: WormBase species to fetch
        :type species: str
        :param project_id: project id associated with the species
        :type project_id: str
        """
        self.raw_files_source = raw_files_source
        self.release_version = release_version
        self.species = species
        self.project_id = project_id

    def get_gene_data(self, include_dead_genes: bool = False) -> Gene:
        """fetch gene data
        :param include_dead_genes: whether to include dead genes in the results
        :type include_dead_genes: bool
        :return: data for one gene at each call, including gene_id and gene_name
        :rtype: Gene"""
        address = self.raw_files_source + '/' + self.release_version + '/species/' + self.species + '/' + \
                  self.project_id + '/annotation/' + self.species + '.' + self.project_id + '.' + self.release_version \
                  + '.geneIDs.txt.gz'
        with urllib.request.urlopen(address) as url:
            gzip_file = gzip.GzipFile(fileobj=url)
            for line in gzip_file:
                fields = line.decode("utf-8").strip().split(',')
                if include_dead_genes or fields[4] != "Dead":
                    yield Gene(fields[1], fields[3])

