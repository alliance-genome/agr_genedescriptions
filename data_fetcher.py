import gzip
import urllib
from collections import namedtuple, defaultdict
from namedlist import namedlist
from typing import List
from descriptions_rules import GO_ASPECT

Gene = namedtuple('Gene', ['id', 'name'])
GOAnnotation = namedtuple('GOAnnotation', ['go_id', 'qualifier', 'paper_reference', 'evidence_code', 'aspect',
                                           'annotation_ext', 'go_name', 'is_obsolete'])
GOOntologyEntry = namedlist('GOOntologyEntry', 'name is_obsolete')


class WBRawDataFetcher:
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
        self.go_data = defaultdict(set)
        self.go_ontology = {}

    def get_gene_data(self, include_dead_genes: bool = False, include_pseudo_genes: bool = False) -> Gene:
        """get all gene data from the fetcher, returning one gene per call

        :param include_dead_genes: whether to include dead genes in the results
        :type include_dead_genes: bool
        :param include_pseudo_genes: whether to include pseudo genes in the results
        :type include_dead_genes: bool
        :return: data for one gene per each call, including gene_id and gene_name
        :rtype: Gene
        """
        address = self.raw_files_source + '/' + self.release_version + '/species/' + self.species + '/' + \
                  self.project_id + '/annotation/' + self.species + '.' + self.project_id + '.' + self.release_version \
                  + '.geneIDs.txt.gz'
        # TODO if include_pseudo_genes == False exclude it - get pseudo gene data from another file
        with urllib.request.urlopen(address) as url:
            gzip_file = gzip.GzipFile(fileobj=url)
            for line in gzip_file:
                fields = line.decode("utf-8").strip().split(',')
                if include_dead_genes or fields[4] != "Dead":
                    name = fields[2] if fields[2] != '' else fields[3]
                    yield Gene(fields[1], name)

    def read_go_data(self) -> None:
        """read go data and gene ontology. After calling this function, go annotations containing mapped go names can
        be retrieved by using the :meth:`data_fetcher.WBRawDataFetcher.get_go_annotations` function
        """
        self._read_go_ontology()
        address = self.raw_files_source + '/' + self.release_version + '/species/' + self.species + '/' + \
                  self.project_id + '/annotation/' + self.species + '.' + self.project_id + '.' + self.release_version \
                  + '.go_annotations.gaf.gz'
        with urllib.request.urlopen(address) as url:
            gzip_file = gzip.GzipFile(fileobj=url)
            for line in gzip_file:
                line = line.decode("utf-8")
                if not line.startswith("!"):
                    fields = line.strip("\n").split('\t')
                    go_aspect = None
                    if fields[8] == 'C':
                        go_aspect = GO_ASPECT.CELLULAR_COMPONENT
                    elif fields[8] == 'F':
                        go_aspect = GO_ASPECT.MOLECULAR_FUNCTION
                    elif fields[8] == 'P':
                        go_aspect = GO_ASPECT.BIOLOGICAL_PROCESS
                    is_obsolete = False
                    if self.go_ontology[fields[4]].is_obsolete == 'true':
                        is_obsolete = True
                    self.go_data[fields[1]].add(GOAnnotation(fields[4], fields[3], fields[5], fields[6], go_aspect,
                                                                fields[15], self.go_ontology[fields[4]].name,
                                                                is_obsolete))

    def _read_go_ontology(self):
        """read go ontology data"""
        address = self.raw_files_source + '/' + self.release_version + '/ONTOLOGY/gene_ontology.' + \
                  self.release_version + '.obo'
        with urllib.request.urlopen(address) as url:
            go_id = None
            go_entry = GOOntologyEntry(None, None)
            for line in url:
                line = line.decode("utf-8")
                if line.strip() == '[Term]':
                    if go_entry.name is not None and go_id is not None:
                        self.go_ontology[go_id] = go_entry
                        go_id = None
                        go_entry = GOOntologyEntry(None, None)
                else:
                    fields = line.strip().split(": ")
                    if fields[0] == "id":
                        go_id = fields[1]
                    elif fields[0] == "name":
                        go_entry.name = fields[1]
                    elif fields[0] == "is_obsolete":
                        go_entry.is_obsolete = fields[1]

    def get_go_annotations(self, geneid: str, go_aspect: GO_ASPECT, include_obsolete: bool = False,
                           priority_list: tuple = ("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "ISS", "ISO", "ISA", "ISM",
                                                   "IGC", "IBA", "IBD", "IKR", "IRD", "RCA",
                                                   "IC")) -> List[GOAnnotation]:
        """retrieve go annotations for a given gene id and for a given aspect. The annotations are unique for each pair
        <gene_id, go_term_id>. This means that when multiple annotations for the same pair are found in the go data, the
        one with the evidence code with highest priority is returned (see the *priority_list* parameter to set the
        priority according to evidence codes)

        :param geneid: the id of the gene related to the annotations to retrieve, in standard format
        :type geneid: str
        :param go_aspect: the aspect of the annotations to retrieve
        :type go_aspect: GO_ASPECT
        :param include_obsolete: whether to include obsolete annotations
        :type include_obsolete: bool
        :param priority_list: the priority list for the evidence codes. If multiple annotations with the same go_term
            are found, only the one with highest priority is returned. The first element in the list has the highest
            priority, whereas the last has the lowest. Only annotations with evidence codes in the priority list are
            returned. All other annotations are ignored
        :type priority_list: List[str]
        :return: the list of go annotations for the given gene
        :rtype: List[GOAnnotation]
        """
        priority_map = dict(zip(priority_list, reversed(range(len(priority_list)))))
        annotations = [annotation for annotation in self.go_data[geneid] if annotation.aspect == go_aspect and
                       (include_obsolete or not annotation.is_obsolete)]
        go_id_selected_annotation = {}
        for annotation in annotations:
            if annotation.evidence_code in priority_map.keys():
                if annotation.go_id in go_id_selected_annotation:
                    if priority_map[annotation.evidence_code] > \
                            priority_map[go_id_selected_annotation[annotation.go_id].evidence_code]:
                        go_id_selected_annotation[annotation.go_id] = annotation
                else:
                    go_id_selected_annotation[annotation.go_id] = annotation

        return [annotation for annotation in go_id_selected_annotation.values()]

