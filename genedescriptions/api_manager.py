import json
import logging
import os
import ssl
import urllib

from urllib import request

logger = logging.getLogger("API manager")


class APIManager(object):
    def __init__(self, textpresso_api_token):
        self.textpresso_api_token = textpresso_api_token
        self.tpc_cache = {}
        self.class_cache = {}
        self.tpc_api_endpoint = "https://textpressocentral.org:18080/v1/textpresso/api/get_documents_count"
        if not os.environ.get('PYTHONHTTPSVERIFY', '') and getattr(ssl, '_create_unverified_context', None):
            ssl._create_default_https_context = ssl._create_unverified_context

    def get_textpresso_popularity(self, keyword: str):
        """get the number of papers in the C. elegans literature that mention a certain keyword from Textpresso Central API

        Args:
            keyword (str): the keyword to search, or any combination of keywords containing AND and OR operators
        Returns:
            int: the popularity of the specified keyword
        """
        if keyword in self.tpc_cache:
            logger.debug("Popularity for keyword found in cache")
            return self.tpc_cache[keyword]
        else:
            data = json.dumps({"token": self.textpresso_api_token, "query": {
                "keywords": keyword, "type": "document", "corpora": ["C. elegans"]}})
            data = data.encode('utf-8')
            req = urllib.request.Request(self.tpc_api_endpoint, data, headers={'Content-type': 'application/json',
                                                                               'Accept': 'application/json'})
            res = urllib.request.urlopen(req)
            logger.debug("Sending request to Textpresso Central API")
            popularity = int(json.loads(res.read().decode('utf-8')))
            self.tpc_cache[keyword] = popularity
            return popularity

    def get_gene_class(self, gene_id: str):
        """get the gene class of a gene from WormBase API

        Args:
            gene_id (str): the Wormbase WBGene ID of the gene
        Returns:
            str: the class of the gene
        """
        if gene_id in self.class_cache:
            logger.debug("Gene class for gene " + gene_id + " found in cache")
            return self.class_cache[gene_id]
        else:
            try:
                logger.debug("Getting gene class for gene " + gene_id)
                gene_class_data = json.loads(urllib.request.urlopen("http://rest.wormbase.org/rest/field/gene/" +
                                                                    gene_id + "/gene_class").read())
                if "gene_class" in gene_class_data and gene_class_data["gene_class"]["data"] and "tag" in \
                        gene_class_data["gene_class"]["data"] and "label" in \
                        gene_class_data["gene_class"]["data"]["tag"]:
                    result = gene_class_data["gene_class"]["data"]["tag"]["label"]
                    self.class_cache[gene_id] = result
                    return result
            except:
                return None
            return None