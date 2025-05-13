import json
import logging
import os
import urllib.request

import requests
from fastapi_okta.okta_utils import get_authentication_token, generate_headers

ATEAM_API = os.environ.get("ATEAM_API", "https://curation.alliancegenome.org/api")


logger = logging.getLogger(__name__)


def get_anatomy_ontologies_roots():
    url = f'{ATEAM_API}/anatomicalterm/rootNodes'
    token = get_authentication_token()
    headers = generate_headers(token)
    try:
        get_request = urllib.request.Request(url=url, method='GET', headers=headers)
        with urllib.request.urlopen(get_request) as get_response:
            if get_response.getcode() == 200:
                logger.debug("Request successful")
                res = get_response.read().decode('utf-8')
                json_res = json.loads(res)
                return [entity for entity in json_res["entities"] if not entity["obsolete"]]
            else:
                logger.error("Request error")
                return False
    except requests.exceptions.RequestException as e:
        logger.error(f"Error occurred: {e}")
        return False


def get_ontology_node_children(node_curie: str):
    url = f'{ATEAM_API}/ontologyterm/{node_curie}/children'
    token = get_authentication_token()
    headers = generate_headers(token)
    try:
        get_request = urllib.request.Request(url=url, method='GET', headers=headers)
        with urllib.request.urlopen(get_request) as get_response:
            if get_response.getcode() == 200:
                logger.debug("Request successful")
                res = get_response.read().decode('utf-8')
                json_res = json.loads(res)
                if "entities" not in json_res:
                    logger.warning("No entities found in response for ontology get children api")
                    return []
                return [entity for entity in json_res["entities"] if not entity["obsolete"]]
            else:
                logger.error("Request error")
                return False
    except requests.exceptions.RequestException as e:
        logger.error(f"Error occurred: {e}")
        return False


def get_expression_annotations_from_api(data_provider: str):
    """Get expression annotations from the A-team API."""
    token = get_authentication_token()
    headers = generate_headers(token)
    data = {"expressionAnnotationSubject.dataProvider.abbreviation": data_provider}
    page = 0
    page_size = 5000
    annotations = []
    try:
        while True:
            url = f'{ATEAM_API}/gene-expression-annotation/findForPublic?limit={page_size}&page={page}&view=ForPublic'
            get_request = urllib.request.Request(url=url, method='POST', headers=headers,
                                                 data=json.dumps(data).encode('utf-8'))
            with urllib.request.urlopen(get_request) as get_response:
                if get_response.getcode() == 200:
                    logger.debug("Request successful")
                    res = get_response.read().decode('utf-8')
                    json_res = json.loads(res)
                    if json_res["returnedRecords"] == 0:
                        break
                    annotations.extend([{"gene_id": row["expressionAnnotationSubject"]["primaryExternalId"],
                                         "gene_symbol": row["expressionAnnotationSubject"]["geneSymbol"]["displayText"],
                                         "anatomy_id": row["expressionPattern"]["whereExpressed"][
                                             "anatomicalStructure"]["curie"]} for row in json_res["results"]
                                        if "expressionPattern" in row
                                        ])
                    page += 1
                else:
                    logger.error("Request error")
                    return False
        return annotations
    except requests.exceptions.RequestException as e:
        logger.error(f"Error occurred: {e}")
        return False


def get_data_providers_from_api():
    """Get data providers from the A-team API."""
    url = f'{ATEAM_API}/species/findForPublic?limit=100&page=0&view=ForPublic'
    token = get_authentication_token()
    headers = generate_headers(token)
    try:
        get_request = urllib.request.Request(url=url, method='GET', headers=headers)
        with urllib.request.urlopen(get_request) as get_response:
            if get_response.getcode() == 200:
                logger.debug("Request successful")
                res = get_response.read().decode('utf-8')
                json_res = json.loads(res)
                return [species["abbreviation"] for species in json_res["results"]]
            else:
                logger.error("Request error")
                return False
    except requests.exceptions.RequestException as e:
        logger.error(f"Error occurred: {e}")
        return False


def get_gene_data_from_api(data_provider: str):
    """Get gene data from the A-team API."""
    page = 0
    page_size = 5000
    token = get_authentication_token()
    headers = generate_headers(token)
    req_data = {"dataProvider.abbreviation": data_provider}
    genes = []
    try:
        while True:
            url = f'{ATEAM_API}/gene/find?limit={page_size}&page={page}'
            request = urllib.request.Request(url=url, method='POST', headers=headers,
                                             data=json.dumps(req_data).encode('utf-8'))
            with urllib.request.urlopen(request) as get_response:
                if get_response.getcode() == 200:
                    logger.debug("Request successful")
                    res = get_response.read().decode('utf-8')
                    json_res = json.loads(res)
                    if json_res["returnedRecords"] == 0:
                        break
                    genes.extend([{"gene_id": row["primaryExternalId"], "gene_symbol": row["geneSymbol"]["displayText"]}
                                  for row in json_res["results"]])
                    page += 1
                else:
                    logger.error("Request error")
                    return False
        return genes
    except requests.exceptions.RequestException as e:
        logger.error(f"Error occurred: {e}")
        return False
