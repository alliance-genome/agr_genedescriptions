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
