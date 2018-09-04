import json
import os
import ssl
import inflect
import re
import urllib.request
import urllib.parse

from namedlist import namedlist
from genedescriptions.config_parser import GenedescConfigParser
from genedescriptions.ontology_tools import *
from ontobio.ontol import Ontology
from collections import namedtuple
from enum import Enum


class DataType(Enum):
    GO = 1
    DO = 2
    EXPR = 3


Gene = namedtuple('Gene', ['id', 'name', 'dead', 'pseudo'])


Sentence = namedlist('Sentence', ['prefix', 'terms_ids', 'postfix', 'text', 'aspect', 'evidence_group', 'terms_merged',
                                  'additional_prefix', 'qualifier', 'ancestors_covering_multiple_terms'])


class SingleDescStats(object):
    """statistics for a single gene description"""
    def __init__(self):
        self.total_number_go_annotations = 0
        self.total_number_do_annotations = 0
        self.number_final_do_term_covering_multiple_initial_do_terms = 0
        self.set_initial_go_ids_f = []
        self.set_initial_go_ids_p = []
        self.set_initial_go_ids_c = []
        self.set_final_go_ids_f = []
        self.set_final_go_ids_p = []
        self.set_final_go_ids_c = []
        self.set_initial_do_ids = []
        self.set_final_do_ids = []
        self.set_best_orthologs = []


class GeneDesc(object):
    """gene description"""
    def __init__(self, gene_id: str, gene_name: str = "", description: str = None, go_description: str = None,
                 go_function_description: str = None, go_process_description: str = None,
                 go_component_description: str = None, do_description: str = None, orthology_description: str = None,
                 stats: SingleDescStats = None, publications: str = "", refs: str = ""):
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.description = description
        self.go_description = go_description
        self.go_function_description = go_function_description
        self.go_process_description = go_process_description
        self.go_component_description = go_component_description
        self.do_description = do_description
        self.orthology_description = orthology_description
        self.publications = publications
        self.refs = refs
        if stats:
            self.stats = stats
        else:
            self.stats = SingleDescStats()


class DescriptionsStats(object):
    """overall statistics for a set of gene descriptions"""
    def __init__(self):
        self.total_number_of_genes = 0
        self.number_genes_with_non_null_description = 0
        self.number_genes_with_non_null_go_description = 0
        self.number_genes_with_non_null_go_function_description = 0
        self.number_genes_with_non_null_go_process_description = 0
        self.number_genes_with_non_null_go_component_description = 0
        self.number_genes_with_null_go_description = 0
        self.number_genes_with_more_than_3_initial_go_terms = 0
        self.number_genes_with_non_null_do_description = 0
        self.number_genes_with_null_do_description = 0
        self.number_genes_with_more_than_3_initial_do_terms = 0
        self.number_genes_with_final_do_terms_covering_multiple_initial_terms = 0
        self.average_number_initial_go_terms_f = 0
        self.average_number_initial_go_terms_p = 0
        self.average_number_initial_go_terms_c = 0
        self.average_number_final_go_terms_f = 0
        self.average_number_final_go_terms_p = 0
        self.average_number_final_go_terms_c = 0
        self.average_number_initial_do_terms = 0
        self.average_number_final_do_terms = 0
        self.average_number_go_annotations = 0
        self.average_number_do_annotations = 0
        self.average_number_orthologs = 0
        self.number_genes_with_non_null_orthology_description = 0
        self.number_genes_with_null_orthology_description = 0
        self.number_genes_with_more_than_3_best_orthologs = 0


class DescriptionsOverallProperties(object):
    def __init__(self, species: str = "", release_version: str = "", date: str = "", go_ontology_url: str = "",
                 go_association_url: str = "", do_ontology_url: str = "", do_association_url: str = ""):
        self.species = species
        self.release_version = release_version
        self.date = date
        self.go_ontology_url = go_ontology_url
        self.go_association_url = go_association_url
        self.do_ontology_url = do_ontology_url
        self.do_association_url = do_association_url


class SentenceMerger(object):
    def __init__(self):
        self.postfix_list = []
        self.terms_ids = set()
        self.term_postfix_dict = {}
        self.evidence_groups = []
        self.term_evgroup_dict = {}
        self.additional_prefix = ""
        self.aspect = ""
        self.qualifier = ""
        self.ancestors_covering_multiple_terms = set()


class SentenceGenerator(object):
    """generates sentences based on description rules"""
    def __init__(self, annotations: List[Dict], ontology: Ontology, evidence_groups_priority_list: List[str],
                 prepostfix_sentences_map: Dict[Tuple[str, str, str], Tuple[str, str]],
                 evidence_codes_groups_map: Dict[str, str],
                 prepostfix_special_cases_sent_map: Dict[Tuple[str, str, str], Tuple[int, str, str, str]] = None):
        """initialize sentence generator object

        Args:
            annotations (List[Dict]): the list of annotations for a given gene
            ontology (Ontology): the ontology linked to the annotations
            evidence_groups_priority_list (List[str]): the list of evidence groups to consider, sorted by priority.
                Sentences of the first group (with highest priority) will be returned in first position and so on
            prepostfix_sentences_map (Dict[Tuple[str, str, str], Tuple[str, str]]): a map with prefix and postfix
                phrases, where keys are tuples of go_aspect, evidence_group and values are tuples prefix, postfix
            prepostfix_special_cases_sent_map (Dict[Tuple[str, str, str], Tuple[int, str, str, str]]): a map for
                special prefix and postfix cases, where keys are tuples of aspect, evidence_group and values are tuples
                of id, match_regex, prefix, postfix. Match_regex is a regular expression that defines the match for the
                special case
            evidence_codes_groups_map (Dict[str, str]): a map between evidence codes and the groups they belong to
        """
        self.evidence_groups_priority_list = evidence_groups_priority_list
        self.prepostfix_sentences_map = prepostfix_sentences_map
        self.ontology = ontology
        self.terms_groups = defaultdict(lambda: defaultdict(set))

        if len(annotations) > 0:
            for annotation in annotations:
                if annotation["evidence"]["type"] in evidence_codes_groups_map:
                    aspect = annotation["aspect"]
                    ev_group = evidence_codes_groups_map[annotation["evidence"]["type"]]
                    qualifier = "_".join(sorted(annotation["qualifiers"])) if "qualifiers" in annotation else ""
                    if prepostfix_special_cases_sent_map and (aspect, ev_group, qualifier) in \
                            prepostfix_special_cases_sent_map:
                        for special_case in prepostfix_special_cases_sent_map[(aspect, ev_group, qualifier)]:
                            if re.match(re.escape(special_case[1]), ontology.label(annotation["object"]["id"],
                                                                                   id_if_null=True)):
                                ev_group = evidence_codes_groups_map[annotation["evidence"]["type"]] + \
                                           str(special_case[0])
                                if ev_group not in self.evidence_groups_priority_list:
                                    self.evidence_groups_priority_list.insert(evidence_groups_priority_list.index(
                                        evidence_codes_groups_map[annotation["evidence"]["type"]]) + 1, ev_group)
                                break
                    self.terms_groups[(aspect, qualifier)][ev_group].add(annotation["object"]["id"])

    def get_sentences(self, aspect: str, qualifier: str = '', keep_only_best_group: bool = False,
                      merge_groups_with_same_prefix: bool = False, remove_parent_terms: bool = True,
                      remove_child_terms: bool = False,
                      merge_num_terms_threshold: int = 3, merge_max_num_terms: int = 3,
                      merge_min_distance_from_root: dict = None,
                      truncate_others_generic_word: str = "several",
                      truncate_others_aspect_words: Dict[str, str] = None,
                      remove_successive_overlapped_terms: bool = True,
                      exclude_terms_ids: List[str] = None,
                      add_multiple_if_covers_more_children: bool = False, rename_cell: bool = False) -> List[Sentence]:
        """generate sentences for specific combination of aspect and qualifier

        Args:
            aspect (str): a data type aspect
            qualifier (str): qualifier
            keep_only_best_group (bool): whether to get only the evidence group with highest priority and discard
                the other evidence groups
            merge_groups_with_same_prefix (bool): whether to merge the phrases for evidence groups with the same prefix
            remove_parent_terms: whether to remove parent terms from the list of terms in each sentence if at least
                one children term is present
            remove_child_terms (bool): whether to remove child terms from the list of terms in each sentence if at least
                one parent term is present
            merge_num_terms_threshold (int): whether to merge terms by common ancestor to reduce the number of terms in
                the set. The trimming algorithm will be applied only if the number of terms is greater than the
                specified number and the specified threshold is greater than 0
            merge_max_num_terms (int): maximum number of terms to display in the final sentence
            merge_min_distance_from_root (dict): minimum distance from root terms for the selection of common ancestors
                during merging operations. Three values must be provided in the form of a dictionary with keys 'F', 'P',
                and 'C' for go aspect names and values integers indicating the threshold for each aspect
            truncate_others_generic_word (str): a generic word to indicate that the set of terms reported in the
                sentence is only a subset of the original terms, e.g., 'several'
            truncate_others_aspect_words (Dict[str, str]): one word for each aspect describing the kind of terms that
                are included in the aspect
            remove_successive_overlapped_terms (bool): whether to remove terms in lower priority evidence groups when
                already present in higher priority groups
            exclude_terms_ids (List[str]): list of term ids to exclude
            add_multiple_if_covers_more_children (bool): whether to add the label '(multiple)' to terms that are
                ancestors covering multiple children
            rename_cell (bool): rename term cell if present
        Returns:
            List[Sentence]: a list of sentences
        """
        if not merge_min_distance_from_root:
            merge_min_distance_from_root = {'F': 1, 'P': 1, 'C': 2, 'D': 3}
        if not truncate_others_aspect_words:
            truncate_others_aspect_words = {'F': 'functions', 'P': 'processes', 'C': 'components', 'D': 'diseases'}
        sentences = []
        terms_already_covered = set()
        evidence_group_priority = {eg: p for p, eg in enumerate(self.evidence_groups_priority_list)}
        for terms, evidence_group, priority in sorted([(t, eg, evidence_group_priority[eg]) for eg, t in
                                                       self.terms_groups[(aspect, qualifier)].items()],
                                                      key=lambda x: x[2]):
            ancestors_covering_multiple_children = set()
            if remove_successive_overlapped_terms:
                terms -= terms_already_covered
            if exclude_terms_ids:
                terms -= set(exclude_terms_ids)
            add_others = False
            if remove_parent_terms:
                terms_no_ancestors = terms - set([ancestor for node_id in terms for ancestor in
                                                  self.ontology.ancestors(node_id)])
                if len(terms) > len(terms_no_ancestors):
                    logging.debug("Removed " + str(len(terms) - len(terms_no_ancestors)) + " parents from terms")
                    terms = terms_no_ancestors
            if 0 < merge_num_terms_threshold <= len(terms):
                merged_terms_coverset = get_merged_nodes_by_common_ancestor(
                    node_ids=list(terms), ontology=self.ontology,
                    min_distance_from_root=merge_min_distance_from_root[aspect],
                    min_number_of_terms=merge_num_terms_threshold)
                if len(merged_terms_coverset.keys()) <= merge_num_terms_threshold:
                    merged_terms = list(merged_terms_coverset.keys())
                    terms_already_covered.update([e for subset in merged_terms_coverset.values() for e in subset])
                else:
                    merged_terms = find_set_covering([(k, self.ontology.label(k, id_if_null=True), v) for k, v in
                                                      merged_terms_coverset.items()],
                                                     max_num_subsets=merge_max_num_terms)
                    for merged_term in merged_terms:
                        terms_already_covered.update(merged_terms_coverset[merged_term])
                    add_others = True
                if add_multiple_if_covers_more_children:
                    ancestors_covering_multiple_children = {self.ontology.label(ancestor, id_if_null=True) for ancestor
                                                            in merged_terms if len(merged_terms_coverset[ancestor]) > 1}
                logging.debug("Reduced number of terms by merging from " + str(len(terms)) + " to " +
                              str(len(merged_terms)))
                terms = merged_terms
            else:
                terms_already_covered.update(terms)
            if remove_child_terms:
                terms = [term for term in terms if
                         len(set(self.ontology.ancestors(term)).intersection(set(terms))) == 0]
            if (aspect, evidence_group, qualifier) in self.prepostfix_sentences_map:
                sentences.append(
                    _get_single_sentence(node_ids=terms, ontology=self.ontology, aspect=aspect,
                                         evidence_group=evidence_group, qualifier=qualifier,
                                         prepostfix_sentences_map=self.prepostfix_sentences_map,
                                         terms_merged=True if 0 < merge_num_terms_threshold < len(terms) else False,
                                         add_others=add_others,
                                         truncate_others_generic_word=truncate_others_generic_word,
                                         truncate_others_aspect_words=truncate_others_aspect_words,
                                         ancestors_with_multiple_children=ancestors_covering_multiple_children,
                                         rename_cell=rename_cell))
                if keep_only_best_group:
                    return sentences
        if merge_groups_with_same_prefix:
            sentences = self.merge_sentences_with_same_prefix(sentences=sentences,
                                                              remove_parent_terms=remove_parent_terms,
                                                              rename_cell=rename_cell)
        return sentences

    def merge_sentences_with_same_prefix(self, sentences: List[Sentence], remove_parent_terms: bool = True,
                                         rename_cell: bool = False):
        """merge sentences with the same prefix

        Args:
            sentences (List[Sentence]): a list of sentences
            remove_parent_terms (bool): whether to remove parent terms if present in the merged set of terms
            rename_cell (bool): whether to rename the term 'cell'
        Returns:
            List[Sentence]: the list of merged sentences, sorted by (merged) evidence group priority
        """
        merged_sentences = defaultdict(SentenceMerger)
        for sentence in sentences:
            prefix = self.prepostfix_sentences_map[(sentence.aspect, sentence.evidence_group, sentence.qualifier)][0]
            merged_sentences[prefix].postfix_list.append(self.prepostfix_sentences_map[(sentence.aspect,
                                                                                        sentence.evidence_group,
                                                                                        sentence.qualifier)][1])
            merged_sentences[prefix].aspect = sentence.aspect
            merged_sentences[prefix].qualifier = sentence.qualifier
            merged_sentences[prefix].terms_ids.update(sentence.terms_ids)
            for term in sentence.terms_ids:
                merged_sentences[prefix].term_postfix_dict[term] = self.prepostfix_sentences_map[
                    (sentence.aspect, sentence.evidence_group, sentence.qualifier)][1]
            merged_sentences[prefix].evidence_groups.append(sentence.evidence_group)
            for term in sentence.terms_ids:
                merged_sentences[prefix].term_evgroup_dict[term] = sentence.evidence_group
            if sentence.additional_prefix:
                merged_sentences[prefix].additional_prefix = sentence.additional_prefix
            merged_sentences[prefix].ancestors_covering_multiple_terms.update(
                sentence.ancestors_covering_multiple_terms)
        if remove_parent_terms:
            for prefix, sent_merger in merged_sentences.items():
                terms_no_ancestors = sent_merger.terms_ids - set([ancestor for node_id in sent_merger.terms_ids for
                                                                  ancestor in self.ontology.ancestors(node_id)])
                if len(sent_merger.terms_ids) > len(terms_no_ancestors):
                    logging.debug("Removed " + str(len(sent_merger.terms_ids) - len(terms_no_ancestors)) +
                                  " parents from terms while merging sentences with same prefix")
                    sent_merger.terms_ids = terms_no_ancestors
        return [Sentence(prefix=prefix, terms_ids=list(sent_merger.terms_ids),
                         postfix=SentenceGenerator.merge_postfix_phrases(sent_merger.postfix_list),
                         text=compose_sentence(prefix=prefix,
                         term_names=[self.ontology.label(node, id_if_null=True) for node in sent_merger.terms_ids],
                         postfix=SentenceGenerator.merge_postfix_phrases(sent_merger.postfix_list),
                         additional_prefix=sent_merger.additional_prefix,
                         ancestors_with_multiple_children=sent_merger.ancestors_covering_multiple_terms,
                         rename_cell=rename_cell),
                         aspect=sent_merger.aspect, evidence_group=", ".join(sent_merger.evidence_groups),
                         terms_merged=True, additional_prefix=sent_merger.additional_prefix,
                         qualifier=sent_merger.qualifier,
                         ancestors_covering_multiple_terms=sent_merger.ancestors_covering_multiple_terms)
                for prefix, sent_merger in merged_sentences.items() if len(sent_merger.terms_ids) > 0]

    @staticmethod
    def merge_postfix_phrases(postfix_phrases: List[str]) -> str:
        """merge postfix phrases and remove possible redundant text at the beginning at at the end of the phrases

        Args:
            postfix_phrases (List[str]): the phrases to merge
        Returns:
            str: the merged postfix phrase
        """
        postfix_phrases = [postfix for postfix in postfix_phrases if postfix]
        if postfix_phrases and len(postfix_phrases) > 0:
            if len(postfix_phrases) > 1:
                inf_engine = inflect.engine()
                shortest_phrase = sorted(zip(postfix_phrases, [len(phrase) for phrase in postfix_phrases]),
                                         key=lambda x: x[1])[0][0]
                first_part = ""
                for idx, letter in enumerate(shortest_phrase):
                    if all(map(lambda x: x[idx] == shortest_phrase[idx], postfix_phrases)):
                        first_part += letter
                    else:
                        break
                last_part = ""
                for idx, letter in zip(range(len(shortest_phrase)), reversed([l for l in shortest_phrase])):
                    if all(map(lambda x: x[len(x) - idx - 1] == shortest_phrase[len(shortest_phrase) - idx - 1],
                               postfix_phrases)):
                        last_part = letter + last_part
                    else:
                        break
                new_phrases = [phrase.replace(first_part, "").replace(last_part, "") for phrase in postfix_phrases]
                if len(last_part.strip().split(" ")) == 1:
                    last_part = inf_engine.plural(last_part)
                if len(new_phrases) > 2:
                    return first_part + ", ".join(new_phrases[0:-1]) + ", and " + new_phrases[-1] + last_part
                elif len(new_phrases) > 1:
                    return first_part + " and ".join(new_phrases) + last_part
                else:
                    return first_part + new_phrases[0] + last_part
            else:
                return postfix_phrases[0]
        else:
            return ""


def compose_sentence(prefix: str, additional_prefix: str, term_names: List[str], postfix: str,
                     ancestors_with_multiple_children: Set[str] = None, rename_cell: bool = False) -> str:
    """compose the text of a sentence given its prefix, terms, and postfix

    Args:
        prefix (str): the prefix of the sentence
        term_names (List[str]): a list of term names
        postfix (str): the postfix of the sentence
        additional_prefix (str): an additional prefix to be used for special cases
        ancestors_with_multiple_children (Set[str]): set containing labels of terms that cover more than one children
            term in the original set and which will appear with the label '(multiple)'
    Returns:
        str: the text of the go sentence
    """
    prefix = prefix + additional_prefix + " "
    term_names = [term_name + " (multiple)" if term_name in ancestors_with_multiple_children else term_name for
                  term_name in sorted(term_names)]

    if postfix != "":
        postfix = " " + postfix
    if rename_cell:
        if "the cell" in term_names or "the Cell" in term_names:
            if len(term_names) == 1:
                prefix = prefix[0:-3]
                term_names = ["widely"]
            else:
                if not additional_prefix:
                    prefix += "several tissues including "
                term_names = [term for term in term_names if term != "the cell" and term != "the Cell"]
    return prefix + concatenate_words_with_oxford_comma(term_names) + postfix


def _get_single_sentence(node_ids: List[str], ontology: Ontology, aspect: str, evidence_group: str, qualifier: str,
                         prepostfix_sentences_map: Dict[Tuple[str, str, str], Tuple[str, str]],
                         terms_merged: bool = False, add_others: bool = False,
                         truncate_others_generic_word: str = "several",
                         truncate_others_aspect_words: Dict[str, str] = None,
                         ancestors_with_multiple_children: Set[str] = None,
                         rename_cell: bool = False) -> Union[Sentence, None]:
    """build a sentence object

    Args:
        node_ids (List[str]): list of ids for the terms to be combined in the sentence
        ontology (Ontology): the ontology containing the nodes
        aspect (str): aspect
        evidence_group (str): evidence group
        qualifier (str): qualifier
        prepostfix_sentences_map (Dict[Tuple[str, str, str], Tuple[str, str]]): map for prefix and postfix phrases
        terms_merged (bool): whether the terms set has been merged to reduce its size
        add_others (bool): whether to say that there are other terms which have been omitted from the sentence
        truncate_others_generic_word (str): a generic word to indicate that the set of terms reported in the sentence is
            only a subset of the original terms, e.g., 'several'
        truncate_others_aspect_words (Dict[str, str]): one word for each aspect describing the kind of terms that are
            included in the aspect
        ancestors_with_multiple_children (Set[str]): set containing labels of terms that cover more than one children
            term in the original set and which will appear with the label '(multiple)'
        rename_cell (bool): whether to rename the term 'cell'
    Returns:
        Union[Sentence,None]: the combined go sentence
    """
    if len(node_ids) > 0:
        prefix = prepostfix_sentences_map[(aspect, evidence_group, qualifier)][0]
        additional_prefix = ""
        others_word = "entities"
        if aspect in truncate_others_aspect_words:
            others_word = truncate_others_aspect_words[aspect]
        if add_others:
            additional_prefix += " " + truncate_others_generic_word + " " + others_word + ", including"
        if aspect == "C":
            additional_prefix += " the"
        postfix = prepostfix_sentences_map[(aspect, evidence_group, qualifier)][1]
        term_labels = [ontology.label(node_id, id_if_null=True) for node_id in node_ids]
        return Sentence(prefix=prefix, terms_ids=node_ids, postfix=postfix,
                        text=compose_sentence(prefix=prefix, term_names=term_labels, postfix=postfix,
                                              additional_prefix=additional_prefix,
                                              ancestors_with_multiple_children=ancestors_with_multiple_children,
                                              rename_cell=rename_cell),
                        aspect=aspect, evidence_group=evidence_group, terms_merged=terms_merged,
                        additional_prefix=additional_prefix, qualifier=qualifier,
                        ancestors_covering_multiple_terms=ancestors_with_multiple_children)
    else:
        return None


def _generate_ortholog_sentence_wormbase_human(orthologs: List[List[str]], human_genes_props: Dict[str, List[str]]):
    """build orthology sentence for WormBase human orthologs

    Args:
        orthologs (List[List[str]]): list of human orthologs, containing gene_id, gene_symbol
        human_genes_props (Dict[str, List[str]]): dictionary containing human gene properties
    Returns:
        Tuple[list, str]: the orthologs and the sentence
    """
    symbol_name_arr = []
    genes_in_families = []
    gene_symbols_wo_family = []
    if len(orthologs) > 3:
        gene_families = defaultdict(list)
        gene_symbols_wo_family = set()
        for ortholog in orthologs:
            if ortholog[0] in human_genes_props and human_genes_props[ortholog[0]] and \
                    len(human_genes_props[ortholog[0]]) == 4:
                gene_families[human_genes_props[ortholog[0]][2]].append(human_genes_props[ortholog[0]])
            else:
                gene_symbols_wo_family.add(ortholog[1])
        if len(list(gene_families.keys())) == 1:
            gene_symbols_wo_family.update(set([human_p[0] + " (" + human_p[1] + ")" for human_p in
                                               gene_families[list(gene_families.keys())[0]]]))
            gene_families = {}
        else:
            for family_symbol, human_ps in gene_families.items():
                if family_symbol == "" or len(human_ps) == 1:
                    for human_p in human_ps:
                        gene_symbols_wo_family.add(human_p[0] + " (" + human_p[1] + ")")
            gene_families = {family_name: human_ps for family_name, human_ps in gene_families.items() if
                             len(human_ps) > 1 and family_name != ""}
        gene_family_names = [human_ps[0][2] + " (" + human_ps[0][3] + ")" for human_ps in gene_families.values()]
        genes_in_families = list(set([hps[0] for gene_list in gene_families.values() for hps in gene_list]))
        gene_symbols_wo_family = list(gene_symbols_wo_family)
        if len(gene_family_names) > 3:
            gene_family_names = gene_family_names[0:3]
        if len(genes_in_families) > 3:
            genes_in_families = genes_in_families[0:3]
        if len(gene_symbols_wo_family) > 3:
            gene_symbols_wo_family = gene_symbols_wo_family[0:3]
        family_word = "family"
        if len(gene_family_names) > 1:
            family_word = "families"
        sentences_arr = []
        if len(gene_symbols_wo_family) > 0:
            sentences_arr.append("human " + concatenate_words_with_oxford_comma(gene_symbols_wo_family))
        if len(gene_family_names) > 0:
            sentences_arr.append("members of the human " + concatenate_words_with_oxford_comma(gene_family_names) +
                                 " gene " + family_word + " including " + concatenate_words_with_oxford_comma(
                genes_in_families))
        orth_sentence = "is an ortholog of " + " and ".join(sentences_arr)
    else:
        symbol_name_arr = sorted([human_genes_props[best_orth[0]][0] + " (" + human_genes_props[best_orth[0]][1] +
                                  ")" if best_orth[0] in human_genes_props and human_genes_props[best_orth[0]] else
                                  best_orth[1] for best_orth in orthologs])
        orth_sentence = "is an ortholog of human " + concatenate_words_with_oxford_comma(symbol_name_arr)

    return [gene.split(" ")[0] for gene in [*genes_in_families, *gene_symbols_wo_family, *symbol_name_arr]], \
           orth_sentence


def generate_orthology_sentence_alliance_human(orthologs: List[List[str]]):
    """build orthology sentence for Alliance human orthologs

    Args:
        orthologs (List[List[str]]): list of human orthologs, containing gene_id, gene_symbol, and gene_name
    Returns:
        str: the orthology sentence
    """
    if len(orthologs) > 0:
        prefix = "human"
        orthologs_display = sorted(orthologs, key=lambda x: x[2])
        if len(orthologs) > 3:
            orthologs_display = orthologs_display[0:3]
            prefix = "several human genes including"
        return "orthologous to " + prefix + " " + concatenate_words_with_oxford_comma(
            [orth[1] + " (" + orth[2] + ")" if orth[2] else orth[1] for orth in orthologs_display])
    else:
        return None


def _generate_ortholog_sentence_wormbase_non_c_elegans(orthologs: List[List[str]], orthologs_sp_fullname: str,
                                                       textpresso_api_token: str):
    """build orthology sentence for WormBase non-human hortologs

        Args:
            orthologs (List[str]): list of human orthologs, containing gene_id, gene_symbol
            orthologs_sp_fullname (str): full name of species from which to extract orthologs
            textpresso_api_token (str): token to access Textpresso Central API
        Returns:
            str: the orthology sentence
        """
    orth_sentence = None
    if len(orthologs) > 0:
        fullname_arr = orthologs_sp_fullname.split(" ")
        if len(fullname_arr[0]) > 2:
            fullname_arr[0] = fullname_arr[0][0] + "."
            orthologs_sp_fullname = " ".join(fullname_arr)
        if len(orthologs) > 3:
            # sort orthologs by tpc popularity and alphabetically (if tied)
            orthologs_pop = [o_p for o_p in sorted([[ortholog, get_textpresso_popularity(
                    textpresso_api_token, ortholog[1])] for ortholog in orthologs], key=lambda x: (x[1], x[0][1]),
                                                  reverse=True)]
            classes_orth_pop = defaultdict(list)
            orthologs_pop_wo_class = []
            for o_p in orthologs_pop:
                gene_class = get_gene_class(o_p[0][0])
                if gene_class:
                    classes_orth_pop[gene_class].append(o_p)
                else:
                    orthologs_pop_wo_class.append(o_p)
            if len(list(classes_orth_pop.keys())) == 1:
                orthologs_pop_wo_class.extend(classes_orth_pop[list(classes_orth_pop.keys())[0]])
                classes_orth_pop = {}
            else:
                for gene_class, orths_with_pop in classes_orth_pop.items():
                    if len(orths_with_pop) == 1:
                        orthologs_pop_wo_class.extend(orths_with_pop)
            classes_orth_pop = {gene_class: ops[0] for gene_class, ops in classes_orth_pop.items() if len(ops) > 1}
            sorted_items = [[o_p, 0] for o_p in orthologs_pop_wo_class]
            sorted_items.extend([[o_p, 1, gene_class] for gene_class, o_p in classes_orth_pop.items()])
            sorted_items.sort(key=lambda x: x[0][1], reverse=True)
            if len(sorted_items) > 3:
                sorted_items = sorted_items[0:3]
            gene_symbols_wo_class = [item[0][0][1] for item in sorted_items if item[1] == 0]
            classes_symbols = [item[2] for item in sorted_items if item[1] == 1]
            genes_symbols_in_classes = [item[0][0][1] for item in sorted_items if item[1] == 1]
            sentences_arr = []
            if len(gene_symbols_wo_class) > 0:
                sentences_arr.append(orthologs_sp_fullname + " " + concatenate_words_with_oxford_comma(
                    gene_symbols_wo_class))
            if len(classes_symbols) > 0:
                genes_symbols_in_classes_sent = concatenate_words_with_oxford_comma(genes_symbols_in_classes)
                classes_symbols_sent = concatenate_words_with_oxford_comma(classes_symbols)
                classes_word = "classes" if len(classes_symbols) > 1 else "class"
                sentences_arr.append("members of the " + orthologs_sp_fullname + " " + classes_symbols_sent +
                                     " gene " + classes_word + " including " + genes_symbols_in_classes_sent)
            orth_sentence = "is an ortholog of " + " and ".join(sentences_arr)
        else:
            # sort orthologs alphabetically
            orthologs_symbols = sorted([orth[1] for orth in orthologs])
            orth_sentence = "is an ortholog of " + orthologs_sp_fullname + " " + \
                            concatenate_words_with_oxford_comma(orthologs_symbols)
    return orth_sentence


def rename_human_ortholog_name(ortholog_name: str):
    new_name = ortholog_name
    if " family member " in ortholog_name:
        new_name = new_name.replace(" family member ", " ")
    new_name = re.sub(r"[,]? kDa$", "", new_name)
    return new_name


def is_human_ortholog_name_valid(ortholog_name: str):
    if "human uncharacterized protein" in ortholog_name.lower():
        return False
    return True


def get_gene_class(gene_id: str):
    """get the gene class of a gene from WormBase API

    Args:
        gene_id (str): the Wormbase WBGene ID of the gene
    Returns:
        str: the class of the gene
    """
    logging.debug("Getting gene class for gene " + gene_id)
    try:
        gene_class_data = json.loads(urllib.request.urlopen("http://rest.wormbase.org/rest/field/gene/" + gene_id +
                                                            "/gene_class").read())
        if "gene_class" in gene_class_data and gene_class_data["gene_class"]["data"] and "tag" in \
                gene_class_data["gene_class"]["data"] and "label" in gene_class_data["gene_class"]["data"]["tag"]:
                return gene_class_data["gene_class"]["data"]["tag"]["label"]
    except:
        return None
    return None


def get_textpresso_popularity(textpresso_api_token: str, keywords: str):
    """get the number of papers in the C. elegans literature that mention a certain keyword from Textpresso Central API

    Args:
        textpresso_api_token (str): a valid token to access Textpresso Central API
        keywords (str): the keyword to search, or any combination of keywords containing AND and OR operators
    Returns:
        int: the popularity of the specified keyword
    """
    if not os.environ.get('PYTHONHTTPSVERIFY', '') and getattr(ssl, '_create_unverified_context', None):
        ssl._create_default_https_context = ssl._create_unverified_context
    api_endpoint = "https://textpressocentral.org:18080/v1/textpresso/api/get_documents_count"
    data = json.dumps({"token": textpresso_api_token, "query": {
        "keywords": keywords, "type": "document", "corpora": ["C. elegans"]}})
    data = data.encode('utf-8')
    req = urllib.request.Request(api_endpoint, data, headers={'Content-type': 'application/json',
                                                              'Accept': 'application/json'})
    res = urllib.request.urlopen(req)
    return int(json.loads(res.read().decode('utf-8')))


def concatenate_words_with_oxford_comma(words: List[str]):
    """concatenate words by separating them with commas and a final oxford comma if more than two or by 'and' if two

    Args:
        words (List[str]): a list of words

    Returns:
        str: a concatenated string representing the list of words
    """
    if len(words) > 2:
        return ", ".join(words[0:-1]) + ", and " + words[-1]
    else:
        return " and ".join(words)


def get_best_human_ortholog_for_info_poor(human_orthologs, ensembl_hgnc_ids_map, evidence_codes, human_df_agr,
                                          go_sent_gen_common_props):
    best_orth = ""
    if len(human_orthologs) > 0:
        exp_orthologs = defaultdict(int)
        predicted_orthologs = defaultdict(int)
        for ortholog in human_orthologs:
            if ortholog[0] in ensembl_hgnc_ids_map and ensembl_hgnc_ids_map[ortholog[0]]:
                ortholog_annotations = human_df_agr.get_annotations_for_gene(
                    gene_id="RGD:" + ensembl_hgnc_ids_map[ortholog[0]], annot_type=DataType.GO,
                    priority_list=evidence_codes)
                for annotation in ortholog_annotations:
                    if annotation['aspect'] == 'F':
                        if go_sent_gen_common_props["evidence_codes_groups_map"][annotation["evidence"]['type']] \
                                in ["EXPERIMENTAL", "HIGH_THROUGHPUT_EXPERIMENTAL"]:
                            exp_orthologs[ensembl_hgnc_ids_map[ortholog[0]]] += 1
                        else:
                            predicted_orthologs[ensembl_hgnc_ids_map[ortholog[0]]] += 1
        if len(exp_orthologs) > 0:
            best_orth = sorted([(key, value) for key, value in exp_orthologs.items()], key=lambda x: x[1],
                               reverse=True)[0][0]
        elif len(predicted_orthologs) > 0:
            best_orth = sorted([(key, value) for key, value in predicted_orthologs.items()], key=lambda x: x[1],
                               reverse=True)[0][0]
    return best_orth


def compose_wormbase_description(gene: Gene, conf_parser: GenedescConfigParser, species, organism, df,
                                 orthologs_sp_fullname, go_sent_gen_common_props, go_sent_common_props,
                                 human_genes_props, do_sent_gen_common_prop, do_sent_common_props, sister_sp_fullname,
                                 sister_df, human_df_agr, desc_writer, ensembl_hgnc_ids_map,
                                 expr_sent_gen_common_props, expr_sent_common_props):
    """compose gene descriptions for WormBase

    Args:
        gene (Gene): a gene object
        conf_parser (GenedescConfigParser): a configuration parser
        species: the species to process
        organism: (str): the organism to process
        df (DataFetcher): the data fetcher containing the data for the description
        orthologs_sp_fullname (str): full name of the organism for orthology
        go_sent_gen_common_props (dict): common properties for go sentences generator
        go_sent_common_props (dict): common properties for go sentences
        human_genes_props (dict): human gene properties
        do_sent_gen_common_prop (dict]): common properties for do sentences generator
        do_sent_common_props (dict): common properties for do sentences
        sister_sp_fullname (str): full name of sister species
        sister_df (DataFetcher): sister species data fetcher
        human_df_agr (DataFetcher): data fetcher to generate GO sentences for information poor genes from AGR human data
        desc_writer (DescriptionWriter): description writer
        ensembl_hgnc_ids_map (Dict[str, str]): ensembl hgnc map
        expr_sent_gen_common_props (dict): common properties for expression sentences generator
        expr_sent_common_props (dict) common properties for expression sentences
    """
    gene_desc = GeneDesc(gene_id=gene.id, gene_name=gene.name,
                         publications=", ".join([annot["publication"] for annot in df.get_annotations_for_gene(
                             gene.id, annot_type=DataType.GO,
                             priority_list=conf_parser.get_go_evidence_groups_priority_list())]),
                         refs=", ".join([annot["refs"] for annot in df.get_annotations_for_gene(
                             gene.id, annot_type=DataType.GO,
                             priority_list=conf_parser.get_go_evidence_groups_priority_list())]))
    joined_sent = []

    best_orthologs, selected_orth_name = df.get_best_orthologs_for_gene(
        gene.id, orth_species_full_name=orthologs_sp_fullname)
    selected_orthologs = []
    if best_orthologs:
        gene_desc.stats.set_best_orthologs = [orth[0] for orth in best_orthologs]
        if len(orthologs_sp_fullname) == 1 and orthologs_sp_fullname[0] == "Homo sapiens":
            sel_orthologs, orth_sent = _generate_ortholog_sentence_wormbase_human(best_orthologs, human_genes_props)
            selected_orthologs = [orth for orth in best_orthologs if orth[1] in sel_orthologs]
        else:
            orth_sent = _generate_ortholog_sentence_wormbase_non_c_elegans(best_orthologs, selected_orth_name,
                                                                           conf_parser.get_textpresso_api_token())
        if orth_sent:
            joined_sent.append(orth_sent)
            gene_desc.orthology_description = orth_sent
    go_annotations = df.get_annotations_for_gene(gene_id=gene.id, annot_type=DataType.GO,
                                                 priority_list=conf_parser.get_go_annotations_priority())
    go_sent_gen_common_props_exp = go_sent_gen_common_props.copy()
    go_sent_gen_common_props_exp["evidence_codes_groups_map"] = {
        evcode: group for evcode, group in go_sent_gen_common_props["evidence_codes_groups_map"].items() if
        "EXPERIMENTAL" in go_sent_gen_common_props["evidence_codes_groups_map"][evcode]}
    go_sent_generator_exp = SentenceGenerator(annotations=go_annotations, ontology=df.go_ontology,
                                              **go_sent_gen_common_props_exp)
    go_sent_generator = SentenceGenerator(annotations=go_annotations, ontology=df.go_ontology,
                                          **go_sent_gen_common_props)
    gene_desc.stats.total_number_go_annotations = len(go_annotations)
    gene_desc.stats.set_initial_go_ids_f = list(set().union(
        [elem for key, sets in go_sent_generator.terms_groups[('F', '')].items() for elem in sets if ('F', key, '') in
         conf_parser.get_go_prepostfix_sentences_map()], [elem for key, sets in go_sent_generator.terms_groups[
            ('F', 'contributes_to')].items() for elem in sets if ('F', key, 'contributes_to') in
                                                          conf_parser.get_go_prepostfix_sentences_map()]))
    gene_desc.stats.set_initial_go_ids_p = [elem for key, sets in go_sent_generator.terms_groups[('P', '')].items() for
                                            elem in sets if ('P', key, '') in
                                            conf_parser.get_go_prepostfix_sentences_map()]
    gene_desc.stats.set_initial_go_ids_c = list(set().union(
        [elem for key, sets in go_sent_generator.terms_groups[('C', '')].items() for elem in sets if ('C', key, '') in
         conf_parser.get_go_prepostfix_sentences_map()],
        [elem for key, sets in go_sent_generator.terms_groups[('C', 'colocalizes_with')].items() for elem in sets if
         ('C', key, 'colocalizes_with') in conf_parser.get_go_prepostfix_sentences_map()]))
    contributes_to_raw_func_sent = go_sent_generator.get_sentences(
        aspect='F', qualifier='contributes_to', merge_groups_with_same_prefix=True, keep_only_best_group=True,
        **go_sent_common_props)
    if contributes_to_raw_func_sent:
        raw_func_sent = go_sent_generator_exp.get_sentences(aspect='F', merge_groups_with_same_prefix=True,
                                                            keep_only_best_group=True, **go_sent_common_props)
    else:
        raw_func_sent = go_sent_generator.get_sentences(aspect='F', merge_groups_with_same_prefix=True,
                                                        keep_only_best_group=True, **go_sent_common_props)
    func_sent = " and ".join([sentence.text for sentence in raw_func_sent])
    if func_sent:
        joined_sent.append(func_sent)
        gene_desc.go_function_description = func_sent
        gene_desc.go_description = func_sent
    gene_desc.stats.set_final_go_ids_f = list(set().union([term_id for sentence in raw_func_sent for
                                                           term_id in sentence.terms_ids],
                                                          [term_id for sentence in contributes_to_raw_func_sent for
                                                           term_id in sentence.terms_ids]))
    contributes_to_func_sent = " and ".join([sentence.text for sentence in contributes_to_raw_func_sent])
    if contributes_to_func_sent:
        joined_sent.append(contributes_to_func_sent)
        if not gene_desc.go_function_description:
            gene_desc.go_function_description = contributes_to_func_sent
        else:
            gene_desc.go_function_description += "; " + contributes_to_func_sent
        if not gene_desc.go_description:
            gene_desc.go_description = contributes_to_func_sent
        else:
            gene_desc.go_description += "; " + contributes_to_func_sent
    raw_proc_sent = go_sent_generator.get_sentences(aspect='P', merge_groups_with_same_prefix=True,
                                                    keep_only_best_group=True, **go_sent_common_props)
    gene_desc.stats.set_final_go_ids_p = [term_id for sentence in raw_proc_sent for term_id in sentence.terms_ids]
    proc_sent = " and ".join([sentence.text for sentence in raw_proc_sent])
    if proc_sent:
        joined_sent.append(proc_sent)
        gene_desc.go_process_description = proc_sent
        if not gene_desc.go_description:
            gene_desc.go_description = proc_sent
        else:
            gene_desc.go_description += "; " + proc_sent
    colocalizes_with_raw_comp_sent = go_sent_generator.get_sentences(
        aspect='C', qualifier='colocalizes_with', merge_groups_with_same_prefix=True,
        keep_only_best_group=True, **go_sent_common_props)
    if colocalizes_with_raw_comp_sent:
        raw_comp_sent = go_sent_generator_exp.get_sentences(aspect='C', merge_groups_with_same_prefix=True,
                                                            keep_only_best_group=True, **go_sent_common_props)
    else:
        raw_comp_sent = go_sent_generator.get_sentences(aspect='C', merge_groups_with_same_prefix=True,
                                                        keep_only_best_group=True, **go_sent_common_props)
    comp_sent = " and ".join([sentence.text for sentence in raw_comp_sent])
    if comp_sent:
        joined_sent.append(comp_sent)
        gene_desc.go_component_description = comp_sent
        if not gene_desc.go_description:
            gene_desc.go_description = comp_sent
        else:
            gene_desc.go_description += " " + comp_sent

    gene_desc.stats.set_final_go_ids_c = list(set().union([term_id for sentence in raw_comp_sent for
                                                           term_id in sentence.terms_ids],
                                                          [term_id for sentence in colocalizes_with_raw_comp_sent for
                                                           term_id in sentence.terms_ids]))
    colocalizes_with_comp_sent = " and ".join([sentence.text for sentence in colocalizes_with_raw_comp_sent])
    if colocalizes_with_comp_sent:
        joined_sent.append(colocalizes_with_comp_sent)
        if not gene_desc.go_component_description:
            gene_desc.go_component_description = colocalizes_with_comp_sent
        else:
            gene_desc.go_component_description += "; " + colocalizes_with_comp_sent
        if not gene_desc.go_description:
            gene_desc.go_description = colocalizes_with_comp_sent
        else:
            gene_desc.go_description += "; " + colocalizes_with_comp_sent
    if conf_parser.get_data_fetcher() == "wb_data_fetcher" and "main_sister_species" in species[organism] and \
            species[organism]["main_sister_species"] and df.get_best_orthologs_for_gene(
        gene.id, orth_species_full_name=[sister_sp_fullname], sister_species_data_fetcher=sister_df,
        ecode_priority_list=["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI",
                             "HEP"])[0]:
        best_ortholog = df.get_best_orthologs_for_gene(
            gene.id, orth_species_full_name=[sister_sp_fullname], sister_species_data_fetcher=sister_df,
            ecode_priority_list=["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI",
                                 "HEP"])[0][0]
        sister_sentences_generator = SentenceGenerator(sister_df.get_annotations_for_gene(
            annot_type=DataType.GO, gene_id="WB:" + best_ortholog[0],
            priority_list=("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP")),
            ontology=df.go_ontology, **go_sent_gen_common_props)
        sister_proc_sent = " and ".join([sentence.text for sentence in sister_sentences_generator.get_sentences(
            aspect='P', merge_groups_with_same_prefix=True, keep_only_best_group=True, **go_sent_common_props)])
        if sister_proc_sent:
            joined_sent.append("in " + species[species[organism]["main_sister_species"]]["name"] + ", " +
                               best_ortholog[1] + " " + sister_proc_sent)
    expr_annotations = df.get_annotations_for_gene(gene_id=gene.id, annot_type=DataType.EXPR,
                                                   priority_list=conf_parser.get_expression_annotations_priority())
    expr_sentence_generator = SentenceGenerator(annotations=expr_annotations, ontology=df.expression_ontology,
                                                **expr_sent_gen_common_props)
    raw_expression_sent = expr_sentence_generator.get_sentences(
        aspect='A', qualifier="Verified", merge_groups_with_same_prefix=True, keep_only_best_group=False,
        **expr_sent_common_props)
    expression_sent = "; ".join([sentence.text for sentence in raw_expression_sent])
    if expression_sent:
        postfix = ""
        if len(df.expression_enriched_extra_data[gene.id[3:]]) > 0:
            postfix = " " + concatenate_words_with_oxford_comma(df.expression_enriched_extra_data[gene.id[3:]]) + \
                      " studies"
        joined_sent.append(expression_sent + postfix)
    if len(joined_sent) == 0:
        raw_expression_sent_enriched = expr_sentence_generator.get_sentences(
            aspect='A', qualifier="Enriched", merge_groups_with_same_prefix=True, keep_only_best_group=False,
            **expr_sent_common_props)
        expression_sent_enriched = ""
        postfix = ""
        if df.expression_ontology is not None:
            expression_sent_enriched = "; ".join([sentence.text for sentence in raw_expression_sent_enriched])
            postfix = " " + concatenate_words_with_oxford_comma(
                df.expression_enriched_extra_data[gene.id[3:]]) + " studies"
        elif df.expression_enriched_bma_data[gene.id[3:]] and len(df.expression_enriched_bma_data[gene.id[3:]][3]) > 0:
            expression_sent_enriched = "is enriched in " + concatenate_words_with_oxford_comma(
                df.expression_enriched_bma_data[gene.id[3:]][2])
            postfix = " based on " + concatenate_words_with_oxford_comma(
                df.expression_enriched_bma_data[gene.id[3:]][3]) + " studies"
        if expression_sent_enriched:
            joined_sent.append(expression_sent_enriched + postfix)
        if df.expression_ontology is None and df.expression_affected_bma_data[gene.id[3:]] and \
                len(df.expression_affected_bma_data[gene.id[3:]][3]) > 0:
            expression_sent_affected = "is affected by " + concatenate_words_with_oxford_comma(
                df.expression_affected_bma_data[gene.id[3:]][2]) + " based on " + \
                                       concatenate_words_with_oxford_comma(
                                           df.expression_affected_bma_data[gene.id[3:]][3]) + " studies"
            joined_sent.append(expression_sent_affected)
    do_annotations = df.get_annotations_for_gene(gene_id=gene.id, annot_type=DataType.DO,
                                                 priority_list=conf_parser.get_do_annotations_priority())
    do_sentence_generator = SentenceGenerator(annotations=do_annotations, ontology=df.do_ontology,
                                              **do_sent_gen_common_prop)
    gene_desc.stats.total_number_do_annotations = len(do_annotations)
    gene_desc.stats.set_initial_do_ids = [term_id for terms in do_sentence_generator.terms_groups.values() for tvalues
                                          in terms.values() for term_id in tvalues]
    raw_disease_sent = do_sentence_generator.get_sentences(
        aspect='D', merge_groups_with_same_prefix=True, keep_only_best_group=False, **do_sent_common_props)
    disease_sent = "; ".join([sentence.text for sentence in raw_disease_sent])
    if disease_sent:
        joined_sent.append(disease_sent)
        gene_desc.do_description = disease_sent
    gene_desc.stats.set_final_do_ids = [term_id for sentence in raw_disease_sent for term_id in sentence.terms_ids]
    if "(multiple)" in disease_sent:
        gene_desc.stats.number_final_do_term_covering_multiple_initial_do_terms = \
            disease_sent.count("(multiple)")

    if not gene_desc.go_description:
        human_func_sent = None
        if len(orthologs_sp_fullname) == 1 and orthologs_sp_fullname[0] == "Homo sapiens":
            # human_orthologs = df.get_all_orthologs_for_gene(gene_id=gene.id, organism="Homo sapiens")
            best_orth = get_best_human_ortholog_for_info_poor(selected_orthologs, ensembl_hgnc_ids_map,
                                                              conf_parser.get_go_annotations_priority(), human_df_agr,
                                                              go_sent_gen_common_props)
            if best_orth:
                best_orth = "RGD:" + best_orth
                human_go_annotations = human_df_agr.get_annotations_for_gene(
                    gene_id=best_orth, annot_type=DataType.GO, priority_list=("EXP", "IDA", "IPI", "IMP", "IGI", "IEP",
                                                                              "HTP", "HDA", "HMP", "HGI", "HEP"))
                human_go_sent_generator = SentenceGenerator(annotations=human_go_annotations,
                                                            ontology=human_df_agr.go_ontology,
                                                            **go_sent_gen_common_props)
                raw_human_func_sent = human_go_sent_generator.get_sentences(aspect='F',
                                                                            merge_groups_with_same_prefix=True,
                                                                            keep_only_best_group=True,
                                                                            **go_sent_common_props)
                human_func_sent = " and ".join([sentence.text for sentence in raw_human_func_sent])
                if human_func_sent:
                    joined_sent.append("human " + human_df_agr.go_associations.subject_label_map[best_orth] + " " +
                                       human_func_sent)
        if not human_func_sent:
            protein_domains = df.protein_domains[gene.id[3:]]
            if protein_domains:
                dom_word = "domain"
                if len(protein_domains) > 1:
                    dom_word = "domains"
                joined_sent.append("is predicted to encode a protein with the following " + dom_word + ": " +
                                   concatenate_words_with_oxford_comma([ptdom[1] if ptdom[1] != "" else ptdom[0] for
                                                                        ptdom in protein_domains]))

    if len(joined_sent) > 0:
        desc = "; ".join(joined_sent) + "."
        if len(desc) > 0:
            gene_desc.description = gene.name + " " + desc
    else:
        gene_desc.description = None
    desc_writer.add_gene_desc(gene_desc)