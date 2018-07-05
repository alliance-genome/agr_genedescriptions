import json
from collections import namedtuple
from enum import Enum

import inflect
import re
import urllib.request

from namedlist import namedlist

from genedescriptions.config_parser import GenedescConfigParser
from genedescriptions.ontology_tools import *
from ontobio.ontol import Ontology


class DataType(Enum):
    GO = 1
    DO = 2


Gene = namedtuple('Gene', ['id', 'name', 'dead', 'pseudo'])


Sentence = namedlist('Sentence', ['prefix', 'terms_ids', 'postfix', 'text', 'aspect', 'evidence_group', 'terms_merged',
                                  'additional_prefix', 'qualifier', 'ancestors_covering_multiple_terms'])


class SingleDescStats(object):
    """statistics for a single gene description"""
    def __init__(self):
        self.total_number_go_annotations = 0
        self.number_initial_go_terms_f = 0
        self.number_initial_go_terms_p = 0
        self.number_initial_go_terms_c = 0
        self.number_final_go_terms_f = 0
        self.number_final_go_terms_p = 0
        self.number_final_go_terms_c = 0
        self.total_number_do_annotations = 0
        self.number_initial_do_terms = 0
        self.number_final_do_terms = 0
        self.number_final_do_term_covering_multiple_initial_do_terms_present = 0


class GeneDesc(object):
    """gene description"""
    def __init__(self, gene_id: str, gene_name: str = "", description: str = None, go_description: str = None,
                 go_function_description: str = None, go_process_description: str = None,
                 go_component_description: str = None, do_description: str = None, stats: SingleDescStats = None,
                 publications: str = "", refs: str = "", species: str = "", release_version: str = ""):
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.description = description
        self.go_description = go_description
        self.go_function_description = go_function_description
        self.go_process_description = go_process_description
        self.go_component_description = go_component_description
        self.do_description = do_description
        self.publications = publications
        self.refs = refs
        self.species = species
        self.release_version = release_version
        if stats:
            self.stats = stats
        else:
            self.stats = SingleDescStats()


class DescriptionsStats(object):
    """overall statistics for a set of gene descriptions"""
    def __init__(self):
        self.total_number__of_genes = 0
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
                            if re.match(re.escape(special_case[1]), ontology.node(annotation["object"]["id"])["label"]):
                                ev_group = evidence_codes_groups_map[annotation["evidence"]["type"]] + \
                                           str(special_case[0])
                                if ev_group not in self.evidence_groups_priority_list:
                                    self.evidence_groups_priority_list.insert(evidence_groups_priority_list.index(
                                        evidence_codes_groups_map[annotation["evidence"]["type"]]) + 1, ev_group)
                                break
                    self.terms_groups[(aspect, qualifier)][ev_group].add(annotation["object"]["id"])

    def get_sentences(self, aspect: str, qualifier: str = '', keep_only_best_group: bool = False,
                      merge_groups_with_same_prefix: bool = False, remove_parent_terms: bool = True,
                      merge_num_terms_threshold: int = 3, merge_min_distance_from_root: dict = None,
                      truncate_others_generic_word: str = "several",
                      truncate_others_aspect_words: Dict[str, str] = None,
                      remove_successive_overlapped_terms: bool = True,
                      exclude_terms_ids: List[str] = None,
                      add_multiple_if_covers_more_children: bool = False) -> List[Sentence]:
        """generate sentences for specific combination of aspect and qualifier

        Args:
            aspect (str): a data type aspect
            qualifier (str): qualifier
            keep_only_best_group (bool): whether to get only the evidence group with highest priority and discard
                the other evidence groups
            merge_groups_with_same_prefix (bool): whether to merge the phrases for evidence groups with the same prefix
            remove_parent_terms: whether to remove parent terms from the list of terms in each sentence if at least
                one children term is present
            merge_num_terms_threshold (int): whether to merge terms by common ancestor to reduce the number of terms in
                the set. The trimming algorithm will be applied only if the number of terms is greater than the
                specified number and the specified threshold is greater than 0
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
                    merged_terms = find_set_covering([(k, self.ontology.node(k)["label"], v) for k, v in
                                                      merged_terms_coverset.items()],
                                                     max_num_subsets=merge_num_terms_threshold)
                    for merged_term in merged_terms:
                        terms_already_covered.update(merged_terms_coverset[merged_term])
                    add_others = True
                if add_multiple_if_covers_more_children:
                    ancestors_covering_multiple_children = {self.ontology.label(ancestor) for ancestor in merged_terms
                                                            if len(merged_terms_coverset[ancestor]) > 1}
                logging.debug("Reduced number of terms by merging from " + str(len(terms)) + " to " +
                              str(len(merged_terms)))
                terms = merged_terms
            else:
                terms_already_covered.update(terms)
            if (aspect, evidence_group, qualifier) in self.prepostfix_sentences_map:
                sentences.append(
                    _get_single_sentence(node_ids=terms, ontology=self.ontology, aspect=aspect,
                                         evidence_group=evidence_group, qualifier=qualifier,
                                         prepostfix_sentences_map=self.prepostfix_sentences_map,
                                         terms_merged=True if 0 < merge_num_terms_threshold < len(terms) else False,
                                         add_others=add_others,
                                         truncate_others_generic_word=truncate_others_generic_word,
                                         truncate_others_aspect_words=truncate_others_aspect_words,
                                         ancestors_with_multiple_children=ancestors_covering_multiple_children))
                if keep_only_best_group:
                    return sentences
        if merge_groups_with_same_prefix:
            sentences = self.merge_sentences_with_same_prefix(sentences=sentences,
                                                              remove_parent_terms=remove_parent_terms)
        return sentences

    def merge_sentences_with_same_prefix(self, sentences: List[Sentence], remove_parent_terms: bool = True):
        """merge sentences with the same prefix

        Args:
            sentences (List[Sentence]): a list of sentences
            remove_parent_terms (bool): whether to remove parent terms if present in the merged set of terms
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
                         term_names=[self.ontology.label(node) for node in sent_merger.terms_ids],
                         postfix=SentenceGenerator.merge_postfix_phrases(sent_merger.postfix_list),
                         additional_prefix=sent_merger.additional_prefix,
                         ancestors_with_multiple_children=sent_merger.ancestors_covering_multiple_terms),
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
                     ancestors_with_multiple_children: Set[str] = None) -> str:
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
    if len(term_names) > 2:
        return prefix + ", ".join(term_names[0:-1]) + ", and " + term_names[len(term_names) - 1] + postfix
    elif len(term_names) > 1:
        return prefix + " and ".join(term_names) + postfix
    else:
        return prefix + term_names[0] + postfix


def _get_single_sentence(node_ids: List[str], ontology: Ontology, aspect: str, evidence_group: str, qualifier: str,
                         prepostfix_sentences_map: Dict[Tuple[str, str, str], Tuple[str, str]],
                         terms_merged: bool = False, add_others: bool = False,
                         truncate_others_generic_word: str = "several",
                         truncate_others_aspect_words: Dict[str, str] = None,
                         ancestors_with_multiple_children: Set[str] = None) -> Union[Sentence, None]:
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
        term_labels = [ontology.label(node_id) for node_id in node_ids]
        return Sentence(prefix=prefix, terms_ids=node_ids, postfix=postfix,
                        text=compose_sentence(prefix=prefix, term_names=term_labels, postfix=postfix,
                                              additional_prefix=additional_prefix,
                                              ancestors_with_multiple_children=ancestors_with_multiple_children),
                        aspect=aspect, evidence_group=evidence_group, terms_merged=terms_merged,
                        additional_prefix=additional_prefix, qualifier=qualifier,
                        ancestors_covering_multiple_terms=ancestors_with_multiple_children)
    else:
        return None


def generate_ortholog_sentence(orthologs: List[List[str]], orthologs_sp_fullname: str, human_genes_props):
    orth_sentence = None
    if orthologs_sp_fullname == "Homo sapiens":
        if len(orthologs) > 3:
            gene_families = defaultdict(list)
            for ortholog in orthologs:
                if human_genes_props[ortholog[0]]:
                    gene_families[human_genes_props[ortholog[0]][2]].append(human_genes_props[ortholog[0]])
            if len(gene_families.values()) > 0:
                gene_family_names = list(gene_families.keys())
                if len(gene_family_names) > 3:
                    gene_family_names = gene_family_names[0:3]
                gene_names = [ortholog[0] + " (" + ortholog[1] + ")" for orthologs in gene_families.values() for
                              ortholog in orthologs]
                if len(gene_names) > 3:
                    gene_names = gene_names[0:3]
                family_word = "family"
                if len(gene_family_names) > 1:
                    family_word = "families"
                if len(gene_family_names) > 2:
                    ortholog_families_str = ", ".join(gene_family_names[0:-1]) + ", and " + gene_family_names[-1]
                else:
                    ortholog_families_str = " and ".join(gene_family_names)
                if len(gene_names) > 2:
                    ortholog_genes_str = ", ".join(gene_names[0:-1]) + ", and " + gene_names[-1]
                else:
                    ortholog_genes_str = " and ".join(gene_names)
                orth_sentence = "is an ortholog of members of the human " + ortholog_families_str + " gene " + \
                                family_word + " including " + ortholog_genes_str
        else:
            symbol_name_arr = sorted([human_genes_props[best_orth[0]][0] + " (" + human_genes_props[best_orth[0]][1] +
                                      ")" for best_orth in orthologs if human_genes_props[best_orth[0]]])
            if len(symbol_name_arr) > 0:
                if len(symbol_name_arr) > 2:
                    orth_sentence = "is an ortholog of human " + ", ".join(symbol_name_arr[0:-1]) + ", and " + \
                                    symbol_name_arr[-1]
                else:
                    orth_sentence = "is an ortholog of human " + " and ".join(symbol_name_arr)
    else:
        fullname_arr = orthologs_sp_fullname.split(" ")
        if len(fullname_arr[0]) > 2:
            fullname_arr[0] = fullname_arr[0][0] + "."
            orthologs_sp_fullname = " ".join(fullname_arr)
        if len(orthologs) > 3:
            gene_classes = defaultdict(list)
            for ortholog in orthologs:
                gene_class_data = json.loads(urllib.request.urlopen("http://rest.wormbase.org/rest/field/gene/" +
                                                                    ortholog[0] + "/gene_class").read())
                if "gene_class" in gene_class_data and gene_class_data["gene_class"]["data"] and "tag" in \
                        gene_class_data["gene_class"]["data"] and "label" in \
                        gene_class_data["gene_class"]["data"]["tag"]:
                    gene_classes[gene_class_data["gene_class"]["data"]["tag"]["label"]].append(ortholog)
            classes_gene_symbols = list(gene_classes.keys())
            if len(classes_gene_symbols) > 0:
                classes_word = "class"
                if len(classes_gene_symbols) > 1:
                    classes_word = "classes"
                if len(classes_gene_symbols) > 2:
                    orth_sentence = "is an ortholog of " + orthologs_sp_fullname + " " + \
                                    ", ".join(classes_gene_symbols[0:-1]) + ", and " + classes_gene_symbols[-1] + \
                                    " gene " + classes_word
                else:
                    orth_sentence = "is an ortholog of " + orthologs_sp_fullname + " " + \
                                    " and ".join(classes_gene_symbols) + " gene " + classes_word
                return orth_sentence
        orthologs_symbols = [orth[1] for orth in orthologs]
        if len(orthologs_symbols) > 2:
            if len(orthologs_symbols) > 3:
                orthologs_symbols = orthologs_symbols[0:3]
            orth_sentence = "is an ortholog of " + orthologs_sp_fullname + " " + ", ".join(orthologs_symbols[0:-1]) + \
                            ", and " + orthologs_symbols[-1]
        else:
            orth_sentence = "is an ortholog of " + orthologs_sp_fullname + " " + " and ".join(orthologs_symbols)
    return orth_sentence


def compose_wormbase_description(gene: Gene, conf_parser: GenedescConfigParser, species, organism, df,
                                 orthologs_sp_fullname, go_sent_gen_common_props, go_sent_common_props,
                                 human_genes_props, do_sent_gen_common_prop, do_sent_common_props, sister_sp_fullname,
                                 sister_df, desc_writer):
    gene_desc = GeneDesc(gene_id=gene.id, gene_name=gene.name,
                         publications=", ".join([annot["publication"] for annot in df.get_annotations_for_gene(
                             gene.id, annot_type=DataType.GO,
                             priority_list=conf_parser.get_go_evidence_groups_priority_list())]),
                         refs=", ".join([annot["refs"] for annot in df.get_annotations_for_gene(
                             gene.id, annot_type=DataType.GO,
                             priority_list=conf_parser.get_go_evidence_groups_priority_list())]),
                         species=species[organism]["full_name"],
                         release_version=conf_parser.get_release("wb_data_fetcher"))
    joined_sent = []

    best_orthologs, selected_orth_name = df.get_best_orthologs_for_gene(
        gene.id, orth_species_full_name=orthologs_sp_fullname)
    if best_orthologs:
        orth_sent = generate_ortholog_sentence(best_orthologs, selected_orth_name, human_genes_props)
        if orth_sent:
            joined_sent.append(orth_sent)
    go_annotations = df.get_annotations_for_gene(gene_id=gene.id, annot_type=DataType.GO,
                                                 priority_list=conf_parser.get_go_annotations_priority())
    go_sent_generator = SentenceGenerator(annotations=go_annotations, ontology=df.go_ontology,
                                          **go_sent_gen_common_props)
    gene_desc.stats.total_number_go_annotations = len(go_annotations)
    gene_desc.stats.number_initial_go_terms = {aspect: len(terms) for aspect, terms in
                                               go_sent_generator.terms_groups.items()}
    raw_func_sent = go_sent_generator.get_sentences(aspect='F', merge_groups_with_same_prefix=True,
                                                    keep_only_best_group=True, **go_sent_common_props)
    gene_desc.stats.number_final_go_terms_f += sum([len(sentence.terms_ids) for sentence in raw_func_sent])
    func_sent = " and ".join([sentence.text for sentence in raw_func_sent])
    if func_sent:
        joined_sent.append(func_sent)
    contributes_to_raw_func_sent = go_sent_generator.get_sentences(
        aspect='F', qualifier='contributes_to', merge_groups_with_same_prefix=True, keep_only_best_group=True,
        **go_sent_common_props)
    gene_desc.stats.number_final_go_terms_f += sum([len(sentence.terms_ids) for sentence in
                                                       contributes_to_raw_func_sent])
    contributes_to_func_sent = " and ".join([sentence.text for sentence in contributes_to_raw_func_sent])
    if contributes_to_func_sent:
        joined_sent.append(contributes_to_func_sent)
    raw_proc_sent = go_sent_generator.get_sentences(aspect='P', merge_groups_with_same_prefix=True,
                                                    keep_only_best_group=True, **go_sent_common_props)
    gene_desc.stats.number_final_go_terms_p += sum([len(sentence.terms_ids) for sentence in raw_proc_sent])
    proc_sent = " and ".join([sentence.text for sentence in raw_proc_sent])
    if proc_sent:
        joined_sent.append(proc_sent)
    raw_comp_sent = go_sent_generator.get_sentences(
        aspect='C', merge_groups_with_same_prefix=True, keep_only_best_group=True, **go_sent_common_props)
    gene_desc.stats.number_final_go_terms_c += sum([len(sentence.terms_ids) for sentence in raw_comp_sent])
    comp_sent = " and ".join([sentence.text for sentence in raw_comp_sent])
    if comp_sent:
        joined_sent.append(comp_sent)
    colocalizes_with_raw_comp_sent = go_sent_generator.get_sentences(
        aspect='C', qualifier='colocalizes_with', merge_groups_with_same_prefix=True,
        keep_only_best_group=True, **go_sent_common_props)
    gene_desc.stats.number_final_go_terms_c += sum([len(sentence.terms_ids) for sentence in
                                                       colocalizes_with_raw_comp_sent])
    colocalizes_with_comp_sent = " and ".join([sentence.text for sentence in colocalizes_with_raw_comp_sent])
    if colocalizes_with_comp_sent:
        joined_sent.append(colocalizes_with_comp_sent)
    do_annotations = df.get_annotations_for_gene(gene_id=gene.id, annot_type=DataType.DO,
                                                 priority_list=conf_parser.get_do_annotations_priority())
    do_sentence_generator = SentenceGenerator(annotations=do_annotations, ontology=df.do_ontology,
                                              **do_sent_gen_common_prop)
    gene_desc.stats.total_number_do_annotations = len(do_annotations)
    gene_desc.stats.number_initial_do_terms = sum([len(terms) for terms in
                                                   do_sentence_generator.terms_groups.values()])
    raw_disease_sent = do_sentence_generator.get_sentences(
        aspect='D', merge_groups_with_same_prefix=True, keep_only_best_group=False, **do_sent_common_props)
    disease_sent = "; ".join([sentence.text for sentence in raw_disease_sent])
    if disease_sent:
        joined_sent.append(disease_sent)
    gene_desc.stats.number_final_do_terms += sum([len(sentence.terms_ids) for sentence in raw_disease_sent])
    if "(multiple)" in disease_sent:
        gene_desc.stats.number_final_do_term_covering_multiple_initial_do_terms_present = \
            disease_sent.count("(multiple)")
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
    if len(joined_sent) > 0:
        desc = "; ".join(joined_sent) + "."
        if len(desc) > 0:
            gene_desc.description = desc[0].upper() + desc[1:]
    else:
        gene_desc.description = None
    desc_writer.add_gene_desc(gene_desc)