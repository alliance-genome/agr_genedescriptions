import copy
import inflect
import re
from namedlist import namedlist
from genedescriptions.ontology_tools import *
from ontobio.ontol import Ontology
from ontobio.assocmodel import AssociationSet

Sentence = namedlist('Sentence', ['prefix', 'terms_ids', 'postfix', 'text', 'aspect', 'evidence_group', 'terms_merged',
                                  'additional_prefix', 'qualifier'])


class SingleDescStats(object):
    """statistics for a single gene description"""
    def __init__(self):
        self.num_terms_notrim_nogroup_priority_nomerge = defaultdict(int)
        self.num_terms_trim_nogroup_priority_nomerge = defaultdict(int)
        self.num_terms_trim_group_priority_merge = defaultdict(int)
        self.total_num_go_annotations = 0
        self.num_prioritized_go_annotations = 0
        self.terms_notrim_nogroup_priority_nomerge = defaultdict(list)
        self.terms_trim_nogroup_priority_nomerge = defaultdict(list)
        self.terms_trim_group_priority_merge = defaultdict(list)


class GeneDesc(object):
    """gene description"""
    def __init__(self, gene_id: str, gene_name: str = "", description: str = "", stats: SingleDescStats = None):
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.description = description
        if stats:
            self.stats = stats
        else:
            self.stats = SingleDescStats()


class DescriptionsStats(object):
    """overall statistics for a set of gene descriptions"""
    def __init__(self):
        self.num_genes_with_go_sentence = 0
        self.average_num_go_terms_if_desc_trim_group_priority_merge = 0
        self.average_num_go_terms_if_desc_trim_nogroup_priority_nomerge = 0
        self.average_num_go_terms_if_desc_notrim_nogroup_priority_nomerge = 0


class SentenceMerger(object):
    def __init__(self):
        self.postfix_list = []
        self.terms_ids: set = set()
        self.term_postfix_dict = {}
        self.evidence_groups = []
        self.term_evgroup_dict = {}
        self.additional_prefix = ""


class SentencesCollection(object):
    """a group of sentences indexed by aspect"""
    def __init__(self, evidence_groups_list, prepostfix_sentences_map, remove_parent_terms: bool, ontology: Ontology):
        self.evidence_groups_list = evidence_groups_list
        self.prepostfix_sentences_map = prepostfix_sentences_map
        self.sentences_map = {}
        self.remove_parent_terms = remove_parent_terms
        self.ontology = ontology

    def set_sentence(self, sentence: Sentence) -> None:
        """add a sentence to the collection

        :param sentence: the sentence to add
        :type sentence: Sentence
        """
        if sentence is not None:
            self.sentences_map[(sentence.aspect, sentence.evidence_group, sentence.qualifier)] = sentence

    def get_sentences(self, aspect: str, qualifier: str = '', keep_only_best_group: bool = False,
                      merge_groups_with_same_prefix: bool = False,
                      desc_stats: SingleDescStats = None) -> List[Sentence]:
        """get all sentences containing the specified aspect

        :param aspect: a GO aspect
        :type aspect: str
        :param qualifier: qualifier
        :type qualifier: str
        :param keep_only_best_group: whether to get only the evidence group with highest priority and discard
            the other evidence groups
        :type keep_only_best_group: bool
        :param merge_groups_with_same_prefix: whether to merge the phrases for evidence groups with the same prefix
        :type merge_groups_with_same_prefix: bool
        :param desc_stats: an object containing the description statistics where to save the total number of annotations
            for the gene
        :type desc_stats: SingleDescStats
        :return: the list of sentences containing the specified GO aspect
        :rtype: List[GOSentence]
        """
        sentences = []
        merged_sentences = defaultdict(SentenceMerger)
        if desc_stats:
            desc_stats.num_terms_trim_nogroup_priority_nomerge[aspect] = 0
            desc_stats.terms_trim_nogroup_priority_nomerge[aspect] = []
            for ga, eg, qu in self.sentences_map.keys():
                if ga == aspect and qu == qualifier:
                    desc_stats.num_terms_trim_nogroup_priority_nomerge[aspect] += \
                        len(self.sentences_map[(aspect, eg, qu)].terms_ids)
                    desc_stats.terms_trim_nogroup_priority_nomerge[aspect].extend(
                        self.sentences_map[(aspect, eg, qu)].terms_ids)
        terms_in_previous_ev_groups = set()
        for eg in self.evidence_groups_list:
            if (aspect, eg, qualifier) in self.sentences_map:
                if merge_groups_with_same_prefix:
                    prefix = self.prepostfix_sentences_map[(aspect, eg, qualifier)][0]
                    merged_sentences[prefix].postfix_list.append(self.prepostfix_sentences_map[(aspect, eg,
                                                                                                qualifier)][1])
                    merged_sentences[prefix].terms_ids.update([term for term in self.sentences_map[
                        (aspect, eg, qualifier)].terms_ids if term not in terms_in_previous_ev_groups])
                    for term in self.sentences_map[(aspect, eg, qualifier)].terms_ids:
                        merged_sentences[prefix].term_postfix_dict[term] = self.prepostfix_sentences_map[
                            (aspect, eg, qualifier)][1]
                    merged_sentences[prefix].evidence_groups.append(eg)
                    for term in self.sentences_map[(aspect, eg, qualifier)].terms_ids:
                        merged_sentences[prefix].term_evgroup_dict[term] = eg
                    terms_in_previous_ev_groups.update(merged_sentences[prefix].terms_ids)
                    if self.sentences_map[(aspect, eg, qualifier)].additional_prefix:
                        merged_sentences[prefix].additional_prefix = \
                            self.sentences_map[(aspect, eg, qualifier)].additional_prefix
                else:
                    sentence_copy = copy.copy(self.sentences_map[(aspect, eg, qualifier)])
                    sentence_copy.terms_ids = [term for term in sentence_copy.terms_ids if term not in
                                               terms_in_previous_ev_groups]
                    sentences.append(sentence_copy)
                    terms_in_previous_ev_groups.update(self.sentences_map[(aspect, eg, qualifier)].terms_ids)
                if keep_only_best_group:
                    break
        if merge_groups_with_same_prefix:
            for prefix, sent_merger in merged_sentences.items():
                # rem parents
                if self.remove_parent_terms:
                    terms_no_ancestors = sent_merger.terms_ids - set([ancestor for node_id in sent_merger.terms_ids for
                                                                      ancestor in self.ontology.ancestors(node_id)])
                    if len(sent_merger.terms_ids) > len(terms_no_ancestors):
                        logging.debug("Removed " + str(len(sent_merger.terms_ids) - len(terms_no_ancestors)) +
                                      " parents from terms while merging sentences with same prefix")
                        sent_merger.terms = terms_no_ancestors
            sentences = [Sentence(prefix=prefix, terms_ids=list(sent_merger.terms_ids),
                                  postfix=SentencesCollection.merge_postfix_phrases(sent_merger.postfix_list),
                                  text=compose_sentence(prefix=prefix,
                                                        term_names=[self.ontology.label(node) for node in
                                                                    sent_merger.terms_ids],
                                                        postfix=SentencesCollection.merge_postfix_phrases(
                                                                 sent_merger.postfix_list),
                                                        additional_prefix=sent_merger.additional_prefix),
                                  aspect=aspect, evidence_group=", ".join(sent_merger.evidence_groups),
                                  terms_merged=True, additional_prefix=sent_merger.additional_prefix,
                                  qualifier=qualifier)
                         for prefix, sent_merger in merged_sentences.items() if len(sent_merger.terms_ids) > 0]
        if desc_stats:
            desc_stats.num_terms_trim_group_priority_merge[aspect] = 0
            desc_stats.terms_trim_group_priority_merge[aspect] = []
            for sentence in sentences:
                desc_stats.num_terms_trim_group_priority_merge[sentence.aspect] += len(sentence.terms_ids)
                desc_stats.terms_trim_group_priority_merge[sentence.aspect].extend(sentence.terms_ids)
        return sentences

    @staticmethod
    def merge_postfix_phrases(postfix_phrases: List[str]) -> str:
        """merge postfix phrases and remove possible redundant text at the beginning at at the end of the phrases

        :param postfix_phrases: the phrases to merge
        :type postfix_phrases: List[str]
        :return: the merged postfix phrase
        :rtype: str
        """
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


def generate_sentences(annotations: List[Dict], ontology: Ontology, evidence_groups_priority_list: List[str],
                       prepostfix_sentences_map: Dict[Tuple[str, str, str], Tuple[str, str]],
                       evidence_codes_groups_map: Dict[str, str],
                       prepostfix_special_cases_sent_map: Dict[Tuple[str, str, str], Tuple[int, str, str, str]] = None,
                       remove_parent_terms: bool = True, merge_num_terms_threshold: int = 3,
                       merge_min_distance_from_root: dict = None, desc_stats: SingleDescStats = None,
                       truncate_others_generic_word: str = "several",
                       truncate_others_aspect_words: Dict[str, str] = None) -> SentencesCollection:
    """generate sentences from a list of annotations

    :param annotations: the list of annotations for a given gene
    :type annotations: List[Dict]
    :param ontology: the ontology linked to the annotations
    :type ontology: Ontology
    :param evidence_groups_priority_list: the list of evidence groups to consider, sorted by priority. Sentences of the
        first group (with highest priority) will be returned in first position and so on
    :type evidence_groups_priority_list: List[str]
    :param prepostfix_sentences_map: a map with prefix and postfix phrases, where keys are tuples of
        go_aspect, evidence_group and values are tuples prefix, postfix
    :type prepostfix_sentences_map: Dict[Tuple[str, str, str], Tuple[str, str]]
    :param prepostfix_special_cases_sent_map: a map for special prefix and postfix cases, where keys are tuples of
        go_aspect, evidence_group and values are tuples of id, match_regex, prefix, postfix. Match_regex is a regular
        expression that defines the match for the special case
    :type prepostfix_special_cases_sent_map: Dict[Tuple[str, str, str], Tuple[int, str, str, str]]
    :param evidence_codes_groups_map: a map between evidence codes and the groups they belong to
    :type evidence_codes_groups_map: Dict[str, str]
    :param remove_parent_terms: whether to remove parent terms from the list of terms in each sentence if at least
        one children term is present
    :type remove_parent_terms: bool
    :param merge_num_terms_threshold: whether to merge terms by common ancestor to
        reduce the number of terms in the set. The trimming algorithm will be applied only if the number of terms is
        greater than the specified number and the specified threshold is greater than 0
    :type merge_num_terms_threshold: int
    :param merge_min_distance_from_root: minimum distance from root terms for the selection of common ancestors
        during merging operations. Three values must be provided in the form of a dictionary with keys 'F', 'P', and
        'C' for go aspect names and values integers indicating the threshold for each aspect
    :type merge_min_distance_from_root: dict
    :param desc_stats: an object containing the description statistics where to save the total number of annotations
        for the gene
    :type desc_stats: SingleDescStats
    :param truncate_others_generic_word: a generic word to indicate that the set of terms reported in the sentence is
        only a subset of the original terms, e.g., 'several'
    :type truncate_others_generic_word: str
    :param truncate_others_aspect_words: one word for each aspect describing the kind of terms that are included in the
        aspect
    :type truncate_others_aspect_words: Dict[str, str]
    :return: a collection of sentences
    :rtype: SentencesCollection
    """
    if len(annotations) > 0:
        if not merge_min_distance_from_root:
            merge_min_distance_from_root = {'F': 1, 'P': 1, 'C': 2, 'D': 3}
        if not truncate_others_aspect_words:
            truncate_others_aspect_words = {'F': 'functions', 'P': 'processes', 'C': 'components'}
        terms_groups = defaultdict(set)
        for annotation in annotations:
            if annotation["evidence"]["type"] in evidence_codes_groups_map:
                map_key = (annotation["aspect"], evidence_codes_groups_map[annotation["evidence"]["type"]],
                           "_".join(sorted(annotation["qualifiers"])).lower())
                if prepostfix_special_cases_sent_map and map_key in prepostfix_special_cases_sent_map:
                    for special_case in prepostfix_special_cases_sent_map[map_key]:
                        if re.match(re.escape(special_case[1]), ontology.node(annotation["object"]["id"])["label"]):
                            map_key = (annotation["aspect"], evidence_codes_groups_map[
                                annotation["evidence"]["type"]] + str(special_case[0]),
                                       "_".join(sorted(annotation["qualifier"] if "qualifier" in annotation
                                                       else "")).lower())
                            if evidence_codes_groups_map[annotation["evidence"]["type"]] + str(special_case[0]) \
                                    not in evidence_groups_priority_list:
                                evidence_groups_priority_list.insert(evidence_groups_priority_list.index(
                                    evidence_codes_groups_map[annotation["evidence"]["type"]]) + 1,
                                                                     evidence_codes_groups_map[
                                                                         annotation["evidence"]["type"]] +
                                                                     str(special_case[0]))
                            break
                terms_groups[map_key].add(annotation["object"]["id"])
        sentences = SentencesCollection(evidence_groups_priority_list, prepostfix_sentences_map,
                                        remove_parent_terms=remove_parent_terms, ontology=ontology)
        for ((aspect, evidence_group, qualifier), terms) in terms_groups.items():
            if desc_stats:
                desc_stats.num_terms_notrim_nogroup_priority_nomerge[aspect] += len(terms)
                desc_stats.terms_notrim_nogroup_priority_nomerge[aspect].extend(terms)
            add_others = False
            if remove_parent_terms:
                terms_no_ancestors = terms - set([ancestor for node_id in terms for ancestor in
                                                  ontology.ancestors(node_id)])
                if len(terms) > len(terms_no_ancestors):
                    logging.debug("Removed " + str(len(terms) - len(terms_no_ancestors)) + " parents from terms")
                    terms = terms_no_ancestors
            if merge_num_terms_threshold > 0:
                merged_terms_coverset = get_merged_nodes_by_common_ancestor(
                    node_ids=terms, ontology=ontology, min_distance_from_root=merge_min_distance_from_root[aspect],
                    min_number_of_terms=merge_num_terms_threshold)
                if len(merged_terms_coverset.keys()) <= merge_num_terms_threshold:
                    merged_terms = list(merged_terms_coverset.keys())
                else:
                    merged_terms = find_set_covering([(k, ontology.node(k)["label"], v) for k, v in
                                                      merged_terms_coverset.items()],
                                                     max_num_subsets=merge_num_terms_threshold)
                    add_others = True
                if 0 < len(merged_terms) < len(terms):
                    logging.debug("Reduced number of terms by merging from " + str(len(terms)) + " to " +
                                  str(len(merged_terms)))
                    terms = merged_terms
            sentences.set_sentence(_get_single_sentence(node_ids=terms, ontology=ontology, aspect=aspect,
                                                        evidence_group=evidence_group, qualifier=qualifier,
                                                        prepostfix_sentences_map=prepostfix_sentences_map,
                                                        terms_merged=True if 0 < merge_num_terms_threshold < len(
                                                               terms) else False, add_others=add_others,
                                                        truncate_others_generic_word=truncate_others_generic_word,
                                                        truncate_others_aspect_words=truncate_others_aspect_words))
        return sentences


def compose_sentence(prefix: str, additional_prefix: str, term_names: List[str], postfix: str) -> str:
    """compose the text of a sentence given its prefix, terms, and postfix

    :param prefix: the prefix of the sentence
    :type prefix: str
    :param term_names: a list of go terms
    :type term_names: List[str]
    :param postfix: the postfix of the sentence
    :type postfix: str
    :param additional_prefix: an additional prefix to be used for special cases
    :type additional_prefix: str
    :return: the text of the go sentence
    :rtype: str"""
    prefix = prefix + additional_prefix + " "
    term_names = sorted(term_names)
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
                         truncate_others_aspect_words: Dict[str, str] = None) -> Union[Sentence, None]:
    """build a sentence

    :param node_ids: list of ids for the terms to be combined in the sentence
    :type node_ids: List[str]
    :param ontology: the ontology containing the nodes
    :type ontology: Ontology
    :param aspect: aspect
    :type aspect: str
    :param evidence_group: evidence group
    :type evidence_group: str
    :param qualifier: qualifier
    :type qualifier: str
    :param prepostfix_sentences_map: map for prefix and postfix phrases
    :type prepostfix_sentences_map: Dict[Tuple[str, str, str], Tuple[str, str]]
    :param terms_merged: whether the terms set has been merged to reduce its size
    :type terms_merged: bool
    :param add_others: whether to say that there are other terms which have been omitted from the sentence
    :type add_others: bool
    :param truncate_others_generic_word: a generic word to indicate that the set of terms reported in the sentence is
        only a subset of the original terms, e.g., 'several'
    :type truncate_others_generic_word: str
    :param truncate_others_aspect_words: one word for each aspect describing the kind of terms that are included in the
        aspect
    :type truncate_others_aspect_words: Dict[str, str]
    :return: the combined go sentence
    :rtype: Union[GOSentence, None]
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
                                              additional_prefix=additional_prefix),
                        aspect=aspect, evidence_group=evidence_group, terms_merged=terms_merged,
                        additional_prefix=additional_prefix, qualifier=qualifier)
    else:
        return None
