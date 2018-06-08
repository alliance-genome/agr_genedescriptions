import inflect
import re
from namedlist import namedlist
from genedescriptions.ontology_tools import *
from ontobio.ontol import Ontology

Sentence = namedlist('Sentence', ['prefix', 'terms_ids', 'postfix', 'text', 'aspect', 'evidence_group', 'terms_merged',
                                  'additional_prefix', 'qualifier', 'ancestors_covering_multiple_terms'])


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
    def __init__(self, gene_id: str, gene_name: str = "", description: str = "", go_description: str = "",
                 disease_description: str = "", stats: SingleDescStats = None):
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.description = description
        self.go_description = go_description
        self.disease_description = disease_description
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
                      desc_stats: SingleDescStats = None,
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
            desc_stats (SingleDescStat): an object containing the description statistics where to save the total number
                of annotations for the gene
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
            truncate_others_aspect_words = {'F': 'functions', 'P': 'processes', 'C': 'components'}
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
            if desc_stats:
                desc_stats.num_terms_notrim_nogroup_priority_nomerge[aspect] += len(terms)
                desc_stats.terms_notrim_nogroup_priority_nomerge[aspect].extend(terms)
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
