import inflect
import re

from genedescriptions.commons import Sentence
from genedescriptions.ontology_tools import *
from ontobio.ontol import Ontology
from genedescriptions.sentence_generation_functions import _get_single_sentence, compose_sentence


logger = logging.getLogger("Description Generator")


class ModuleSentences(object):
    def __init__(self, sentences):
        self.sentences = sentences

    def get_description(self):
        return " and ".join([sentence.text for sentence in self.sentences])

    def get_ids(self, experimental_only: bool = False):
        return [term_id for sentence in self.sentences for term_id in sentence.terms_ids if not experimental_only or
                sentence.evidence_group.startswith("EXPERIMENTAL")]

    def contains_sentences(self):
        return len(self.sentences) > 0


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

    def get_module_sentences(self, aspect: str, qualifier: str = '', keep_only_best_group: bool = False,
                             merge_groups_with_same_prefix: bool = False, remove_parent_terms: bool = True,
                             remove_child_terms: bool = False, merge_num_terms_threshold: int = 3,
                             merge_max_num_terms: int = 3, merge_min_distance_from_root: dict = None,
                             truncate_others_generic_word: str = "several",
                             truncate_others_aspect_words: Dict[str, str] = None,
                             remove_successive_overlapped_terms: bool = True, exclude_terms_ids: List[str] = None,
                             add_multiple_if_covers_more_children: bool = False, rename_cell: bool = False,
                             blacklisted_ancestors: List[str] = None) -> ModuleSentences:
        """generate description for a specific combination of aspect and qualifier

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
            blacklisted_ancestors (List[str]): list of common ancestor terms to be blacklisted
        Returns:
            ModuleSentences: the module sentences
        """
        if not merge_min_distance_from_root:
            merge_min_distance_from_root = {'F': 1, 'P': 1, 'C': 2, 'D': 3, 'A': 3}
        if not truncate_others_aspect_words:
            truncate_others_aspect_words = {'F': 'functions', 'P': 'processes', 'C': 'components', 'D': 'diseases',
                                            'A': 'tissues'}
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
                    logger.debug("Removed " + str(len(terms) - len(terms_no_ancestors)) + " parents from terms")
                    terms = terms_no_ancestors
            if 0 < merge_num_terms_threshold <= len(terms):
                merged_terms_coverset = get_merged_nodes_by_common_ancestor(
                    node_ids=list(terms), ontology=self.ontology,
                    min_distance_from_root=merge_min_distance_from_root[aspect],
                    min_number_of_terms=merge_num_terms_threshold, nodeids_blacklist=blacklisted_ancestors)
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
                logger.debug("Reduced number of terms by merging from " + str(len(terms)) + " to " +
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
                    return ModuleSentences(sentences)
        if merge_groups_with_same_prefix:
            sentences = self.merge_sentences_with_same_prefix(sentences=sentences,
                                                              remove_parent_terms=remove_parent_terms,
                                                              rename_cell=rename_cell)
        return ModuleSentences(sentences)

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
                    logger.debug("Removed " + str(len(sent_merger.terms_ids) - len(terms_no_ancestors)) +
                                 " parents from terms while merging sentences with same prefix")
                    sent_merger.terms_ids = terms_no_ancestors
        return [Sentence(prefix=prefix, terms_ids=list(sent_merger.terms_ids),
                         postfix=SentenceGenerator.merge_postfix_phrases(sent_merger.postfix_list),
                         text=compose_sentence(prefix=prefix,
                                               term_names=[self.ontology.label(node, id_if_null=True) for node in
                                                           sent_merger.terms_ids],
                                               postfix=SentenceGenerator.merge_postfix_phrases(
                                                   sent_merger.postfix_list),
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


