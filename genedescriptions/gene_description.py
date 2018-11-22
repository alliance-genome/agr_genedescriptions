import inflect

from typing import List

from genedescriptions.commons import Module
from genedescriptions.descriptions_generator import OntologySentenceGenerator, ModuleSentences
from genedescriptions.sentence_generation_functions import concatenate_words_with_oxford_comma
from genedescriptions.stats import SingleDescStats


class GeneDescription(object):
    """gene description"""

    def __init__(self, gene_id: str, gene_name: str = "", add_gene_name: bool = False):
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.description = None
        self.go_description = None
        self.go_function_description = None
        self.go_process_description = None
        self.go_component_description = None
        self.do_description = None
        self.do_experimental_description = None
        self.do_biomarker_description = None
        self.do_orthology_description = None
        self.orthology_description = None
        self.tissue_expression_description = None
        self.gene_expression_cluster_description = None
        self.molecule_expression_cluster_description = None
        self.anatomy_expression_cluster_description = None
        self.protein_domain_description = None
        self.human_gene_function_description = None
        self.sister_species_description = None
        self.publications = None
        self.refs = None
        self.stats = SingleDescStats()
        self.paper_evidences = []
        self.accession_evidences = []
        self.add_gene_name = add_gene_name

    @staticmethod
    def _get_description(desc, desc_destination):
        if desc_destination:
            return desc_destination + "; " + desc
        else:
            return desc

    @staticmethod
    def _get_merged_ids(ids, ids_destination):
        if len(ids_destination) > 0:
            return list(set().union(set(ids_destination), set(ids)))
        else:
            return list(set(ids))

    def set_or_extend_module_description_and_final_stats(self, module: Module,
                                                         module_sentences: ModuleSentences = None,
                                                         description: str = None,
                                                         additional_postfix_terms_list: List[str] = None,
                                                         additional_postfix_final_word: str = None,
                                                         use_single_form: bool = False):
        """set description text and stats for a specific module

        if previous data is present in the specified module, the provided description and the stats are merged with the
        existing ones

        Args:
            module (Module): the description module to update
            module_sentences (ModuleSentences): optional - module sentences object from which to take the description
                and stats
            description (str): optional - description text to be added
            additional_postfix_terms_list (List[str]): optional - list of terms to be merged and added as postfix to the
                description
            additional_postfix_final_word: optional - word to be added at the end of the postfix (automatically
                converted to plural if the list of terms has more than one element)
            use_single_form (bool): whether to use a single form for the final word without transforming it to plural
        """
        desc = ""
        if module_sentences:
            desc = module_sentences.get_description()
        elif description:
            inflect_engine = inflect.engine()
            desc = description
            if additional_postfix_terms_list and len(additional_postfix_terms_list) > 0:
                desc += " " + concatenate_words_with_oxford_comma(additional_postfix_terms_list) + " " + \
                        (additional_postfix_final_word if use_single_form or len(additional_postfix_terms_list) == 1
                         else inflect_engine.plural_noun(additional_postfix_final_word))
        if desc:
            if self.description and self.description != self.gene_name:
                self.description = self.description[0:-1] + "; " + desc + "."
            else:
                self.description = self.gene_name + " " + desc + "." if self.add_gene_name else desc + "."
            if module == Module.GO_FUNCTION:
                self.go_function_description = self._get_description(desc, self.go_function_description)
                self.stats.set_final_go_ids_f = self._get_merged_ids(module_sentences.get_ids(experimental_only=False),
                                                                     self.stats.set_final_go_ids_f)
                self.stats.set_final_experimental_go_ids_f = self._get_merged_ids(module_sentences.get_ids(
                    experimental_only=True), self.stats.set_final_experimental_go_ids_f)
            elif module == Module.GO_PROCESS:
                self.go_process_description = self._get_description(desc, self.go_process_description)
                self.stats.set_final_go_ids_p = self._get_merged_ids(module_sentences.get_ids(experimental_only=False),
                                                                     self.stats.set_final_go_ids_p)
                self.stats.set_final_experimental_go_ids_p = self._get_merged_ids(module_sentences.get_ids(
                    experimental_only=True), self.stats.set_final_experimental_go_ids_p)
            elif module == Module.GO_COMPONENT:
                self.go_component_description = self._get_description(desc, self.go_component_description)
                self.stats.set_final_go_ids_c = self._get_merged_ids(module_sentences.get_ids(experimental_only=False),
                                                                     self.stats.set_final_go_ids_c)
                self.stats.set_final_experimental_go_ids_c = self._get_merged_ids(module_sentences.get_ids(
                    experimental_only=True), self.stats.set_final_experimental_go_ids_c)
            elif module == Module.EXPRESSION:
                self.tissue_expression_description = self._get_description(desc, self.tissue_expression_description)
                self.stats.set_final_expression_ids = self._get_merged_ids(
                    module_sentences.get_ids(experimental_only=False), self.stats.set_final_expression_ids)
                self.stats.set_final_experimental_go_ids_c = self._get_merged_ids(module_sentences.get_ids(
                    experimental_only=True), self.stats.set_final_experimental_go_ids_c)
            elif module == Module.EXPRESSION_CLUSTER_GENE:
                self.gene_expression_cluster_description = self._get_description(
                    desc, self.gene_expression_cluster_description)
            elif module == Module.EXPRESSION_CLUSTER_ANATOMY:
                self.anatomy_expression_cluster_description = self._get_description(
                    desc, self.anatomy_expression_cluster_description)
            elif module == Module.EXPRESSION_CLUSTER_MOLECULE:
                self.molecule_expression_cluster_description = self._get_description(
                    desc, self.molecule_expression_cluster_description)
            elif module == Module.DO_EXPERIMENTAL:
                self.do_experimental_description = self._get_description(desc, self.do_experimental_description)
                self.stats.set_final_do_ids = self._get_merged_ids(module_sentences.get_ids(experimental_only=False),
                                                                   self.stats.set_final_do_ids)
            elif module == Module.DO_BIOMARKER:
                self.do_biomarker_description = self._get_description(desc, self.do_biomarker_description)
                self.stats.set_final_do_ids = self._get_merged_ids(module_sentences.get_ids(experimental_only=False),
                                                                   self.stats.set_final_do_ids)
            elif module == Module.DO_ORTHOLOGY:
                self.do_orthology_description = self._get_description(desc, self.do_orthology_description)
                self.stats.set_final_do_ids = self._get_merged_ids(module_sentences.get_ids(experimental_only=False),
                                                                   self.stats.set_final_do_ids)
            elif module == Module.SISTER_SP:
                self.sister_species_description = self._get_description(desc, self.sister_species_description)
            elif module == Module.ORTHOLOGY:
                self.orthology_description = self._get_description(desc, self.orthology_description)
            elif module == Module.INFO_POOR_HUMAN_FUNCTION:
                self.human_gene_function_description = self._get_description(desc, self.human_gene_function_description)
            elif module == Module.PROTEIN_DOMAIN:
                self.protein_domain_description = self._get_description(desc, self.protein_domain_description)
            # Multimodule fields
            if module == Module.GO_PROCESS or module == Module.GO_FUNCTION or module == Module.GO_COMPONENT:
                self.go_description = "; ".join([sent for sent in [self.go_function_description,
                                                                   self.go_process_description,
                                                                   self.go_component_description] if sent])
            if module == Module.DO_EXPERIMENTAL or module == Module.DO_BIOMARKER or module == Module.DO_ORTHOLOGY:
                self.do_description = "; ".join([sent for sent in [self.do_experimental_description,
                                                                   self.do_biomarker_description,
                                                                   self.do_orthology_description] if sent])
                self.stats.number_final_do_term_covering_multiple_initial_do_terms = self.do_description.count(
                    "(multiple)")

    @staticmethod
    def _get_module_initial_stats(aspect: str, sentence_generator: OntologySentenceGenerator, main_qualifier: str = "",
                                  additional_qualifier: str = None):
        if not additional_qualifier:
            return len([elem for key, sets in sentence_generator.terms_groups[(aspect, main_qualifier)].items() for elem
                        in sets if (aspect, key, main_qualifier) in sentence_generator.prepostfix_sentences_map])
        else:
            return len(list(set().union(
                [elem for key, sets in sentence_generator.terms_groups[(aspect, main_qualifier)].items() for elem in
                 sets if (aspect, key, main_qualifier) in sentence_generator.prepostfix_sentences_map],
                [elem for key, sets in sentence_generator.terms_groups[
                    (aspect, additional_qualifier)].items() for elem in sets if (aspect, key, additional_qualifier) in
                 sentence_generator.prepostfix_sentences_map])))

    def set_initial_stats(self, module: Module, sentence_generator: OntologySentenceGenerator,
                          sentence_generator_exp_only: OntologySentenceGenerator = None):
        """set initial stats for a specific module

        Args:
            module: the module
            sentence_generator: the main sentence generator
            sentence_generator_exp_only: sentence generator with experimental evidence codes only

        Returns:

        """
        if module == Module.GO_FUNCTION:
            self.stats.num_initial_go_ids_f = self._get_module_initial_stats(
                aspect="F", additional_qualifier="contributes_to", sentence_generator=sentence_generator)
            self.stats.num_initial_experimental_go_ids_f = self._get_module_initial_stats(
                aspect="F", additional_qualifier="contributes_to", sentence_generator=sentence_generator_exp_only)
        elif module == Module.GO_COMPONENT:
            self.stats.num_initial_go_ids_c = self._get_module_initial_stats(
                aspect="C", additional_qualifier="colocalizes_with", sentence_generator=sentence_generator)
            self.stats.num_initial_experimental_go_ids_c = self._get_module_initial_stats(
                aspect="C", additional_qualifier="colocalizes_with", sentence_generator=sentence_generator_exp_only)
        elif module == Module.GO_PROCESS:
            self.stats.num_initial_go_ids_p = self._get_module_initial_stats(
                aspect="P", sentence_generator=sentence_generator)
            self.stats.num_initial_experimental_go_ids_p = self._get_module_initial_stats(
                aspect="P", sentence_generator=sentence_generator_exp_only)
        elif module == Module.EXPRESSION:
            self.stats.num_initial_expression_ids = self._get_module_initial_stats(
                aspect="A", main_qualifier="Verified", sentence_generator=sentence_generator)
        elif module == Module.DO_EXPERIMENTAL:
            self.stats.total_number_do_exp_bio_annotations += len(sentence_generator.annotations)
            self.stats.set_initial_do_ids = self._get_merged_ids(
                [term_id for terms in sentence_generator.terms_groups.values() for tvalues in terms.values() for
                 term_id in tvalues], self.stats.set_initial_do_ids)
        elif module == Module.DO_BIOMARKER:
            self.stats.total_number_do_exp_bio_annotations += len(sentence_generator.annotations)
            self.stats.set_initial_do_ids = self._get_merged_ids(
                [term_id for terms in sentence_generator.terms_groups.values() for tvalues in terms.values() for term_id
                 in tvalues], self.stats.set_initial_do_ids)
        elif module == Module.DO_ORTHOLOGY:
            self.stats.total_number_do_via_orth_annotations = len(sentence_generator.annotations)
            self.stats.set_initial_do_ids = self._get_merged_ids(
                [term_id for terms in sentence_generator.terms_groups.values() for tvalues in terms.values() for term_id
                 in tvalues], self.stats.set_initial_do_ids)
        self.stats.total_number_do_annotations = self.stats.total_number_do_exp_bio_annotations + \
                                                 self.stats.total_number_do_via_orth_annotations
        if module == Module.GO_PROCESS or module == Module.GO_FUNCTION or module == Module.GO_COMPONENT:
            self.stats.total_number_go_annotations = len(sentence_generator.annotations)
