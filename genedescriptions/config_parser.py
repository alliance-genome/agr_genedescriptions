import yaml
from typing import Dict, Union, Tuple, List


class GenedescConfigParser(object):
    def __init__(self, file_path):
        with open(file_path) as conf_file:
            self.config = yaml.load(conf_file)

    def get_textpresso_api_token(self) -> str:
        """get textpresso API token

        Returns:
            str: the textpresso API token
        """
        return self.config["generic_data_fetcher"]["textpresso_api_token"]

    def get_data_fetcher(self) -> str:
        """get the data fetcher type from the configuration file

        Returns:
            str: the data fetcher name
        """
        return self.config["generic_data_fetcher"]["data_fetcher"]

    def get_cache_location(self) -> str:
        """get the root location of the cache files

        Returns:
            str: the path to the cache location
        """
        return self.config["generic_data_fetcher"]["cache_location"]

    def get_wb_species(self) -> Dict[str, Dict[str, str]]:
        """get species list for WormBase data fetcher

        Returns:
            Dict[str, Dict[str, str]]: the species dictionary, with species name as key and a dictionary with
                project_id as value
        """
        return self.config["wb_data_fetcher"]["organisms"]

    def get_go_evidence_codes(self) -> Dict[str, Dict[str, Union[str, id]]]:
        """get the configured evidence codes

        Returns:
            Dict[str, Dict[str, Union[str, id]]]: a dictionary of the configured evidence codes, with evidence code
                names as keys and a dictionary with 'evidence group' and 'priority' and their values
        """
        return self.config["go_sentences_options"]["evidence_codes"]

    def get_do_evidence_codes(self) -> Dict[str, Dict[str, Union[str, id]]]:
        """get the configured evidence codes for do

        Returns:
            Dict[str, Dict[str, Union[str, id]]]: a dictionary of the configured evidence codes, with evidence code
                names as keys and a dictionary with 'evidence group' and 'priority' and their values
        """
        return self.config["do_sentences_options"]["evidence_codes"]

    def get_do_via_orth_evidence_codes(self) -> Dict[str, Dict[str, Union[str, id]]]:
        """get the configured evidence codes for do via orthology

        Returns:
            Dict[str, Dict[str, Union[str, id]]]: a dictionary of the configured evidence codes, with evidence code
                names as keys and a dictionary with 'evidence group' and 'priority' and their values
        """
        return self.config["do_via_orth_sentences_options"]["evidence_codes"]

    def get_expression_evidence_codes(self) -> Dict[str, Dict[str, Union[str, id]]]:
        """get the configured evidence codes for expression

        Returns:
            Dict[str, Dict[str, Union[str, id]]]: a dictionary of the configured evidence codes, with evidence code
                names as keys and a dictionary with 'evidence group' and 'priority' and their values
        """
        return self.config["expression_sentences_options"]["evidence_codes"]

    def get_go_prepostfix_sentences_map(self) -> Dict[Tuple[str, str, str], Tuple[str, str]]:
        """get the map that links go aspects and evidence groups with their pre- and postfix phrases, including special
        cases

        Returns:
            Dict[Tuple[str, str, str], Tuple[str, str]]: a dictionary that maps a tuple of aspect and group with prefix
                and postfix phrases to be used to build the automatically generated sentences. Special cases are
                transformed and treated as additional groups by appending their id at the end of the group they belong
                to forming thus a new group name. Group priorities for special cases remain equal to their root group
        """
        prepost_map = {(prepost["aspect"], prepost["group"], prepost["qualifier"]): (prepost["prefix"],
                                                                                     prepost["postfix"])
                       for prepost in self.config["go_sentences_options"]["go_prepostfix_sentences_map"]}
        special_cases = self.get_go_prepostfix_special_cases_sent_map()
        for key, scs in special_cases.items():
            for special_case in scs:
                prepost_map[(key[0], key[1] + str(special_case[0]), key[2])] = (special_case[2], special_case[3])
        return prepost_map

    def get_do_prepostfix_sentences_map(self) -> Dict[Tuple[str, str, str], Tuple[str, str]]:
        """get the map that links do aspects and evidence groups with their pre- and postfix phrases

        Returns:
            Dict[Tuple[str, str, str], Tuple[str, str]]: a dictionary that maps a tuple of aspect and group with prefix
                and postfix phrases to be used to build the automatically generated sentences
        """
        prepost_map = {(prepost["aspect"], prepost["group"], prepost["qualifier"]): (prepost["prefix"],
                                                                                     prepost["postfix"])
                       for prepost in self.config["do_sentences_options"]["do_prepostfix_sentences_map"]}
        return prepost_map

    def get_do_via_orth_prepostfix_sentences_map(self) -> Dict[Tuple[str, str, str], Tuple[str, str]]:
        """get the map that links do via orthology aspects and evidence groups with their pre- and postfix phrases

        Returns:
            Dict[Tuple[str, str, str], Tuple[str, str]]: a dictionary that maps a tuple of aspect and group with prefix
                and postfix phrases to be used to build the automatically generated sentences
        """
        prepost_map = {(prepost["aspect"], prepost["group"], prepost["qualifier"]): (prepost["prefix"],
                                                                                     prepost["postfix"])
                       for prepost in self.config["do_via_orth_sentences_options"]["do_prepostfix_sentences_map"]}
        return prepost_map

    def get_expression_prepostfix_sentences_map(self) -> Dict[Tuple[str, str, str], Tuple[str, str]]:
        """get the map that links expression aspects and evidence groups with their pre- and postfix phrases

        Returns:
            Dict[Tuple[str, str, str], Tuple[str, str]]: a dictionary that maps a tuple of aspect and group with prefix
                and postfix phrases to be used to build the automatically generated sentences
        """
        prepost_map = {(prepost["aspect"], prepost["group"], prepost["qualifier"]): (prepost["prefix"],
                                                                                     prepost["postfix"])
                       for prepost in self.config["expression_sentences_options"]
                       ["expression_prepostfix_sentences_map"]}
        return prepost_map

    def get_do_prepostfix_sentences_map_humans(self) -> Dict[Tuple[str, str, str], Tuple[str, str]]:
        """get the map that links do aspects and evidence groups with their pre- and postfix phrases

        Returns:
            Dict[Tuple[str, str, str], Tuple[str, str]]: a dictionary that maps a tuple of aspect and group with prefix
                and postfix phrases to be used to build the automatically generated sentences
        """
        prepost_map = None
        if "do_prepostfix_sentences_map_humans" in self.config["do_sentences_options"]:
            prepost_map = {(prepost["aspect"], prepost["group"], prepost["qualifier"]): (prepost["prefix"],
                                                                                         prepost["postfix"])
                           for prepost in self.config["do_sentences_options"]["do_prepostfix_sentences_map_humans"]}
        return prepost_map

    def get_do_via_orth_prepostfix_sentences_map_humans(self) -> Dict[Tuple[str, str, str], Tuple[str, str]]:
        """get the map that links do aspects and evidence groups with their pre- and postfix phrases

        Returns:
            Dict[Tuple[str, str, str], Tuple[str, str]]: a dictionary that maps a tuple of aspect and group with prefix
                and postfix phrases to be used to build the automatically generated sentences
        """
        prepost_map = None
        if "do_prepostfix_sentences_map_humans" in self.config["do_sentences_options"]:
            prepost_map = {(prepost["aspect"], prepost["group"], prepost["qualifier"]): (prepost["prefix"],
                                                                                         prepost["postfix"])
                           for prepost in self.config["do_via_orth_sentences_options"][
                               "do_prepostfix_sentences_map_humans"]}
        return prepost_map

    def get_go_prepostfix_special_cases_sent_map(self) -> Dict[Tuple[str, str, str], Tuple[int, str, str, str]]:
        """get a map of pre- and postfix phrases for special cases

        Returns:
            Dict[Tuple[str, str, str], Tuple[int, str, str, str]]: a map between aspect, group and qualifier, and
            special cases properties
        """
        return {(prepost["aspect"], prepost["group"], prepost["qualifier"]): [(sp_case["id"], sp_case["match_regex"],
                                                                               sp_case["prefix"], sp_case["postfix"])
                                                                              for sp_case in prepost["special_cases"]]
                for prepost in self.config["go_sentences_options"]["go_prepostfix_sentences_map"] if
                prepost["special_cases"]}

    def get_go_annotations_priority(self) -> List[str]:
        """get the priority list for evidence codes

        Returns:
            List[str]: a list of evidence codes, sorted by priority. The first element has the highest priority
        """
        return [key for key, priority in sorted([(key, ec["priority"]) for key, ec in
                                                 self.get_go_evidence_codes().items()], key=lambda x: x[1])]

    def get_do_annotations_priority(self) -> List[str]:
        """get the priority list for do evidence codes

        Returns:
            List[str]: a list of evidence codes, sorted by priority. The first element has the highest priority
        """
        return [key for key, priority in sorted([(key, ec["priority"]) for key, ec in
                                                 self.get_do_evidence_codes().items()], key=lambda x: x[1])]

    def get_do_via_orth_annotations_priority(self) -> List[str]:
        """get the priority list for do evidence codes

        Returns:
            List[str]: a list of evidence codes, sorted by priority. The first element has the highest priority
        """
        return [key for key, priority in sorted([(key, ec["priority"]) for key, ec in
                                                 self.get_do_via_orth_evidence_codes().items()], key=lambda x: x[1])]

    def get_expression_annotations_priority(self) -> List[str]:
        """get the priority list for expression evidence codes

        Returns:
            List[str]: a list of evidence codes, sorted by priority. The first element has the highest priority
        """
        return [key for key, priority in sorted([(key, ec["priority"]) for key, ec in
                                                 self.get_expression_evidence_codes().items()], key=lambda x: x[1])]

    def get_go_evidence_groups_priority_list(self) -> List[str]:
        """get the priority list for go evidence groups

        Returns:
            List[str]: the priority list for evidence groups
        """
        return [group for group, p in sorted([(g, p) for g, p in
                                              self.config["go_sentences_options"]["group_priority"].items()],
                                             key=lambda x: x[1])]

    def get_do_evidence_groups_priority_list(self) -> List[str]:
        """get the priority list for do evidence groups

        Returns:
            List[str]: the priority list for evidence groups
        """
        return [group for group, p in sorted([(g, p) for g, p in
                                              self.config["do_sentences_options"]["group_priority"].items()],
                                             key=lambda x: x[1])]

    def get_do_via_orth_evidence_groups_priority_list(self) -> List[str]:
        """get the priority list for do evidence groups

        Returns:
            List[str]: the priority list for evidence groups
        """
        return [group for group, p in sorted([(g, p) for g, p in
                                              self.config["do_via_orth_sentences_options"]["group_priority"].items()],
                                             key=lambda x: x[1])]

    def get_expression_evidence_groups_priority_list(self) -> List[str]:
        """get the priority list for expression evidence groups

        Returns:
            List[str]: the priority list for evidence groups
        """
        return [group for group, p in sorted([(g, p) for g, p in
                                              self.config["expression_sentences_options"]["group_priority"].items()],
                                             key=lambda x: x[1])]

    def get_go_evidence_codes_groups_map(self) -> Dict[str, str]:
        """get the map between evidence codes and evidence groups for go

        Returns:
            Dict[str, str]: the map between codes and groups
        """
        return {name: evidence["group"] for name, evidence in self.get_go_evidence_codes().items()}

    def get_do_evidence_codes_groups_map(self) -> Dict[str, str]:
        """get the map between evidence codes and evidence groups for do

        Returns:
            Dict[str, str]: the map between codes and groups
        """
        return {name: evidence["group"] for name, evidence in self.get_do_evidence_codes().items()}

    def get_do_via_orth_evidence_codes_groups_map(self) -> Dict[str, str]:
        """get the map between evidence codes and evidence groups for do

        Returns:
            Dict[str, str]: the map between codes and groups
        """
        return {name: evidence["group"] for name, evidence in self.get_do_via_orth_evidence_codes().items()}

    def get_expression_evidence_codes_groups_map(self) -> Dict[str, str]:
        """get the map between evidence codes and evidence groups for expression

        Returns:
            Dict[str, str]: the map between codes and groups
        """
        return {name: evidence["group"] for name, evidence in self.get_expression_evidence_codes().items()}

    def get_go_terms_exclusion_list(self) -> List[str]:
        """get the list of go terms to exclude from the gene descriptions

        Returns:
            List[str]: the exclusion list
        """
        return self.config["go_sentences_options"]["exclude_terms"]

    def get_do_terms_exclusion_list(self) -> List[str]:
        """get the list of go terms to exclude from the gene descriptions

        Returns:
            List[str]: the exclusion list
        """
        return self.config["do_sentences_options"]["exclude_terms"]

    def get_do_via_orth_terms_exclusion_list(self) -> List[str]:
        """get the list of go terms to exclude from the gene descriptions

        Returns:
            List[str]: the exclusion list
        """
        return self.config["do_via_orth_sentences_options"]["exclude_terms"]

    def get_expression_terms_exclusion_list(self) -> List[str]:
        """get the list of expression terms to exclude from the gene descriptions

        Returns:
            List[str]: the exclusion list
        """
        return self.config["expression_sentences_options"]["exclude_terms"] \
            if "exclude_terms" in self.config["expression_sentences_options"] else None

    def get_raw_file_sources(self, data_fetcher: str) -> str:
        """get the url pointing to the raw files source

        Args:
            data_fetcher (str): the data fetcher to use (e.g., agr_data_fetcher, wb_data_fetcher)
        Returns:
            str: the location of the root directory of the raw files
        """
        return self.config[data_fetcher]["raw_files_source"]

    def get_release(self, data_fetcher):
        """get the release code related to the input data to use

        Returns:
            str: the release code
        """
        return self.config[data_fetcher]["release"]

    def get_agr_mod_property(self, mod: str, property_name: str):
        """get a specific property for the defined organism is AGR

        Args:
            mod (str): the name of the mod
            property_name (str): the name of the property to retrieve
        Returns:
            str: the value of the requested property
        """
        return self.config["agr_data_fetcher"]["organisms"][mod][property_name]

    def get_agr_organisms_to_process(self) -> List[str]:
        """get the list of organisms to process

        Returns:
            List[str]: the list of organisms to process
        """
        return self.config["agr_data_fetcher"]["organisms_to_process"]

    def get_wb_organisms_to_process(self) -> List[str]:
        """get the list of organisms to process for WormBase

        Returns:
            List[str]: the list of organisms to process
        """
        return self.config["wb_data_fetcher"]["organisms_to_process"]

    def get_go_rename_terms(self) -> Dict[str, str]:
        """get the regexp to rename GO terms and their replacement strings

        Returns:
            Dict[str, str]: a map of replacements for GO terms
        """
        return self.config["go_sentences_options"]["rename_terms"]

    def get_expression_rename_terms(self) -> Dict[str, str]:
        """get the regexp to rename expression terms and their replacement strings

        Returns:
            Dict[str, str]: a map of replacements for expression terms
        """
        return self.config["expression_sentences_options"]["rename_terms"]

    def get_go_remove_parents_if_children_are_present(self) -> bool:
        """get the value of the option to remove parent terms from sentences if children are present in the term set

        Returns:
            bool: the value of the option
        """
        return self.config["go_sentences_options"]["remove_parents_if_children_are_present"]

    def get_go_remove_children_if_parent_is_present(self) -> bool:
        """get the value of the option to remove child terms from sentences if parent is present in the term set

        Returns:
            bool: the value of the option
        """
        return self.config["go_sentences_options"]["remove_children_if_parent_is_present"]

    def get_do_remove_parents_if_children_are_present(self) -> bool:
        """get the value of the option to remove parent terms from sentences if children are present in the term set

        Returns:
            bool: the value of the option
        """
        return self.config["do_sentences_options"]["remove_parents_if_children_are_present"]

    def get_do_remove_children_if_parent_is_present(self) -> bool:
        """get the value of the option to remove child terms from sentences if parent is present in the term set

        Returns:
            bool: the value of the option
        """
        return self.config["do_sentences_options"]["remove_children_if_parent_is_present"]

    def get_do_via_orth_remove_parents_if_children_are_present(self) -> bool:
        """get the value of the option to remove parent terms from sentences if children are present in the term set

        Returns:
            bool: the value of the option
        """
        return self.config["do_via_orth_sentences_options"]["remove_parents_if_children_are_present"]

    def get_do_via_orth_remove_children_if_parent_is_present(self) -> bool:
        """get the value of the option to remove child terms from sentences if parent is present in the term set

        Returns:
            bool: the value of the option
        """
        return self.config["do_via_orth_sentences_options"]["remove_children_if_parent_is_present"]

    def get_expression_remove_parents_if_children_are_present(self) -> bool:
        """get the value of the option to remove parent terms from sentences if children are present in the term set

        Returns:
            bool: the value of the option
        """
        return self.config["expression_sentences_options"]["remove_parents_if_children_are_present"]

    def get_expression_remove_children_if_parent_is_present(self) -> bool:
        """get the value of the option to remove child terms from sentences if parent is present in the term set

        Returns:
            bool: the value of the option
        """
        return self.config["expression_sentences_options"]["remove_children_if_parent_is_present"]

    def get_go_trim_terms_by_common_ancestors(self) -> bool:
        """get the value of the option to trim terms by common ancestors

        Returns:
            bool: the value of the option
        """
        return self.config["go_sentences_options"]["trim_terms_by_common_ancestors"]

    def get_go_trim_min_num_terms(self) -> int:
        """get the threshold value that indicates the minimum number of terms per go aspect for which trimming has to
        be applied

        Returns:
            int: the value of the option
        """
        return self.config["go_sentences_options"]["trim_if_more_than_terms"]

    def get_go_max_num_terms(self) -> int:
        """get the maximum number of terms per go aspect to be displayed in the final description

        Returns:
            int: the value of the option
        """
        return self.config["go_sentences_options"]["max_num_terms"]

    def get_go_trim_min_distance_from_root(self):
        """get the minimum distance from root in the GO ontology to be considered while looking for common ancestors
        between terms

        :return: the distances for all go aspects
        :rtype: Dict[str, int]
        """
        return self.config["go_sentences_options"]["trim_min_distance_from_root"]

    def get_go_truncate_others_aggregation_word(self) -> str:
        """get the generic word used to indicate that one or more terms have been omitted from the sentence, e.g.,
        'several'

        Returns:
           str: the aggregation word
        """
        return self.config["go_sentences_options"]["go_truncate_others_aggregation_word"]

    def get_go_truncate_others_terms(self) -> Dict[str, str]:
        """get the specific words used for each aspect to indicate that one or more terms have been omitted from the
        sentence, e.g., 'functions'

        Returns:
            Dict[str, int]: a dictionary containing one word for each aspect
        """
        return self.config["go_sentences_options"]["go_truncate_others_terms"]

    def get_do_trim_min_num_terms(self) -> int:
        """get the threshold value that indicates the minimum number of terms per go aspect for which trimming has to
        be applied

        Returns:
            int: the value of the option
        """
        return self.config["do_sentences_options"]["trim_if_more_than_terms"]

    def get_do_max_num_terms(self) -> int:
        """get the maximum number of terms to be displayed in the final description

        Returns:
            int: the value of the option
        """
        return self.config["do_sentences_options"]["max_num_terms"]

    def get_do_trim_min_distance_from_root(self):
        """get the minimum distance from root in the GO ontology to be considered while looking for common ancestors
        between terms

        Returns:
            Dict[str, int]: the distances for all go aspects
        """
        return self.config["do_sentences_options"]["trim_min_distance_from_root"]

    def get_do_truncate_others_aggregation_word(self) -> str:
        """get the generic word used to indicate that one or more terms have been omitted from the sentence, e.g.,
        'several'

        Returns:
            str: the aggregation word
        """
        return self.config["do_sentences_options"]["do_truncate_others_aggregation_word"]

    def get_do_truncate_others_terms(self) -> Dict[str, str]:
        """get the specific words used for each aspect to indicate that one or more terms have been omitted from the
        sentence

        Returns:
            str: a dictionary containing one word for each aspect
        """
        return self.config["do_sentences_options"]["do_truncate_others_terms"]

    def get_do_via_orth_trim_min_num_terms(self) -> int:
        """get the threshold value that indicates the minimum number of terms per do via orth aspect for which trimming
        has to be applied

        Returns:
            int: the value of the option
        """
        return self.config["do_via_orth_sentences_options"]["trim_if_more_than_terms"]

    def get_do_via_orth_max_num_terms(self) -> int:
        """get the maximum number of terms to be displayed in the final description

        Returns:
            int: the value of the option
        """
        return self.config["do_via_orth_sentences_options"]["max_num_terms"]

    def get_do_via_orth_trim_min_distance_from_root(self):
        """get the minimum distance from root in the GO ontology to be considered while looking for common ancestors
        between terms

        Returns:
            Dict[str, int]: the distances for all go aspects
        """
        return self.config["do_via_orth_sentences_options"]["trim_min_distance_from_root"]

    def get_do_via_orth_truncate_others_aggregation_word(self) -> str:
        """get the generic word used to indicate that one or more terms have been omitted from the sentence, e.g.,
        'several'

        Returns:
            str: the aggregation word
        """
        return self.config["do_via_orth_sentences_options"]["do_truncate_others_aggregation_word"]

    def get_do_via_orth_truncate_others_terms(self) -> Dict[str, str]:
        """get the specific words used for each aspect to indicate that one or more terms have been omitted from the
        sentence

        Returns:
            str: a dictionary containing one word for each aspect
        """
        return self.config["do_via_orth_sentences_options"]["do_truncate_others_terms"]

    def get_expression_trim_min_num_terms(self) -> int:
        """get the threshold value that indicates the minimum number of terms per go aspect for which trimming has to
        be applied

        Returns:
            int: the value of the option
        """
        return self.config["expression_sentences_options"]["trim_if_more_than_terms"]

    def get_expression_trim_min_distance_from_root(self):
        """get the minimum distance from root in the GO ontology to be considered while looking for common ancestors
        between terms

        Returns:
            Dict[str, int]: the distances for all go aspects
        """
        return self.config["expression_sentences_options"]["trim_min_distance_from_root"]

    def get_expression_max_num_terms(self) -> int:
        """get the maximum number of terms to be displayed in the final description

        Returns:
            int: the value of the option
        """
        return self.config["expression_sentences_options"]["max_num_terms"]

    def get_expression_truncate_others_aggregation_word(self) -> str:
        """get the generic word used to indicate that one or more terms have been omitted from the sentence, e.g.,
        'several'

        Returns:
            str: the aggregation word
        """
        return self.config["expression_sentences_options"]["expression_truncate_others_aggregation_word"]

    def get_expression_truncate_others_terms(self) -> Dict[str, str]:
        """get the specific words used for each aspect to indicate that one or more terms have been omitted from the
        sentence

        Returns:
            str: a dictionary containing one word for each aspect
        """
        return self.config["expression_sentences_options"]["expression_truncate_others_terms"]

    def get_genedesc_writer(self):
        """get the type of writer

        Returns:
            str: the type of writer
        """
        return self.config["generic_genedesc_writer"]["genedesc_writer"]

    def get_genedesc_output_dir(self, genedesc_writer: str):
        return self.config[genedesc_writer + "_options"]["output_dir"]

    def get_ortholog_species(self):
        return

    def get_wb_human_orthologs_go_ontology(self):
        return self.config["wb_data_fetcher"]["agr_go_ontology"]

    def get_wb_human_orthologs_go_associations(self):
        return self.config["wb_data_fetcher"]["agr_human_go_associations"]

