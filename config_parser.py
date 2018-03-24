import yaml
from typing import Dict, Union, Tuple, List


class GenedescConfigParser(object):
    def __init__(self, file_path):
        with open(file_path) as conf_file:
            self.config = yaml.load(conf_file)

    def get_data_fetcher(self) -> str:
        """get the data fetcher type from the configuration file

        :return: the data fetcher name
        :rtype: str
        """
        return self.config["generic_data_fetcher"]["data_fetcher"]

    def get_cache_location(self) -> str:
        """get the root location of the cache files

        :return: the path to the cache location
        :rtype: str
        """
        return self.config["generic_data_fetcher"]["cache_location"]

    def get_wb_species(self) -> Dict[str, Dict[str, str]]:
        """get species list for WormBase data fetcher

        :return: the species dictionary, with species name as key and a dictionary with project_id as value
        :rtype: Dict[str, Dict[str, str]]
        """
        return self.config["wb_data_fetcher"]["organisms"]

    def get_go_evidence_codes(self) -> Dict[str, Dict[str, Union[str, id]]]:
        """get the configured evidence codes

        :return: a dictionary of the configured evidence codes, with evidence code names as keys and a dictionary with
            'evidence group' and 'priority' and their values
        :rtype: Dict[str, Dict[str, Union[str, id]]]
        """
        return self.config["go_sentences_options"]["evidence_codes"]

    def get_go_prepostfix_sentences_map(self) -> Dict[Tuple[str, str], Tuple[str, str]]:
        """get the map that links go aspects and evidence groups with their pre- and postfix phrases, including special
        cases

        :return: a dictionary that maps a tuple of aspect and group with prefix and postfix phrases to be used to
            build the automatically generated sentences. Special cases are transformed and treated as additional groups
            by appending their id at the end of the group they belong to forming thus a new group name. Group priorities
            for special cases remain equal to their root group
        :rtype: Dict[Tuple[str, str], Tuple[str, str]]
        """
        prepost_map = {(prepost["aspect"], prepost["group"]): (prepost["prefix"], prepost["postfix"])
                       for prepost in self.config["go_sentences_options"]["go_prepostfix_sentences_map"]}
        special_cases = self.get_go_prepostfix_special_cases_sent_map()
        for key, scs in special_cases.items():
            for special_case in scs:
                prepost_map[(key[0], key[1] + str(special_case[0]))] = (special_case[2], special_case[3])
        return prepost_map

    def get_go_prepostfix_special_cases_sent_map(self) -> Dict[Tuple[str, str], Tuple[int, str, str, str]]:
        """get a map of pre- and postfix phrases for special cases

        :return: a map between aspect and group, and special cases properties
        :rtype: Dict[Tuple[str, str], Tuple[int, str, str, str]]
        """
        return {(prepost["aspect"], prepost["group"]): [(sp_case["id"], sp_case["match_regex"], sp_case["prefix"],
                                                         sp_case["postfix"]) for sp_case in
                                                        prepost["special_cases"]] for prepost in
                self.config["go_sentences_options"]["go_prepostfix_sentences_map"] if prepost["special_cases"]}

    def get_go_annotations_priority(self) -> List[str]:
        """get the priority list for evidence codes

        :return: a list of evidence codes, sorted by priority. The first element has the highest priority
        :rtype: List[str]"""
        return [key for key, priority in sorted([(key, ec["priority"]) for key, ec in
                                                 self.get_go_evidence_codes().items()], key=lambda x: x[1])]

    def get_evidence_groups_priority_list(self) -> List[str]:
        """get the priority list for evidence groups

        :return: the priority list for evidence groups
        :rtype: List[str]
        """
        return [group for group, p in sorted([(g, p) for g, p in
                                              self.config["go_sentences_options"]["group_priority"].items()],
                                             key=lambda x: x[1])]

    def get_evidence_codes_groups_map(self) -> Dict[str, str]:
        """get the map between evidence codes and evidence groups

        :return: the map between codes and groups
        :rtype: Dict[str, str]
        """
        return {name: evidence["group"] for name, evidence in self.get_go_evidence_codes().items()}

    def get_go_terms_exclusion_list(self) -> List[str]:
        """get the list of go terms to exclude from the gene descriptions

        :return: the exclusion list
        :rtype: List[str]
        """
        return self.config["go_sentences_options"]["exclude_terms"]

    def get_raw_file_sources(self, data_fetcher: str) -> str:
        """get the url pointing to the raw files source

        :param data_fetcher: the data fetcher to use (e.g., agr_data_fetcher, wb_data_fetcher)
        :type data_fetcher: str
        :return: the location of the root directory of the raw files
        :rtype: str
        """
        return self.config[data_fetcher]["raw_files_source"]

    def get_chebi_file_source(self) -> str:
        """get the url pointing to the chebi file

        :return: the location of the chebi file
        :rtype: str
        """
        return self.config["generic_data_fetcher"]["chebi_file_source"]

    def get_release(self, data_fetcher):
        """get the release code related to the input data to use

        :return: the release code
        :rtype: str
        """
        return self.config[data_fetcher]["release"]

    def get_agr_mod_property(self, mod: str, property_name: str):
        """get a specific property for the defined organism is AGR

        :param mod: the name of the mod
        :type mod: str
        :param property_name: the name of the property to retrieve
        :type property_name: str
        :return: the value of the requested property
        :rtype: str
        """
        return self.config["agr_data_fetcher"]["organisms"][mod][property_name]

    def get_agr_organisms_to_process(self) -> List[str]:
        """get the list of organisms to process

        :return: the list of organisms to process
        :rtype: List[str]
        """
        return self.config["agr_data_fetcher"]["organisms_to_process"]

    def get_wb_organisms_to_process(self) -> List[str]:
        """get the list of organisms to process for WormBase

        :return: the list of organisms to process
        :rtype: List[str]
        """
        return self.config["wb_data_fetcher"]["organisms_to_process"]

    def get_go_rename_terms(self) -> Dict[str, str]:
        """get the regexp to rename GO terms and their replacement strings

        :return: a map of replacements for GO terms
        :rtype: Dict[str, str]
        """
        return self.config["go_sentences_options"]["rename_terms"]

    def get_go_remove_parents_if_children_are_present(self) -> bool:
        return self.config["go_sentences_options"]["remove_parents_if_children_are_present"]

    def get_go_merge_terms_by_common_ancestors(self) -> bool:
        return self.config["go_sentences_options"]["merge_terms_by_common_ancestors"]

    def get_go_merge_min_num_terms(self) -> int:
        return self.config["go_sentences_options"]["merge_if_more_than_terms"]

    def get_go_merge_min_distance_from_root(self):
        return self.config["go_sentences_options"]["merge_min_distance_from_root"]

    def get_genedesc_writer(self):
        return self.config["generic_genedesc_writer"]["genedesc_writer"]

    def get_genedesc_output_dir(self, genedesc_writer: str):
        return self.config[genedesc_writer + "_options"]["output_dir"]


