import yaml

from enum import Enum
from typing import Dict, List
from genedescriptions.commons import Module


class ConfigModuleProperty(Enum):
    RENAME_TERMS = 1
    DEL_PARENTS_IF_CHILD = 2
    DEL_CHILDREN_IF_PARENT = 3
    APPLY_TRIMMING = 4
    MAX_NUM_TERMS_BEFORE_TRIMMING = 5
    MAX_NUM_TERMS_IN_SENTENCE = 6
    DISTANCE_FROM_ROOT = 7
    CUTOFF_SEVERAL_WORD = 8
    CUTOFF_SEVERAL_CATEGORY_WORD = 9
    EXCLUDE_TERMS = 10
    ADD_MULTIPLE_TO_COMMON_ANCEST = 11
    RENAME_CELL = 12
    REMOVE_OVERLAP = 13


class GenedescConfigParser(object):
    def __init__(self, file_path):
        with open(file_path) as conf_file:
            self.config = yaml.load(conf_file)

    def get_module_property(self, module: Module, prop: ConfigModuleProperty):
        module_name = self._get_module_name(module)
        property_name = self._get_module_property_name(prop)
        if module_name in self.config and property_name in self.config[module_name]:
            return self.config[module_name][property_name]
        else:
            return None

    @staticmethod
    def _get_module_name(module: Module):
        module_name = ""
        if module == Module.GO:
            module_name = "go_sentences_options"
        elif module == Module.DO_EXP_AND_BIO:
            module_name = "do_sentences_options"
        elif module == Module.DO_ORTHOLOGY:
            module_name = "do_via_orth_sentences_options"
        elif module == Module.EXPRESSION:
            module_name = "expression_sentences_options"
        return module_name

    @staticmethod
    def _get_module_property_name(prop: ConfigModuleProperty):
        property_name = ""
        if prop == ConfigModuleProperty.RENAME_TERMS:
            property_name = "rename_terms"
        if prop == ConfigModuleProperty.EXCLUDE_TERMS:
            property_name = "exclude_terms"
        elif prop == ConfigModuleProperty.DEL_PARENTS_IF_CHILD:
            property_name = "remove_parents_if_children_are_present"
        elif prop == ConfigModuleProperty.DEL_CHILDREN_IF_PARENT:
            property_name = "remove_children_if_parent_is_present"
        elif prop == ConfigModuleProperty.APPLY_TRIMMING:
            property_name = "trim_terms_by_common_ancestors"
        elif prop == ConfigModuleProperty.MAX_NUM_TERMS_BEFORE_TRIMMING:
            property_name = "trim_if_more_than_terms"
        elif prop == ConfigModuleProperty.MAX_NUM_TERMS_IN_SENTENCE:
            property_name = "max_num_terms"
        elif prop == ConfigModuleProperty.DISTANCE_FROM_ROOT:
            property_name = "trim_min_distance_from_root"
        elif prop == ConfigModuleProperty.CUTOFF_SEVERAL_WORD:
            property_name = "truncate_others_aggregation_word"
        elif prop == ConfigModuleProperty.CUTOFF_SEVERAL_CATEGORY_WORD:
            property_name = "truncate_others_terms"
        elif prop == ConfigModuleProperty.ADD_MULTIPLE_TO_COMMON_ANCEST:
            property_name = "add_multiple_if_covers_more_children"
        elif prop == ConfigModuleProperty.RENAME_CELL:
            property_name = "rename_cell"
        elif prop == ConfigModuleProperty.REMOVE_OVERLAP:
            property_name = "remove_overlapped_terms"
        return property_name

    def get_prepostfix_sentence_map(self, module: Module, special_cases_only: bool = False, humans: bool = False):
        module_name = self._get_module_name(module)
        if special_cases_only:
            return {(prepost["aspect"], prepost["group"], prepost["qualifier"]): [
                (sp_case["id"], sp_case["match_regex"], sp_case["prefix"], sp_case["postfix"])
                for sp_case in prepost["special_cases"]]
                for prepost in self.config[module_name]["prepostfix_sentences_map"] if
                "special_cases" in prepost and prepost["special_cases"]}
        else:
            prepost_map = {(prepost["aspect"], prepost["group"], prepost["qualifier"]): (
                prepost["prefix"], prepost["postfix"]) for prepost in self.config[module_name][
                "prepostfix_sentences_map_humans" if humans else "prepostfix_sentences_map"]}
            special_cases_only = self.get_prepostfix_sentence_map(module=module, special_cases_only=True, humans=humans)
            for key, scs in special_cases_only.items():
                for special_case in scs:
                    prepost_map[(key[0], key[1] + str(special_case[0]), key[2])] = (special_case[2], special_case[3])
            return prepost_map

    def get_annotations_priority(self, module: Module) -> List[str]:
        module_name = self._get_module_name(module)
        return [key for key, priority in sorted(
            [(key, ec["priority"]) for key, ec in self.config[module_name]["evidence_codes"].items()],
            key=lambda x: x[1])]

    def get_evidence_groups_priority_list(self, module: Module) -> List[str]:
        module_name = self._get_module_name(module)
        return [group for group, p in sorted([(g, p) for g, p in self.config[module_name]["group_priority"].items()],
                                             key=lambda x: x[1])]

    def get_evidence_codes_groups_map(self, module: Module) -> Dict[str, str]:
        module_name = self._get_module_name(module)
        return {name: evidence["group"] for name, evidence in
                self.config[module_name]["evidence_codes"].items()}

    def get_out_dir(self) -> str:
        return self.config["generic"]["output_dir"]

    def get_cache_dir(self) -> str:
        return self.config["generic"]["cache_location"]

    def get_wb_raw_file_sources(self) -> str:
        return self.config["wb_options"]["raw_files_source"]

    def get_wb_release(self):
        return self.config["wb_options"]["release"]

    def get_wb_organisms_to_process(self) -> List[str]:
        return self.config["wb_options"]["organisms_to_process"]

    def get_wb_human_orthologs_go_ontology(self):
        return self.config["wb_options"]["agr_go_ontology"]

    def get_wb_human_orthologs_go_associations(self):
        return self.config["wb_options"]["agr_human_go_associations"]

    def get_wb_organisms_info(self):
        return self.config["wb_options"]["organisms"]
