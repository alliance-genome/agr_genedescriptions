from genedescriptions.commons import Gene, Module
from genedescriptions.config_parser import GenedescConfigParser
from genedescriptions.data_manager import DataManager
from genedescriptions.descriptions_generator import OntologySentenceGenerator
from genedescriptions.gene_description import GeneDescription


def set_gene_ontology_module(dm: DataManager, conf_parser: GenedescConfigParser, gene_desc: GeneDescription,
                             gene: Gene):
    go_sent_generator_exp = OntologySentenceGenerator(gene_id=gene.id, module=Module.GO, data_manager=dm,
                                                      config=conf_parser, limit_to_group="EXPERIMENTAL")
    go_sent_generator = OntologySentenceGenerator(gene_id=gene.id, module=Module.GO, data_manager=dm,
                                                  config=conf_parser)
    contributes_to_module_sentences = go_sent_generator.get_module_sentences(
        config=conf_parser, aspect='F', qualifier='contributes_to', merge_groups_with_same_prefix=True,
        keep_only_best_group=True)
    if contributes_to_module_sentences.contains_sentences():
        func_module_sentences = go_sent_generator_exp.get_module_sentences(
            config=conf_parser, aspect='F', merge_groups_with_same_prefix=True, keep_only_best_group=True)
        gene_desc.set_or_extend_module_description_and_final_stats(module_sentences=func_module_sentences,
                                                                   module=Module.GO_FUNCTION)
    else:
        func_module_sentences = go_sent_generator.get_module_sentences(
            config=conf_parser, aspect='F', merge_groups_with_same_prefix=True, keep_only_best_group=True)
        gene_desc.set_or_extend_module_description_and_final_stats(
            module_sentences=func_module_sentences, module=Module.GO_FUNCTION)
    gene_desc.set_or_extend_module_description_and_final_stats(
        module_sentences=contributes_to_module_sentences, module=Module.GO_FUNCTION)
    proc_module_sentences = go_sent_generator.get_module_sentences(
        config=conf_parser, aspect='P', merge_groups_with_same_prefix=True, keep_only_best_group=True)
    gene_desc.set_or_extend_module_description_and_final_stats(
        module_sentences=proc_module_sentences, module=Module.GO_PROCESS)
    colocalizes_with_module_sentences = go_sent_generator.get_module_sentences(
        config=conf_parser, aspect='C', qualifier='colocalizes_with', merge_groups_with_same_prefix=True,
        keep_only_best_group=True)
    if colocalizes_with_module_sentences.contains_sentences():
        comp_module_sentences = go_sent_generator_exp.get_module_sentences(
            config=conf_parser, aspect='C', merge_groups_with_same_prefix=True, keep_only_best_group=True)
        gene_desc.set_or_extend_module_description_and_final_stats(
            module_sentences=comp_module_sentences, module=Module.GO_COMPONENT)
    else:
        comp_module_sentences = go_sent_generator.get_module_sentences(
            config=conf_parser, aspect='C', merge_groups_with_same_prefix=True, keep_only_best_group=True)
        gene_desc.set_or_extend_module_description_and_final_stats(
            module_sentences=comp_module_sentences, module=Module.GO_COMPONENT)
    gene_desc.set_or_extend_module_description_and_final_stats(module_sentences=colocalizes_with_module_sentences,
                                                               module=Module.GO_COMPONENT)
    gene_desc.set_initial_stats(module=Module.GO_FUNCTION, sentence_generator=go_sent_generator,
                                sentence_generator_exp_only=go_sent_generator_exp)
    gene_desc.set_initial_stats(module=Module.GO_PROCESS, sentence_generator=go_sent_generator,
                                sentence_generator_exp_only=go_sent_generator_exp)
    gene_desc.set_initial_stats(module=Module.GO_COMPONENT, sentence_generator=go_sent_generator,
                                sentence_generator_exp_only=go_sent_generator_exp)


def set_disease_module(df: DataManager, conf_parser: GenedescConfigParser, gene_desc: GeneDescription, gene: Gene):
    do_sentence_exp_generator = OntologySentenceGenerator(gene_id=gene.id,
                                                          module=Module.DO_EXPERIMENTAL, data_manager=df,
                                                          config=conf_parser, limit_to_group="EXPERIMENTAL")
    disease_exp_module_sentences = do_sentence_exp_generator.get_module_sentences(
        config=conf_parser, aspect='D', merge_groups_with_same_prefix=True, keep_only_best_group=False)
    gene_desc.set_or_extend_module_description_and_final_stats(module=Module.DO_EXPERIMENTAL,
                                                               module_sentences=disease_exp_module_sentences)
    do_sentence_bio_generator = OntologySentenceGenerator(gene_id=gene.id, module=Module.DO_BIOMARKER,
                                                          data_manager=df, config=conf_parser,
                                                          limit_to_group="BIOMARKER")
    disease_bio_module_sentences = do_sentence_bio_generator.get_module_sentences(
        config=conf_parser, aspect='D', merge_groups_with_same_prefix=True, keep_only_best_group=False)
    gene_desc.set_or_extend_module_description_and_final_stats(module=Module.DO_BIOMARKER,
                                                               module_sentences=disease_bio_module_sentences)
    do_via_orth_sentence_generator = OntologySentenceGenerator(
        gene_id=gene.id, module=Module.DO_ORTHOLOGY, data_manager=df, config=conf_parser)
    disease_via_orth_module_sentences = do_via_orth_sentence_generator.get_module_sentences(
        config=conf_parser, aspect='D', merge_groups_with_same_prefix=True, keep_only_best_group=False)
    gene_desc.set_or_extend_module_description_and_final_stats(module=Module.DO_ORTHOLOGY,
                                                               module_sentences=disease_via_orth_module_sentences)
    gene_desc.set_initial_stats(module=Module.DO_EXPERIMENTAL, sentence_generator=do_sentence_exp_generator)
    gene_desc.set_initial_stats(module=Module.DO_BIOMARKER, sentence_generator=do_sentence_bio_generator)
    gene_desc.set_initial_stats(module=Module.DO_ORTHOLOGY, sentence_generator=do_via_orth_sentence_generator)

