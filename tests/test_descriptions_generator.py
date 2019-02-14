import logging
import unittest
import os

from ontobio import AssociationSetFactory
from genedescriptions.commons import Module, Gene
from genedescriptions.config_parser import GenedescConfigParser
from genedescriptions.data_manager import DataManager, DataType
from genedescriptions.descriptions_generator import OntologySentenceGenerator
from genedescriptions.gene_description import GeneDescription
from genedescriptions.precanned_modules import set_gene_ontology_module

logger = logging.getLogger("Gene Ontology Module tests")


class TestDescriptionsGenerator(unittest.TestCase):

    def setUp(self):
        logger.info("Starting Ontology Tools tests")
        self.this_dir = os.path.split(__file__)[0]
        self.conf_parser = GenedescConfigParser(os.path.join(self.this_dir, os.path.pardir, "tests", "config_test.yml"))
        self.df = DataManager(do_relations=None, go_relations=["subClassOf", "BFO:0000050"])
        logger.info("Loading go ontology from file")
        logging.basicConfig(filename=None, level="ERROR", format='%(asctime)s - %(name)s - %(levelname)s: %(message)s')
        self.df.load_ontology_from_file(ontology_type=DataType.GO, ontology_url="file://" + os.path.join(
            self.this_dir, "data", "go_gd_test.obo"),
                                        ontology_cache_path=os.path.join(self.this_dir, "cache", "go_gd_test.obo"),
                                        config=self.conf_parser)
        logger.info("Loading go associations from file")
        self.df.load_associations_from_file(associations_type=DataType.GO, associations_url="file://" + os.path.join(
            self.this_dir, "data", "gene_association_1.7.wb.partial"),
                                            associations_cache_path=os.path.join(self.this_dir, "cache",
                                                                                 "gene_association_1.7.wb.partial"),
                                            config=self.conf_parser)
        logging.basicConfig(filename=None, level="INFO", format='%(asctime)s - %(name)s - %(levelname)s: %(message)s')
        logger.info("Loading do ontology from file")
        logging.basicConfig(filename=None, level="ERROR",
                            format='%(asctime)s - %(name)s - %(levelname)s: %(message)s')

    def test_trimming_with_high_priority(self):
        generator = OntologySentenceGenerator(gene_id="WB:WBGene00000912", module=Module.GO,
                                              data_manager=self.df, config=self.conf_parser)
        sentences = generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                   qualifier='', merge_groups_with_same_prefix=True,
                                                   keep_only_best_group=True, high_priority_term_ids=['GO:0007568',
                                                                                                      'GO:1900426'])
        self.assertTrue("several processes" in sentences.get_description())
        self.assertTrue("aging" in sentences.get_description())
        self.assertTrue("positive regulation of defense response to bacterium" in sentences.get_description())
        self.assertTrue("regulation of cellular biosynthetic process" in sentences.get_description())

    def test_generate_sentence_wb(self):
        go_sent_generator = OntologySentenceGenerator(gene_id="WB:WBGene00000018", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue("several processes" in sentences.get_description())
        self.assertTrue("DNA damage checkpoint" in sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="WB:WBGene00000001", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue("several processes" not in sentences.get_description())
        self.assertTrue("dauer larval development, determination of adult lifespan, and insulin receptor "
                        "signaling pathway" in sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="WB:WBGene00002335", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='F',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        print(sentences.get_description())

    def test_generate_sentence_fb(self):
        self.df.load_associations_from_file(associations_type=DataType.GO, associations_url="file://" + os.path.join(
            self.this_dir, "data", "gene_association_1.7.fb.partial"),
                                            associations_cache_path=os.path.join(self.this_dir, "cache",
                                                                                 "gene_association_1.7.fb.partial"),
                                            config=self.conf_parser)
        go_sent_generator = OntologySentenceGenerator(gene_id="FB:FBgn0027655", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in several processes, including axo-dendritic transport, establishment of mitotic '
                        'spindle orientation, and positive regulation of extent of heterochromatin assembly',
                        sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="FB:FBgn0045035", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in several processes, including DNA metabolic process, chromatin organization, '
                        'and negative regulation of cell death', sentences.get_description())

    def test_generate_sentence_human(self):
        self.df.load_associations_from_file(associations_type=DataType.GO, associations_url="file://" + os.path.join(
            self.this_dir, "data", "gene_association_1.7.human.partial"),
                                            associations_cache_path=os.path.join(self.this_dir, "cache",
                                                                                 "gene_association_1.7.human.partial"),
                                            config=self.conf_parser)
        go_sent_generator = OntologySentenceGenerator(gene_id="RGD:HGNC:4851", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in several processes, including positive regulation of autophagy, regulation of '
                        'phosphate metabolic process, and regulation of signal transduction',
                        sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="RGD:HGNC:1884", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in several processes, including anion transport, cellular response to forskolin, '
                        'and positive regulation of transport', sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="RGD:HGNC:795", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in several processes, including cellular protein modification process, negative '
                        'regulation of cell cycle, and regulation of nucleobase-containing compound metabolic process',
                        sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="RGD:HGNC:11291", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in several processes, including cell migration involved in sprouting '
                        'angiogenesis, cellular senescence, and regulation of cellular biosynthetic process',
                        sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="RGD:HGNC:348", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in several processes, including cAMP-mediated signaling, cellular response to '
                        'organic cyclic compound, and regulation of B cell proliferation', sentences.get_description())

    def test_generate_sentence_mgi(self):
        self.df.load_associations_from_file(associations_type=DataType.GO, associations_url="file://" + os.path.join(
            self.this_dir, "data", "gene_association_1.7.mgi.partial"),
                                            associations_cache_path=os.path.join(self.this_dir, "cache",
                                                                                 "gene_association_1.7.mgi.partial"),
                                            config=self.conf_parser)
        go_sent_generator = OntologySentenceGenerator(gene_id="MGI:MGI:96067", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in several processes, including carboxylic acid metabolic process, '
                        'cytoskeleton-dependent intracellular transport, and nervous system development',
                        sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="MGI:MGI:88388", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in several processes, including cellular response to cAMP, ion homeostasis, and '
                        'positive regulation of transport', sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="MGI:MGI:107202", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in several processes, including DNA metabolic process, animal organ development, '
                        'and meiotic nuclear division', sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="MGI:MGI:106658", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in several processes, including actin cytoskeleton organization, animal organ '
                        'development, and morphogenesis of an epithelium', sentences.get_description())

    def test_generate_sentence_zfin(self):
        self.df.load_associations_from_file(associations_type=DataType.GO, associations_url="file://" + os.path.join(
            self.this_dir, "data", "gene_association_1.7.zfin.partial"),
                                            associations_cache_path=os.path.join(self.this_dir, "cache",
                                                                                 "gene_association_1.7.zfin.partial"),
                                            config=self.conf_parser)
        go_sent_generator = OntologySentenceGenerator(gene_id="ZFIN:ZDB-GENE-010608-2", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        print(sentences.get_description())

    def test_generate_sentence_rgd(self):
        self.df.load_associations_from_file(associations_type=DataType.GO, associations_url="file://" + os.path.join(
            self.this_dir, "data", "gene_association_1.7.rgd.partial"),
                                            associations_cache_path=os.path.join(self.this_dir, "cache",
                                                                                 "gene_association_1.7.rgd.partial"),
                                            config=self.conf_parser)
        go_sent_generator = OntologySentenceGenerator(gene_id="RGD:68337", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in several processes, including central nervous system development, mRNA '
                        'transport, and negative regulation of cysteine-type endopeptidase activity',
                        sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="RGD:2332", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in several processes, including animal organ development, anion transport, and '
                        'regulation of cell development', sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="RGD:1593265", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in several processes, including apoptotic process, cellular response to DNA '
                        'damage stimulus, and positive regulation of cellular metabolic process',
                        sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="RGD:1559787", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in several processes, including cellular response to glucose stimulus, '
                        'long-term memory, and negative regulation of cell migration', sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="RGD:61995", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in several processes, including cellular response to hormone stimulus, '
                        'dephosphorylation, and negative regulation of transport', sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="RGD:2074", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in several processes, including cellular response to organic cyclic compound, '
                        'ovarian follicle development, and regulation of transcription, DNA-templated',
                        sentences.get_description())

    def test_generate_sentence_sgd(self):
        self.df.load_associations_from_file(associations_type=DataType.GO, associations_url="file://" + os.path.join(
            self.this_dir, "data", "gene_association_1.7.sgd.partial"),
                                            associations_cache_path=os.path.join(self.this_dir, "cache",
                                                                                 "gene_association_1.7.sgd.partial"),
                                            config=self.conf_parser)
        go_sent_generator = OntologySentenceGenerator(gene_id="SGD:S000004695", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in cellular protein-containing complex assembly and signal peptide processing',
                        sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="SGD:S000004916", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in iron import into the mitochondrion and iron-sulfur cluster assembly',
                        sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="SGD:S000004646", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in DNA replication initiation, regulation of cellular biosynthetic process, and '
                        'regulation of mating type switching', sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="SGD:S000000253", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in regulation of transcription by RNA polymerase II, termination of RNA '
                        'polymerase II transcription, and transcriptional start site selection at RNA polymerase II '
                        'promoter', sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="SGD:S000000364", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in several processes, including cellular protein metabolic process, negative '
                        'regulation of cell cycle, and regulation of nucleobase-containing compound metabolic process',
                        sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="SGD:S000002284", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in several processes, including SCF complex disassembly in response to cadmium '
                        'stress, cellular macromolecule catabolic process, and cellular protein complex disassembly',
                        sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="SGD:S000004603", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in several processes, including exit from mitosis, meiotic nuclear division, and '
                        'regulation of organelle organization', sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="SGD:S000004802", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in several processes, including DNA metabolic process, meiotic chromosome '
                        'segregation, and negative regulation of cell cycle', sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="SGD:S000005707", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in actin cortical patch localization and positive regulation of cellular '
                        'component organization', sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="SGD:S000001596", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in nucleic acid metabolic process and replicative cell aging',
                        sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="SGD:S000004777", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in meiotic heteroduplex formation, meiotic mismatch repair, and reciprocal '
                        'meiotic recombination', sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="SGD:S000006074", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in several processes, including DNA damage checkpoint, DNA metabolic process, '
                        'and deoxyribonucleoside triphosphate biosynthetic process', sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="SGD:S000002678", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in cellular iron ion homeostasis and copper ion export',
                        sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="SGD:S000003487", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in ubiquinone biosynthetic process',
                        sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="SGD:S000000458", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in cellular protein-containing complex assembly and endoplasmic reticulum to '
                        'Golgi vesicle-mediated transport', sentences.get_description())
        go_sent_generator = OntologySentenceGenerator(gene_id="SGD:S000006068", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue('is involved in fatty acid transport', sentences.get_description())
        gene = Gene(id="SGD:S000000482", name="DPB3", dead=False, pseudo=False)
        gene_desc = GeneDescription(gene_id=gene.id, gene_name=gene.name, add_gene_name=True)
        set_gene_ontology_module(dm=self.df, conf_parser=self.conf_parser, gene_desc=gene_desc, gene=gene)
        self.assertTrue('is involved in fatty acid transport', gene_desc.description)

    def test_information_content_sentence_generation(self):
        self.conf_parser.config["go_sentences_options"]["trimming_algorithm"] = "ic"
        go_sent_generator = OntologySentenceGenerator(gene_id="WB:WBGene00000912", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue(sentences.get_description() != "", "Description is empty")
        go_sent_generator = OntologySentenceGenerator(gene_id="WB:WBGene00003931", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue(sentences.get_description() != "", "Description is empty")

    def test_naive_algorithm_sentence_generation(self):
        self.conf_parser.config["go_sentences_options"]["trimming_algorithm"] = "naive"
        go_sent_generator = OntologySentenceGenerator(gene_id="WB:WBGene00000912", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue(sentences.get_description() != "", "Description is empty")
        go_sent_generator = OntologySentenceGenerator(gene_id="WB:WBGene00003931", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue(sentences.get_description() != "", "Description is empty")
        self.conf_parser.config["go_sentences_options"]["trimming_algorithm"] = "lca"
        go_sent_generator = OntologySentenceGenerator(gene_id="WB:WBGene00000912", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue(sentences.get_description() != "", "Description is empty")
        go_sent_generator = OntologySentenceGenerator(gene_id="WB:WBGene00003931", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue(sentences.get_description() != "", "Description is empty")

    def test_blacklist_during_trimming(self):
        associations = [association for subj_associations in self.df.go_associations.associations_by_subj.values() for
                        association in subj_associations]
        associations.append(DataManager.create_annotation_record(source_line="", gene_id="WB:WBGene00003931",
                                                                 gene_symbol="", gene_type="gene", taxon_id="",
                                                                 object_id="GO:0043055", qualifiers="", aspect="P",
                                                                 ecode="EXP", references="", prvdr="WB", date=""))
        associations.append(DataManager.create_annotation_record(source_line="", gene_id="WB:WBGene00003931",
                                                                 gene_symbol="", gene_type="gene", taxon_id="",
                                                                 object_id="GO:0061065", qualifiers="", aspect="P",
                                                                 ecode="EXP", references="", prvdr="WB", date=""))
        associations.append(DataManager.create_annotation_record(source_line="", gene_id="WB:WBGene00003931",
                                                                 gene_symbol="", gene_type="gene", taxon_id="",
                                                                 object_id="GO:0043054", qualifiers="", aspect="P",
                                                                 ecode="EXP", references="", prvdr="WB", date=""))
        associations.append(DataManager.create_annotation_record(source_line="", gene_id="WB:WBGene00003931",
                                                                 gene_symbol="", gene_type="gene", taxon_id="",
                                                                 object_id="GO:0043053", qualifiers="", aspect="P",
                                                                 ecode="EXP", references="", prvdr="WB", date=""))
        self.df.go_associations = AssociationSetFactory().create_from_assocs(assocs=associations,
                                                                             ontology=self.df.go_ontology)
        self.conf_parser.config["go_sentences_options"]["exclude_terms"].append("GO:0040024")
        go_sent_generator = OntologySentenceGenerator(gene_id="WB:WBGene00003931", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue("dauer larval development" not in sentences.get_description(), "Blacklist not working")
        self.conf_parser.config["go_sentences_options"]["trimming_algorithm"] = "ic"
        go_sent_generator = OntologySentenceGenerator(gene_id="WB:WBGene00003931", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue("dauer larval development" not in sentences.get_description(), "Blacklist not working")
        self.conf_parser.config["go_sentences_options"]["trimming_algorithm"] = "lca"
        go_sent_generator = OntologySentenceGenerator(gene_id="WB:WBGene00003931", module=Module.GO,
                                                      data_manager=self.df, config=self.conf_parser)
        sentences = go_sent_generator.get_module_sentences(config=self.conf_parser, aspect='P',
                                                           qualifier='', merge_groups_with_same_prefix=True,
                                                           keep_only_best_group=True)
        self.assertTrue("dauer larval development" not in sentences.get_description(), "Blacklist not working")

    def test_disease_trimming(self):
        self.df.load_ontology_from_file(ontology_type=DataType.DO, ontology_url="file://" + os.path.join(
            self.this_dir, "data", "doid.obo"),
                                        ontology_cache_path=os.path.join(self.this_dir, "cache", "doid.obo"),
                                        config=self.conf_parser)
        associations = [DataManager.create_annotation_record(source_line="", gene_id="RGD:HGNC:7225",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="DOID:0110200", qualifiers="", aspect="D",
                                                             ecode="IAGP", references="", prvdr="RGD", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="RGD:HGNC:7225",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="DOID:0110152", qualifiers="", aspect="D",
                                                             ecode="IAGP", references="", prvdr="RGD", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="RGD:HGNC:7225",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="DOID:0110158", qualifiers="", aspect="D",
                                                             ecode="IAGP", references="", prvdr="RGD", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="RGD:HGNC:7225",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="DOID:0110157", qualifiers="", aspect="D",
                                                             ecode="IAGP", references="", prvdr="RGD", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="RGD:HGNC:7225",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="DOID:0050540", qualifiers="", aspect="D",
                                                             ecode="IAGP", references="", prvdr="RGD", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="RGD:HGNC:7225",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="DOID:0110195", qualifiers="", aspect="D",
                                                             ecode="IAGP", references="", prvdr="RGD", date=""),
                        DataManager.create_annotation_record(source_line="", gene_id="RGD:HGNC:7225",
                                                             gene_symbol="", gene_type="gene", taxon_id="",
                                                             object_id="DOID:10595", qualifiers="", aspect="D",
                                                             ecode="IAGP", references="", prvdr="RGD", date="")]
        self.df.do_associations = AssociationSetFactory().create_from_assocs(assocs=associations,
                                                                             ontology=self.df.do_ontology)
        generator = OntologySentenceGenerator(gene_id="RGD:HGNC:7225", module=Module.DO_EXPERIMENTAL,
                                              data_manager=self.df, config=self.conf_parser, humans=True)
        sentences = generator.get_module_sentences(
            config=self.conf_parser, aspect='D', qualifier='', merge_groups_with_same_prefix=True,
            keep_only_best_group=True)
        self.assertTrue("Charcot-Marie-Tooth disease (multiple)" in sentences.get_description())



