from os import environ

from sqlalchemy import create_engine, text
from sqlalchemy.orm import sessionmaker


def create_ateam_db_session():
    """Create and return a SQLAlchemy session connected to the A-team database."""
    USER = environ.get('PERSISTENT_STORE_DB_USERNAME', 'unknown')
    PASSWORD = environ.get('PERSISTENT_STORE_DB_PASSWORD', 'unknown')
    SERVER = environ.get('PERSISTENT_STORE_DB_HOST', 'localhost')
    PORT = environ.get('PERSISTENT_STORE_DB_PORT', '5432')
    DB = environ.get('PERSISTENT_STORE_DB_NAME', 'unknown')

    engine_var = 'postgresql://' + USER + ":" + PASSWORD + '@' + SERVER + ':' + PORT + '/' + DB
    engine = create_engine(engine_var)

    SessionClass = sessionmaker(bind=engine, autoflush=False, autocommit=False)
    session = SessionClass()
    return session


def get_gene_data(species_taxon: str):
    """Get genes from the A-team database."""
    session = create_ateam_db_session()
    try:
        sql_query = text("""
        SELECT
            be.primaryexternalid geneId,
            slota.displaytext geneSymbol
        FROM
            biologicalentity be JOIN slotannotation slota ON be.id = slota.singlegene_id
            JOIN ontologyterm taxon ON be.taxon_id = taxon.id
        WHERE
            slota.obsolete = false
        AND
            be.obsolete = false
        AND
            slota.slotannotationtype = 'GeneSymbolSlotAnnotation'
        AND
            taxon.curie = :species_taxon;
        """)
        rows = session.execute(sql_query, {'species_taxon': species_taxon}).fetchall()
        return [{"gene_id": row[0], "gene_symbol": row[1]} for row in rows]
    finally:
        session.close()


def get_expression_annotations(taxon_id: str):
    """Get expression annotations from the A-team database."""
    session = create_ateam_db_session()
    try:
        sql_query = text("""
        SELECT
            be.primaryexternalid geneId,
            slota.displaytext geneSymbol,
            ot.curie anatomyId
        FROM
            geneexpressionannotation gea JOIN expressionpattern ep ON gea.expressionpattern_id = ep.id
                                         JOIN anatomicalsite asi ON ep.whereexpressed_id = asi.id
                                         JOIN ontologyterm ot ON asi.anatomicalstructure_id = ot.id
                                         JOIN gene g ON gea.expressionannotationsubject_id = g.id
                                         JOIN biologicalentity be ON g.id = be.id
                                         JOIN ontologyterm ot_taxon ON be.taxon_id = ot_taxon.id
                                         JOIN slotannotation slota ON g.id = slota.singlegene_id
        WHERE
            slota.obsolete = false
        AND
            be.obsolete = false
        AND
            slota.slotannotationtype = 'GeneSymbolSlotAnnotation'
        AND
            ot.curie <> 'WBbt:0000100'
        AND ot_taxon.curie = :taxon_id;
        """)
        rows = session.execute(sql_query, {'taxon_id': taxon_id}).fetchall()
        return [{"gene_id": row[0], "gene_symbol": row[1], "anatomy_id": row[2]} for row in rows]
    finally:
        session.close()


def get_ontology_pairs(curie_prefix: str):
    session = create_ateam_db_session()
    try:
        sql_query = text("""
        SELECT DISTINCT
            otp.curie parentCurie,
            otp.name parentName,
            otp.namespace parentType,
            otp.obsolete parentIsObsolete,
            otc.curie childCurie,
            otc.name childName,
            otc.namespace childType,
            otc.obsolete childIsObsolete,
            jsonb_array_elements_text(otpc.closuretypes) AS relType
        FROM
            ontologyterm otc JOIN ontologytermclosure otpc ON otc.id = otpc.closuresubject_id
                             JOIN ontologyterm otp ON otpc.closureobject_id = otp.id
        WHERE
            otp.curie LIKE :curieprefix
        AND
            otpc.distance = 1
        AND
            otpc.closuretypes in ('["part_of"]', '["is_a"]')
        """)
        rows = session.execute(sql_query, {'curieprefix': f"{curie_prefix}%"}).fetchall()
        return [
            {
                "parent_curie": row[0],
                "parent_name": row[1],
                "parent_type": row[2],
                "parent_is_obsolete": row[3],
                "child_curie": row[4],
                "child_name": row[5],
                "child_type": row[6],
                "child_is_obsolete": row[7],
                "rel_type": row[8]
            } for row in rows]
    finally:
        session.close()


def get_data_providers():
    """Get data providers from the A-team database."""
    session = create_ateam_db_session()
    try:
        sql_query = text("""
        SELECT
            s.displayName, t.curie
        FROM
            species s
        JOIN
            ontologyterm t ON s.taxon_id = t.id
        WHERE
            s.obsolete = false
        AND
            s.assembly_curie is not null
        """)
        rows = session.execute(sql_query).fetchall()
        return [(row[0], row[1]) for row in rows]
    finally:
        session.close()


def get_disease_annotations(taxon_id: str):
    """Get direct and indirect disease ontology (DO) annotations from the A-team database.
    - Direct: gene -> DO term (via diseaseannotation_gene)
    - Indirect: gene -> allele -> DO term (via allelediseaseannotation_gene, only if gene has a single allele)
    """
    session = create_ateam_db_session()
    try:
        # Direct gene -> DO term annotations
        direct_query = text("""
            SELECT
                be.primaryexternalid AS geneId,
                slota.displaytext geneSymbol,
                ot.curie AS doId,
                rel.name as relationshipType
            FROM
                diseaseannotation da
            JOIN genediseaseannotation gda ON da.id = gda.id
            JOIN gene g ON gda.diseaseannotationsubject_id = g.id
            JOIN biologicalentity be ON g.id = be.id
            JOIN ontologyterm ot ON da.diseaseannotationobject_id = ot.id
            JOIN slotannotation slota ON g.id = slota.singlegene_id
            JOIN vocabularyterm rel ON da.relation_id = rel.id
            WHERE
                da.obsolete = false
            AND da.negated = false
            AND ot.namespace = 'disease_ontology'
            AND be.taxon_id = (SELECT id FROM ontologyterm WHERE curie = :taxon_id)
        """)
        direct_rows = session.execute(direct_query, {"taxon_id": taxon_id}).mappings().all()

        # Indirect gene -> allele -> DO term annotations (only if allele has a single inferred gene)
        # TODO: there are no inferred genes for WB - check that the query is ok
        indirect_query = text("""
            SELECT
                be.primaryexternalid AS geneId,
                slota.displaytext geneSymbol,
                ot.curie AS doId,
                rel.name as relationshipType
            FROM
                allelediseaseannotation ada
            JOIN diseaseannotation da ON ada.id = da.id
            JOIN biologicalentity be ON ada.inferredgene_id = be.id
            JOIN slotannotation slota ON be.id = slota.singlegene_id
            JOIN ontologyterm ot ON da.diseaseannotationobject_id = ot.id
            JOIN vocabularyterm rel ON da.relation_id = rel.id
            WHERE
                da.obsolete = false
            AND da.negated = false
            AND ot.namespace = 'disease_ontology'
            AND be.taxon_id = (SELECT id FROM ontologyterm WHERE curie = :taxon_id)
            AND ada.diseaseannotationsubject_id IN (
                SELECT ada2.diseaseannotationsubject_id
                FROM allelediseaseannotation ada2
                GROUP BY ada2.diseaseannotationsubject_id
                HAVING COUNT(*) = 1
            )
        """)
        indirect_rows = session.execute(indirect_query, {"taxon_id": taxon_id}).mappings().all()

        # Combine and deduplicate
        seen = set()
        results = []
        for row in list(direct_rows) + list(indirect_rows):
            key = (row["geneId"], row["geneSymbol"], row["doId"], row["relationshipType"])
            if key not in seen:
                results.append({
                    "gene_id": row["geneId"],
                    "gene_symbol": row["geneSymbol"],
                    "do_id": row["doId"],
                    "relationship_type": row["relationshipType"]
                })
                seen.add(key)
        return results
    finally:
        session.close()


def get_best_human_orthologs_for_taxon(taxon_curie: str):
    """
    Get best human orthologs for all genes from a given species taxon curie.
    Only consider orthologs predicted by methods that pass the strict filter (strictfilter=true).
    For each gene, return a tuple: (list of best orthologs, bool indicating if any orthologs were excluded).
    """
    session = create_ateam_db_session()
    try:
        sql_query = text("""
        SELECT
            subj_be.primaryexternalid AS gene_id,
            subj_slota.displaytext AS gene_symbol,
            obj_be.primaryexternalid AS ortho_id,
            obj_slota.displaytext AS ortho_symbol,
            obj_full_name_slota.displaytext AS ortho_full_name,
            COUNT(DISTINCT pm.predictionmethodsmatched_id) AS method_count
        FROM genetogeneorthology gto
        JOIN genetogeneorthologygenerated gtog ON gto.id = gtog.id AND gtog.strictfilter = true
        JOIN genetogeneorthologygenerated_predictionmethodsmatched pm ON gtog.id = pm.genetogeneorthologygenerated_id
        JOIN gene subj_gene ON gto.subjectgene_id = subj_gene.id
        JOIN biologicalentity subj_be ON subj_gene.id = subj_be.id
        JOIN slotannotation subj_slota ON subj_gene.id = subj_slota.singlegene_id AND subj_slota.slotannotationtype = 'GeneSymbolSlotAnnotation' AND subj_slota.obsolete = false
        JOIN gene obj_gene ON gto.objectgene_id = obj_gene.id
        JOIN biologicalentity obj_be ON obj_gene.id = obj_be.id
        JOIN slotannotation obj_slota ON obj_gene.id = obj_slota.singlegene_id AND obj_slota.slotannotationtype = 'GeneSymbolSlotAnnotation' AND obj_slota.obsolete = false
        JOIN slotannotation obj_full_name_slota ON obj_gene.id = obj_full_name_slota.singlegene_id AND obj_full_name_slota.slotannotationtype = 'GeneFullNameSlotAnnotation' AND obj_full_name_slota.obsolete = false
        JOIN ontologyterm obj_taxon ON obj_be.taxon_id = obj_taxon.id
        JOIN ontologyterm subj_taxon ON subj_be.taxon_id = subj_taxon.id
        WHERE subj_taxon.curie = :taxon_curie
          AND obj_taxon.curie = 'NCBITaxon:9606'
          AND subj_slota.obsolete = false
          AND obj_slota.obsolete = false
          AND subj_be.obsolete = false
          AND obj_be.obsolete = false
        GROUP BY gto.subjectgene_id, gto.objectgene_id, subj_be.primaryexternalid, subj_slota.displaytext, obj_be.primaryexternalid, obj_slota.displaytext, obj_full_name_slota.displaytext
        """)
        rows = session.execute(sql_query, {'taxon_curie': taxon_curie}).mappings().all()
        from collections import defaultdict
        gene_orthologs = defaultdict(list)
        for row in rows:
            gene_id = row['gene_id']
            ortho_info = [row['ortho_id'], row['ortho_symbol'], row['ortho_full_name']]
            method_count = row['method_count']
            gene_orthologs[gene_id].append((ortho_info, method_count))
        result = {}
        for gene_id, ortho_list in gene_orthologs.items():
            if not ortho_list:
                continue
            max_count = max(x[1] for x in ortho_list)
            best_orthos = [x[0] for x in ortho_list if x[1] == max_count]
            excluded = len(ortho_list) > len(best_orthos)
            result[gene_id] = (best_orthos, excluded)
        return result
    finally:
        session.close()

