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
                g.primaryexternalid AS geneId,
                ot.curie AS doId
            FROM
                diseaseannotation da
            JOIN diseaseannotation_gene dag ON da.id = dag.diseaseannotation_id
            JOIN gene g ON dag.with_id = g.id
            JOIN diseaseannotation_ontologyterm daot ON da.id = daot.diseaseannotation_id
            JOIN ontologyterm ot ON daot.ontologyterm_id = ot.id
            WHERE
                da.obsolete = false
            AND ot.namespace = 'disease_ontology'
            AND g.taxon_id = (SELECT id FROM ontologyterm WHERE curie = :taxon_id)
        """)
        direct_rows = session.execute(direct_query, {"taxon_id": taxon_id}).fetchall()

        # Indirect gene -> allele -> DO term annotations (only if gene has a single allele)
        indirect_query = text("""
            SELECT
                g.primaryexternalid AS geneId,
                ot.curie AS doId
            FROM
                allelediseaseannotation ada
            JOIN allelediseaseannotation_gene adag ON ada.id = adag.allelediseaseannotation_id
            JOIN gene g ON adag.inferredgene_id = g.id
            JOIN allelediseaseannotation_ontologyterm adaot ON ada.id = adaot.allelediseaseannotation_id
            JOIN ontologyterm ot ON adaot.ontologyterm_id = ot.id
            WHERE
                ada.obsolete = false
            AND ot.namespace = 'disease_ontology'
            AND g.taxon_id = (SELECT id FROM ontologyterm WHERE curie = :taxon_id)
            AND (
                SELECT COUNT(*) FROM allelediseaseannotation_gene adag2 WHERE adag2.inferredgene_id = g.id
            ) = 1
        """)
        indirect_rows = session.execute(indirect_query, {"taxon_id": taxon_id}).fetchall()

        # Combine and deduplicate
        seen = set()
        results = []
        for row in list(direct_rows) + list(indirect_rows):
            key = (row["geneId"], row["doId"])
            if key not in seen:
                results.append({"gene_id": row["geneId"], "do_id": row["doId"]})
                seen.add(key)
        return results
    finally:
        session.close()

