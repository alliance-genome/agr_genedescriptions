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


def get_gene_data(provider: str):
    """Get genes from the A-team database."""
    session = create_ateam_db_session()
    try:
        sql_query = text("""
        SELECT
            be.primaryexternalid geneId,
            slota.displaytext geneSymbol
        FROM
            biologicalentity be JOIN slotannotation slota ON be.id = slota.singlegene_id
            JOIN organization org ON be.dataprovider_id = org.id
        WHERE
            slota.obsolete = false
        AND
            slota.slotannotationtype = 'GeneSymbolSlotAnnotation'
        AND
            org.abbreviation = :provider;
        """)
        rows = session.execute(sql_query, {'provider': provider}).fetchall()
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
        SELECT
            otc.curie parentCurie,
            otc.name parentName,
            otc.namespace parentType,
            otc.obsolete parentIsObsolete,
            otp.curie childCurie,
            otp.name childName,
            otp.namespace childType,
            otp.obsolete childIsObsolete,
            'IS_A' relType
        FROM
            ontologyterm otp JOIN ontologyterm_isa_parent_children otpc ON otp.id = otpc.isaparents_id
                             JOIN ontologyterm otc ON otpc.isachildren_id = otc.id
        WHERE
            otp.curie LIKE :curieprefix;
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
