from GoGraph.classes.decorators import exception, create_logger
from GoGraph.classes.searches import Searches


class CheckCathData:

    version = ""

    def __init__(self):
        self.version = "1.0a"

    @exception(create_logger())
    def start_check(self, output):
        from utilities.mysqlconnectionpool import MySQLConnectionPool
        cnx = MySQLConnectionPool.get_instance().get_connection()

        try:
            with cnx.cursor() as db_cursor:
                sql_count_superfamilies = """
                  SELECT count(DISTINCT cf.cathfamily_id)
                  FROM GoGraph_cathfamilies cf;
                """
                sql_count_funfams = """
                  SELECT count(DISTINCT cff.cathfunfamilyfull_id)
                  FROM GoGraph_cathfunfamilies cff; 
                """

                sql_count_superfamilies_used_in_scoring = """
                    select count(distinct pcf.cathfamily_id)
                    from GoGraph_proteincathfamily pcf
                """

                sql_count_funfams_used_in_scoring = """
                    select count(distinct pcff.cathfunfamilyfull_id)
                    from GoGraph_proteincathfunfamily pcff;
                """

                output.append("Basic Statistics")
                output.append("================")

                db_cursor.execute(sql_count_superfamilies)
                output.append("Number of super families: {:,}".format(int(db_cursor.fetchone()[0])))

                db_cursor.execute(sql_count_funfams)
                output.append("Number of functional families: {:,}".format(int(db_cursor.fetchone()[0])))

                db_cursor.execute(sql_count_superfamilies_used_in_scoring)
                output.append("Number of super families used in scoring: {:,}".format(int(db_cursor.fetchone()[0])))

                db_cursor.execute(sql_count_funfams_used_in_scoring)
                output.append("Number of functional families used in scoring: {:,}".format(int(db_cursor.fetchone()[0])))

                output.append("Checking the quality of annotations in database")
                output.append("===============================================")

                sql_get_distinct_proteins = """
                    select distinct p.protein_name
                    from GoGraph_proteingoannotation pga
                    join GoGraph_proteins p on pga.protein_name = p.protein_name;
                """

                sql_get_distinct_go_terms_for_protein = """
                    select distinct pga.go_term
                    from GoGraph_proteingoannotation pga
                    join GoGraph_proteins p on pga.protein_name = p.protein_name
                    where p.protein_name = '{}'
                """

                db_cursor.execute(sql_get_distinct_proteins)
                all_proteins = db_cursor.fetchall()
                count = 0

                for protein in all_proteins:
                    search = Searches()
                    expected_go_terms = set(search.quickgo_terms_by_protein_id(protein[0]))
                    # expected_go_terms = expected_go_terms.union(set(search.go_terms_by_protein_id(protein_name)))
                    if 'Gene ontology (GO)\n' in expected_go_terms:
                        expected_go_terms.remove('Gene ontology (GO)\n')
                    output.append("{}: {}".format(protein[0], sorted(expected_go_terms)))

                    db_cursor.execute(sql_get_distinct_go_terms_for_protein)
                    all_proteins = db_cursor.fetchall()

                    count += 1
                    if count > 10:
                        break
        finally:
            MySQLConnectionPool.get_instance().close_connection(cnx)
