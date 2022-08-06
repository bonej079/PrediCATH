import traceback

import pymysql


class PfamFunfamIndexes:

    def __init__(self):
        import os
        self.phdPrefix = os.getenv("PHD_DB_PREFIX")

    def compute_index(self, indexes, term_annotation, output):

        cnx = None
        cursor = None

        if len(indexes) == 0:
            raise Exception("No indexes selected")

        from utilities.mysqlconnectionpool import MySQLConnectionPool
        cnx = MySQLConnectionPool.get_instance().get_connection()

        try:
            cursor = cnx.cursor()

            if term_annotation == "ALLPROTEINS":
                output.append("ALLPROTEINS is not implemented yet!")
            elif term_annotation == "SUPERFAMPROTEINS":
                output.append("SUPERFAMPROTEINS is not implemented yet!")
            elif term_annotation == "FUNFAMPROTEINS":
                fn_calls = ''

                for index in indexes:
                    if index.upper() == 'JACCARD':
                        fn_calls += ",Calculate_Jaccard_Pfam_CMP_Same_Funfam(a.pfamfunfamilyfull_id, a.go_term, 0)  AS Jaccard"
                    elif index.upper() == 'SORENSEN':
                        fn_calls += ",Calculate_Sorensen_Pfam_CMP_Same_Funfam(a.pfamfunfamilyfull_id, a.go_term, 0) AS Sorensen"
                    elif index.upper() == 'OVERLAP':
                        fn_calls += ",Calculate_Overlap_Pfam_CMP_Same_Funfam(a.pfamfunfamilyfull_id, a.go_term, 0)  AS Overlap"
                    else:
                        fn_calls += ",NULL"

                sqlTruncateTable = f"TRUNCATE TABLE {self.phdPrefix}.pfamfunfamiliessimilarity_samefunfam_similarity_allEC"

                cursor.execute(sqlTruncateTable)

                sqlToExecute = f"""
                    INSERT INTO {self.phdPrefix}.pfamfunfamiliessimilarity_samefunfam_similarity_allEC (pfamfunfam, goterm, jaccard, sorensen, overlap)
                      (
                        SELECT
                          a.pfamfunfamilyfull_id
                          , a.go_term
                          {fn_calls}                         
                        FROM (
                               SELECT DISTINCT
                                 pga.go_term,
                                 pcff.pfamfunfamilyfull_id
                               FROM GoGraph_proteinpfamfunfamily pcff
                                 JOIN GoGraph_proteingoannotation pga ON pcff.protein_name = pga.protein_name
                                WHERE pcff.meets_inclusion_threshold = 1
                                AND pga.evidence_code IN (SELECT evidence_code FROM evidence_codes)
                             ) a
                      );
                """

                cursor.execute(sqlToExecute)
                cnx.commit()
            elif term_annotation == "SIMILARFUNFAMPROTEINS":
                output.append("SIMILARFUNFAMPROTEINS is not implemented yet!")

            output.append("Indexes calculated and MySQL Tables are populated")
        except (pymysql.ProgrammingError, pymysql.DatabaseError, pymysql.DataError, pymysql.IntegrityError,
                pymysql.NotSupportedError, pymysql.OperationalError, pymysql.MySQLError) as pe:
            output.append(f"Something is wrong: {str(pe)}")
            output.append(f"Exception Stack Trace: {traceback.format_exc()}")
        finally:
            if cursor is not None:
                cursor.close()
            if cnx is not None:
                MySQLConnectionPool.get_instance().close_connection(cnx)

        return None
