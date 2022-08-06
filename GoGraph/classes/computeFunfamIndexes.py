import traceback
from datetime import datetime

import pymysql


class FunfamIndexes:
    def __init__(self):
        import os
        self.phdPrefix = os.getenv('PHD_DB_PREFIX')

    def compute_index_neo4j(self, indexes, term_annotation, similarity_threshold, output):
        # evidenceCodes = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'TAS', 'IEP', 'IC']
        evidenceCodes = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'TAS', 'IEP', 'IC', 'ISS', ' ISO', ' ISA', ' ISM', ' IGC', ' IBA', ' IBD', ' IKR', ' IRD', ' RCA']

        from utilities.neo4jconnectionpool import Neo4JConnectionPool
        cp = Neo4JConnectionPool()

        graph = None

        try:
            graph = cp.get_connection(none_if_unable=True)

            if graph is None:
                raise Exception("Unable to get a Neo4J Connection")

            start_time = datetime.now()

            count = 0
            countFunFams = 0

            import re
            superFamRe = re.compile("^([0-9\.]+)\.FF")

            files = {}
            for index in indexes:
                files[index] = open("/tmp/FF-{}-Index.txt".format(index), "w")

            cypherQuery = """MATCH (pg:Protein)-[cb:`classified by`]->(d:CathFunFam) WHERE d.fullname = "{ff}"
                        WITH collect(DISTINCT pg) as ffprotein, count(DISTINCT pg) as ffcount
                    """

            if term_annotation == "ALLPROTEINS":
                cypherQuery += """MATCH (ptonly:Protein)-[ab:`annotated by`]->(t:Term) WHERE t.unique_id = "{gt}"
                        AND ab.evidence in {e}
                        WITH count(DISTINCT ptonly) as tcount, ffcount
                    MATCH (pt:Protein)-[ab:`annotated by`]->(t:Term) WHERE t.unique_id = "{gt}" AND pt IN ffprotein
                        AND ab.evidence in {e}
                        WITH count(DISTINCT pt) as cntIntersection, ffcount, tcount
                    """
            elif term_annotation == "SUPERFAMPROTEINS":
                cypherQuery += """MATCH (c:CathFunFam)<-[cb:`classified by`]-(ptonly:Protein)-[ab:`annotated by`]->(t:Term)
                            WHERE t.unique_id = "{gt}" AND c.fullname =~ "{sf}.*"
                            AND ab.evidence in {e}
                        WITH count(DISTINCT ptonly) as tcount, ffcount
                    MATCH (c:CathFunFam)<-[cb:`classified by`]-(pt:Protein)-[ab:`annotated by`]->(t:Term)
                        WHERE t.unique_id = "{gt}" AND c.fullname = "{ff}"
                        AND ab.evidence in {e}
                        WITH count(DISTINCT pt) as cntIntersection, ffcount, tcount
                    """
            elif term_annotation == "FUNFAMPROTEINS":
                cypherQuery += """MATCH (c:CathFunFam)<-[cb:`classified by`]-(ptonly:Protein)-[ab:`annotated by`]->(t:Term)
                            WHERE t.unique_id = "{gt}" AND c.fullname = "{ff}"
                            AND ab.evidence in {e}
                        WITH count(DISTINCT ptonly) as tcount, count(DISTINCT ptonly) as cntIntersection, ffcount
                    """
            elif term_annotation == "SIMILARFUNFAMPROTEINS":
                cypherQuery += """MATCH x=(rff:CathFunFam)<-[r:`is related to`]-(c:CathFunFam)
                          WHERE (c.fullname = "{ff}") AND r.score > {threshold}
                        WITH COLLECT(DISTINCT rff.fullname) AS related1, ffcount
                        MATCH x=(rff:CathFunFam)<-[r:`is related to`]-(c:CathFunFam)
                          WHERE (rff.fullname = "{ff}") AND r.score > {threshold}
                        WITH COLLECT(DISTINCT c.fullname) + related1 AS related, ffcount
                        OPTIONAL MATCH (c:CathFunFam)<-[cb:`classified by`]-(p:Protein)-[ab:`annotated by`]-(t:Term)
                              WHERE (c.fullname in related) AND t.unique_id="{gt}"
                              AND ab.evidence in {e}
                            WITH count(DISTINCT p) AS tcount, ffcount, COLLECT(DISTINCT p) AS relatedProteins
                    OPTIONAL MATCH (c:CathFunFam)<-[cb:`classified by`]-(ptonly:Protein)-[ab:`annotated by`]->(t:Term)
                            WHERE t.unique_id = "{gt}" AND c.fullname = "{ff}" AND ptonly in relatedProteins
                            AND ab.evidence in {e}
                            WITH count(ptonly) as cntIntersection, ffcount, tcount
                """

            cypherQuery += """
            RETURN tcount, ffcount, cntIntersection"""

            for funfam in graph.find(label="CathFunFam"):
                # output.append(funfam["fullname"])
                countFunFams += 1

                filled_query = ""

                terms = graph.run(
                    """MATCH (t:Term)<-[ab:`annotated by`]-(p:Protein)-[r:`classified by`]->(d:CathFunFam)                    
                    WHERE d.fullname = {ff} 
                    AND ab.evidence in {e}
                    RETURN DISTINCT t.unique_id""",
                    ff=funfam["fullname"], e=evidenceCodes)

                for term in terms:
                    count += 1

                    if term_annotation == "ALLPROTEINS" or term_annotation == "FUNFAMPROTEINS":
                        filled_query = cypherQuery.format(ff=funfam["fullname"], gt=term[0], e=evidenceCodes)
                    elif term_annotation == "SUPERFAMPROTEINS":
                        superFam = superFamRe.match(funfam["fullname"]).group(1)
                        filled_query = cypherQuery.format(ff=funfam["fullname"], gt=term[0], sf=superFam,
                                                          e=evidenceCodes)
                    elif term_annotation == "SIMILARFUNFAMPROTEINS":
                        filled_query = cypherQuery.format(ff=funfam["fullname"], gt=term[0],
                                                          threshold=similarity_threshold, e=evidenceCodes)

                    if "JACCARD" in indexes:
                        filled_query += """, ((cntIntersection * 1.0) / (ffcount + tcount - cntIntersection)) as jaccard"""

                        # jaccardIndex = graph.run(jaccardQuery)
                        #
                        # for calc in jaccardIndex:
                        #     jaccardOut.write("{}|{}|{:.4f}\n".format(funfam["fullname"], term[0], calc["jaccard"]))
                        #     # print("{}|{}|{:.4f}" .format(funfam["fullname"], term[0]["unique_id"], calc["jaccard"]))

                    if "SORENSEN" in indexes:
                        filled_query += """, ((cntIntersection * 2.0) / (ffcount + tcount)) as sorensen"""

                        # sorensenIndex = graph.run(sorensenQuery)
                        #
                        # for calc in sorensenIndex:
                        #     sorensenOut.write("{}|{}|{:.4f}\n".format(funfam["fullname"], term[0], calc["sorensen"]))
                        #     # print("{}|{}|{:.4f}" .format(funfam["fullname"], term[0]["unique_id"], calc["sorensen"]))

                    if "OVERLAP" in indexes:
                        # Formula for min: http://goo.gl/3i7m4c
                        filled_query += """, CASE WHEN tcount = 0 THEN 0 WHEN ffcount = 0 THEN 0 ELSE
                                        ( (cntIntersection * 1.0) / (0.5 * (ffcount + tcount - abs(ffcount - tcount)))) END
                                        as overlap"""

                        # overlapIndex = graph.run(overlapQuery)
                        #
                        # for calc in overlapIndex:
                        #     overlapOut.write("{}|{}|{:.4f}\n".format(funfam["fullname"], term[0], calc["overlap"]))
                        #     # print("{}|{}|{:.4f}" .format(funfam["fullname"], term[0]["unique_id"], calc["overlap"]))

                    # print(filled_query)
                    index_calculation = graph.run(filled_query)

                    for calc in index_calculation:
                        if "JACCARD" in indexes:
                            files["JACCARD"].write(
                                "{}|{}|{:.4f}\n".format(funfam["fullname"], term[0], calc["jaccard"]))
                        if "SORENSEN" in indexes:
                            files["SORENSEN"].write(
                                "{}|{}|{:.4f}\n".format(funfam["fullname"], term[0], calc["sorensen"]))
                        if "OVERLAP" in indexes:
                            files["OVERLAP"].write(
                                "{}|{}|{:.4f}\n".format(funfam["fullname"], term[0], calc["overlap"]))
                            # print("{}|{}|{:.4f}" .format(funfam["fullname"], term[0]["unique_id"], calc["jaccard"]))

                    if (count % 100) == 0:
                        timeDiff = (datetime.now() - start_time)
                        if "JACCARD" in indexes:
                            files["JACCARD"].flush()
                        if "SORENSEN" in indexes:
                            files["SORENSEN"].flush()
                        if "OVERLAP" in indexes:
                            files["OVERLAP"].flush()
                        print("Processed {} terms from {} funfams so far. Time taken is: {} seconds".format(count,
                                                                                                            countFunFams,
                                                                                                            timeDiff.seconds))

                    output.append("Total funfams retrieved is %d." % count)
        except OSError as err:
            print("OS error: {0}".format(err))
            output.append("OS error: {0}".format(err))
        except:
            import sys
            print("Unexpected error:", sys.exc_info()[0])
            output.append("Unexpected error:", sys.exc_info()[0])
        finally:
            if "JACCARD" in indexes:
                files["JACCARD"].close()
            if "SORENSEN" in indexes:
                files["SORENSEN"].close()
            if "OVERLAP" in indexes:
                files["OVERLAP"].close()
            if graph != None:
                cp.close_connection(graph)

        return None

    def compute_index(self, indexes, term_annotation, similarity_threshold, output):

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
                        fn_calls += ",Calculate_Jaccard_CMP_Same_Funfam(a.cathfunfamilyfull_id, a.go_term, 0)  AS Jaccard"
                    elif index.upper() == 'SORENSEN':
                        fn_calls += ",Calculate_Sorensen_CMP_Same_Funfam(a.cathfunfamilyfull_id, a.go_term, 0) AS Sorensen"
                    elif index.upper() == 'OVERLAP':
                        fn_calls += ",Calculate_Overlap_CMP_Same_Funfam(a.cathfunfamilyfull_id, a.go_term, 0)  AS Overlap"
                    else:
                        fn_calls += ",NULL"

                sqlTruncateTable = f"TRUNCATE TABLE {self.phdPrefix}.cathfunfamiliessimilarity_samefunfam_similarity_allEC"

                cursor.execute(sqlTruncateTable)

                sqlToExecute = f"""
                    INSERT INTO {self.phdPrefix}.cathfunfamiliessimilarity_samefunfam_similarity_allEC (cathfunfam, goterm, jaccard, sorensen, overlap)
                      (
                        SELECT
                          a.cathfunfamilyfull_id
                          , a.go_term
                          {fn_calls}                         
                        FROM (
                               SELECT DISTINCT
                                 pga.go_term,
                                 pcff.cathfunfamilyfull_id
                               FROM GoGraph_proteincathfunfamily pcff
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

    def compute_index_enriched(self, indexes, term_annotation, similarity_threshold, output):

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
                protein_go_terms_dict = dict()
                go_terms_proteins_dict = dict()
                go_terms_parents_dict = dict()

                sql_fetch_cath_funfams = f"SELECT cff.cathfunfamilyfull_id FROM {self.phdPrefix}.GoGraph_cathfunfamilies cff"
                sql_fetch_protein_annotations = """
                    
                """

                rows = cursor.fetchall(sql_fetch_cath_funfams)

                for row in rows:
                    print(row)
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
