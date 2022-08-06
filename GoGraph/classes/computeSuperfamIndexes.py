import os

from datetime import datetime
import pymysql

import traceback


class SuperfamilyIndexes:
    def __init__(self):
        import os
        self.phdPrefix = os.getenv('PHD_DB_PREFIX')

    def compute_index_neo4j(self, indexes, output):

        # evidenceCodes = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'TAS', 'IEP', 'IC']
        evidenceCodes = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'TAS', 'IEP', 'IC', 'ISS', ' ISO', ' ISA', ' ISM', ' IGC', ' IBA', ' IBD', ' IKR', ' IRD', ' RCA']

        try:
            # authenticate("localhost", "neo4j", "neo4j-password")
            # graph = Graph(host="localhost", user="neo4j", password="neo4j-password", bolt=True)
            # http.socket_timeout = 9999999

            graph_ips = [ip.strip() for ip in os.getenv('GRAPH_IP').split(',')]

            from utilities.neo4jconnectionpool import Neo4JConnectionPool
            cp = Neo4JConnectionPool(graph_ips, os.getenv('GRAPH_USER'), os.getenv('GRAPH_PASSWORD'), max_connections=3)

            graph = cp.get_connection()

            startTime = datetime.now()

            text_file = {}

            for indexUsed in indexes:
                output.append("Starting Computation of {} Index".format(indexUsed))

                text_file[indexUsed] = open("/tmp/SF-{}-Index.txt".format(indexUsed.capitalize()), "w")

                count = 0
                countCathFams = 0

            for cathfam in graph.find(label="CathFam"):
                countCathFams += 1
                # output.append(funfam["fullname"])

                terms = graph.run(
                    """MATCH (t:Term)<-[ab:`annotated by`]-(p:Protein)-[r:`classified by`]->(d:CathFam)
                    WHERE d.name = {f}
                    AND ab.evidence in {e} 
                    RETURN distinct t""",
                    f=cathfam["name"], e=evidenceCodes)

                for term in terms:
                    count += 1

                    if "JACCARD" in indexes:
                        jaccardIndex = graph.run(
                            """MATCH (pg:Protein)-[cb:`classified by`]->(d:CathFam) WHERE d.name = "{f}"
                            WITH collect(pg) as ffprotein, count(DISTINCT pg) as ffcount
                            MATCH (pt:Protein)-[ab:`annotated by`]->(t:Term) WHERE t.unique_id = "{gt}" AND pt IN ffprotein
                            AND ab.evidence in {e}
                            WITH count(pt) as cntIntersection, ffcount
                            MATCH (ptonly:Protein)-[ab:`annotated by`]->(t:Term) WHERE t.unique_id = "{gt}"
                            AND ab.evidence in {e}
                            WITH count(ptonly) as tcount, cntIntersection, ffcount
                            RETURN ((cntIntersection * 1.0) / (ffcount + tcount - cntIntersection)) as jaccard""".format(
                                f=cathfam["name"], gt=term[0]["unique_id"], e=evidenceCodes))

                        for calc in jaccardIndex:
                            text_file["JACCARD"].write("{}|{}|{:.4f}\n".format(cathfam["name"], term[0]["unique_id"], calc["jaccard"]))
                            # print("{}|{}|{:.4f}" .format(funfam["fullname"], term[0]["unique_id"], calc["jaccard"]))

                    if "SORENSEN" in indexes:
                        sorensenIndex = graph.run(
                            """MATCH (pg:Protein)-[cb:`classified by`]->(d:CathFam) WHERE d.name = "{f}"
                            WITH collect(pg) as ffprotein, count(pg) as ffcount
                            MATCH (pt:Protein)-[ab:`annotated by`]->(t:Term) WHERE t.unique_id = "{gt}" AND pt IN ffprotein
                            AND ab.evidence in {e}
                            WITH count(pt) as cntIntersection, ffcount
                            MATCH (ptonly:Protein)-[ab:`annotated by`]->(t:Term) WHERE t.unique_id = "{gt}"
                            AND ab.evidence in {e}
                            WITH count(ptonly) as tcount, cntIntersection, ffcount
                            RETURN ((cntIntersection * 2.0) / (ffcount + tcount)) as sorensen""".format(
                                f=cathfam["name"], gt=term[0]["unique_id"], e=evidenceCodes))

                        for calc in sorensenIndex:
                            text_file["SORENSEN"].write(
                                "{}|{}|{:.4f}\n".format(cathfam["name"], term[0]["unique_id"], calc["sorensen"]))
                            # print("{}|{}|{:.4f}" .format(funfam["fullname"], term[0]["unique_id"], calc["sorensen"]))

                    if "OVERLAP"  in indexes:
                        overlapIndex = graph.run(
                            """MATCH (pg:Protein)-[cb:`classified by`]->(d:CathFam) WHERE d.name = "{f}"
                            WITH collect(pg) as ffprotein, count(pg) as ffcount
                            MATCH (pt:Protein)-[ab:`annotated by`]->(t:Term) WHERE t.unique_id = "{gt}" AND pt IN ffprotein
                            AND ab.evidence in {e}
                            WITH count(pt) as cntIntersection, ffcount
                            MATCH (ptonly:Protein)-[ab:`annotated by`]->(t:Term) WHERE t.unique_id = "{gt}"
                            AND ab.evidence in {e}
                            WITH count(ptonly) as tcount, cntIntersection, ffcount
                            RETURN ( (cntIntersection * 1.0) / (0.5 * (ffcount + tcount - abs(ffcount - tcount))) ) as overlap""".format(
                                f=cathfam["name"], gt=term[0]["unique_id"], e=evidenceCodes))

                        for calc in overlapIndex:
                            text_file["OVERLAP"].write("{}|{}|{:.4f}\n".format(cathfam["name"], term[0]["unique_id"], calc["overlap"]))
                            # print("{}|{}|{:.4f}" .format(funfam["fullname"], term[0]["unique_id"], calc["overlap"]))

                    if (count % 1000) == 0:
                        timeDiff = (datetime.now() - startTime)
                        for indexUsed in indexes:
                            text_file[indexUsed].flush()
                        print(
                            "Processed {} terms from {} fams so far. Time taken is: {} seconds".format(count,
                                                                                                       countCathFams,
                                                                                                       timeDiff.seconds))

            output.append("Total funfams scored is {}.".format(count))
        except OSError as err:
            print("OS error: {0}".format(err))
        except:
            import sys
            print("Unexpected error:", sys.exc_info()[0])

        return None

    def compute_index(self, indexes, output):
        cnx = None
        cursor = None
        fn_calls = ''

        if len(indexes) == 0:
            raise Exception("No indexes selected")

        from utilities.mysqlconnectionpool import MySQLConnectionPool
        cnx = MySQLConnectionPool.get_instance().get_connection()

        try:
            cursor = cnx.cursor()

            for index in indexes:
                if index.upper() == 'JACCARD':
                    fn_calls += ",Calculate_Jaccard_SuperFamily(a.cathfamily_id, a.go_term, 0)  AS Jaccard"
                elif index.upper() == 'SORENSEN':
                    fn_calls += ",Calculate_Sorensen_SuperFamily(a.cathfamily_id, a.go_term, 0) AS Sorensen"
                elif index.upper() == 'OVERLAP':
                    fn_calls += ",Calculate_Overlap_SuperFamily(a.cathfamily_id, a.go_term, 0)  AS Overlap"
                else:
                    fn_calls += ",NULL"

            sqlTruncateTable = f"TRUNCATE TABLE {self.phdPrefix}.cathfamiliessimilarity"

            cursor.execute(sqlTruncateTable)

            sqlToExecute = f"""
                INSERT INTO {self.phdPrefix}.cathfamiliessimilarity (cathfam, goterm, jaccard, sorensen, overlap)
                  (
                    SELECT
                      a.cathfamily_id
                      , a.go_term
                      {fn_calls}                         
                    FROM (
                           SELECT DISTINCT
                             pga.go_term,
                             pcff.cathfamily_id
                           FROM GoGraph_proteincathfamily pcff
                             JOIN GoGraph_proteingoannotation pga ON pcff.protein_name = pga.protein_name
                            WHERE pcff.meets_inclusion_threshold = 1
                            AND pga.evidence_code IN (SELECT evidence_code FROM evidence_codes)
                         ) a
                  );
            """

            cursor.execute(sqlToExecute)
            cnx.commit()

            output.append("Indexes calculated and MySQL Tables are populated")
        except (pymysql.ProgrammingError, pymysql.DatabaseError, pymysql.DataError, pymysql.IntegrityError, pymysql.NotSupportedError, pymysql.OperationalError,
                pymysql.MySQLError) as pe:
            output.append(f"Something is wrong: {str(pe)}")
            output.append(f"Exception Stack Trace: {traceback.format_exc()}")
        finally:
            if cursor is not None:
                cursor.close()
            if cnx is not None:
                MySQLConnectionPool.get_instance().close_connection(cnx)

        return None
