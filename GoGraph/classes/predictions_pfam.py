import io
import os
import pickle
import socket
import tempfile
import traceback
from hashlib import blake2b
from math import sqrt

import pymysql
from Bio import SeqIO
from Bio.Seq import Seq
from joblib import Parallel, delayed

from GoGraph.classes.blast import Blast
from GoGraph.classes.geneontology import GeneOntologyUtility
from GoGraph.classes.pfam_funfams import PFAMmer
from GoGraph.classes.searches import Searches


def unwrap_self_predict(arg, **kwarg):
    return Predictions_PFAM.predict(*arg, **kwarg)


class Predictions_PFAM:
    __slots__ = ['threshold', 'inherit_go_terms', 'separate_ontologies', 'indexes_to_use', 'ts', 'evidenceCodes', 'go_tool', 'logger', 'connection', 'cache_connection',
                 'sql_term_probability', 'use_rollback', 'phdPrefix', 'phdCachePrefix']

    def __init__(self, threshold: float, inherit_go_terms: bool = False, separate_ontologies: bool = False, indexes_to_use: list = [], use_rollback=False,
                 log_to_file: bool = True):

        # evidenceCodes = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'TAS', 'IEP', 'IC']
        self.evidenceCodes = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'TAS', 'IEP', 'IC', 'ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'IBA', 'IBD', 'IKR', 'IRD', 'RCA']  # , 'IEA']

        self.go_tool = GeneOntologyUtility()

        self.inherit_go_terms = inherit_go_terms
        self.separate_ontologies = separate_ontologies
        self.indexes_to_use = indexes_to_use
        self.threshold = threshold

        import datetime
        self.ts = datetime.datetime.now().strftime("%Y%m%d%H%M%S")

        tmp_path = os.getenv('PREDICTION_SAVE_PATH_' + socket.gethostname().upper())

        from GoGraph.classes.decorators import create_logger
        import logging
        if log_to_file:
            self.logger = create_logger(os.path.join(tmp_path, 'predictions', 'phd.log'), loglevel=logging.DEBUG)
        else:
            self.logger = create_logger(loglevel=logging.DEBUG)

        self.phdPrefix = os.getenv("PHD_DB_PREFIX")
        self.phdCachePrefix = os.getenv("PHDCACHE_DB_PREFIX")

        self.sql_term_probability = f"""
                SELECT format(ic.probability,2) as "probability"
                FROM {self.phdPrefix}.infocontentcafa ic
                WHERE ic.goterm = %s;
            """

        self.use_rollback = use_rollback

        from utilities.mysqlconnectionpool import MySQLConnectionPool
        pool = MySQLConnectionPool.get_instance()
        self.connection = pool.get_connection()

        self.cache_connection = MySQLConnectionPool().get_instance().get_connection(pool_name="CACHE")

    def __exit__(self, exc_type, exc_val, exc_tb):
        from utilities.mysqlconnectionpool import MySQLConnectionPool
        MySQLConnectionPool.get_instance().close_connection(self.connection)
        MySQLConnectionPool.get_instance().close_connection(self.cache_connection, pool_name="CACHE")

    def search_expected_terms_quickgo(self, protein_id: str, searchable_sequence: str, search: Searches,
                                      only_admit_evidence_codes=bool):
        if only_admit_evidence_codes:
            expected_go_terms = search.quickgo_terms_by_protein_id(protein_id=protein_id,
                                                                   only_admit_evidence_codes=self.evidenceCodes)
        else:
            expected_go_terms = search.quickgo_terms_by_protein_id(protein_id=protein_id)

        if expected_go_terms is not None:
            if 'Gene ontology (GO)\n' in expected_go_terms:
                expected_go_terms.remove('Gene ontology (GO)\n')
            if '' in expected_go_terms and len(expected_go_terms) == 1:
                blast = Blast()
                results = blast.search_by_blast(str(searchable_sequence))
                if results[0] is not None:
                    protein_id = "{}({})".format(protein_id, results[0])
                    if only_admit_evidence_codes:
                        expected_go_terms = search.quickgo_terms_by_protein_id(results[0], self.evidenceCodes)
                    else:
                        expected_go_terms = search.quickgo_terms_by_protein_id(results[0])
                else:
                    expected_go_terms = None
            else:
                self.logger.info("Expected GO Terms for protein: {}".format(expected_go_terms))

        return expected_go_terms

    def search_expected_terms_from_rollback(self, protein_id: str, cnx: pymysql.Connection,
                                            only_admit_evidence_codes=bool):
        expected_go_terms = []

        terms_sql = f"""
            SELECT DISTINCT go_term
            FROM {self.phdPrefix}.uniprot_goa_125 ug125
            WHERE ug125.protein_name='{protein_id}'
        """ if not only_admit_evidence_codes else f"""
            SELECT DISTINCT go_term
            FROM {self.phdPrefix}.uniprot_goa_125 ug125
            WHERE ug125.protein_name='{protein_id}'
            AND ug125.evidence_code in ({str(self.evidenceCodes).replace(']', '').replace('[', '')})
        """

        with cnx.cursor() as cursor:
            cursor.execute(terms_sql)

            for row in cursor:
                expected_go_terms.append(row['go_term'])

        return expected_go_terms

    def get_expected_go_terms(self, protein_id: str, searchable_sequence: Seq, search: Searches):
        cnx = self.connection

        expected_go_terms = set()
        bp_terms = []
        mf_terms = []
        cc_terms = []

        # Step 1C: Get the expected Go Terms
        self.logger.debug("Step 1C")
        # expected_go_terms = set(search.go_terms_by_protein_id(protein_id))
        # expected_go_terms = expected_go_terms.union(set(search.go_terms_by_protein_id(protein_name)))
        if protein_id is not None:
            if self.use_rollback:
                expected_go_terms = self.search_expected_terms_from_rollback(protein_id, cnx)
            else:
                expected_go_terms = self.search_expected_terms_quickgo(protein_id,
                                                                       searchable_sequence,
                                                                       search)

            # Step 1D: Since we have expected terms, check what ontology they belong to
            self.logger.debug("Step 1D")

            if expected_go_terms is not None and len(expected_go_terms) > 0:
                term_ontologies = self.go_tool.get_ontology_for_term_list(expected_go_terms)

                for term in term_ontologies:
                    ttype = term_ontologies[term]
                    if ttype == 'biological_process':
                        bp_terms.append(term)
                    if ttype == 'molecular_function':
                        mf_terms.append(term)
                    if ttype == 'cellular_component':
                        cc_terms.append(term)

        return expected_go_terms, bp_terms, cc_terms, mf_terms

    def step0_identify_sequence(self, sequence: str):
        self.logger.debug("Step 0")

        fasta_handle = io.StringIO(sequence)

        target_id = None
        protein_description = None
        searchable_sequence = None

        for seq_record in SeqIO.parse(fasta_handle, "fasta"):
            target_id = seq_record.id
            protein_description = seq_record.description
            searchable_sequence = seq_record.seq

        return target_id, protein_description, searchable_sequence

    def step1A_protein_name_and_id(self, searchable_sequence: Seq, protein_description: str):
        # Step 1A

        # If the FASTA does not have
        # if protein_id is None:
        self.logger.debug("Blasting Sequence: {}".format(searchable_sequence))
        blast = Blast()
        results = blast.search_by_blast(str(searchable_sequence), verbose=True)

        if results is not None:
            protein_id = results[0]
            protein_name = results[1]
        else:
            protein_id = protein_description.split(' ')[0].strip()
            if '(' in protein_description:
                protein_name = protein_description.split(' ')[1].strip(')').strip('(').strip()
            else:
                protein_name = protein_description[protein_description.find(' '):].strip()

            self.logger.warn("Warning: Cannot get details from Uniprot. Reverted to provided details.")

        # protein_gene_access = search.gene_access_from_protein_id(protein_id)

        self.logger.debug("Processing {}".format(protein_description))
        self.logger.debug(
            "Processing {} (Blasted protein id {})".format(protein_description, protein_id))

        return protein_id, protein_name

    def step1B_pfam_funfamilies(self, pfamer: PFAMmer, protein_description: str, searchable_sequence: Seq,
                                output: list):
        self.logger.debug("Step 1B")
        related_families = pfamer.pfams_search(protein_description,
                                               sequence=str(searchable_sequence),
                                               output=output)
        # , verbose=False, timeout=300, output=output)
        self.logger.info("Related families for protein: {}".format(related_families))

        superfamilies = set()
        funfams = set()

        if related_families is not None:
            for family in related_families:
                superfamily = family[0:family.find('.FF')]
                superfamilies.add(superfamily)
                funfams.add(family)

        self.logger.debug(superfamilies)
        self.logger.debug(funfams)

        if related_families is None or len(related_families) == 0:
            self.logger.error(f"ERROR: Pfamer returned None for: {protein_description}")
            self.logger.error(f"ERROR: Faulty sequence is: \n{str(searchable_sequence)}")

        return superfamilies, funfams

    def predict(self, sequence: str, output: list, cafa_set: str, file_dict: dict, term_file_dict: dict):
        db_cursor = None

        ontologies_list = ['BP', 'CC', 'MF']

        sql_retrieve_protein = f"""
            SELECT * FROM {self.phdCachePrefix}PFAM.pred_! WHERE hashed_sequence = %s
        """

        sql_insert_protein = f"""
                                REPLACE INTO {self.phdCachePrefix}PFAM.pred_!(hashed_sequence, protein_description, protein_id, protein_name, related_sfamilies,
                                related_ffamilies, expected_go_terms, expected_go_terms_bp, expected_go_terms_cc,
                                expected_go_terms_mf) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s);
                            """

        sql_insert_protein_rollback = f"""
                                            REPLACE INTO {self.phdCachePrefix}PFAM.pred_!(hashed_sequence, protein_description, protein_id, protein_name, related_sfamilies,
                                            related_ffamilies, expected_go_terms_rollback, expected_go_terms_rollback_bp, expected_go_terms_rollback_cc,
                                            expected_go_terms_rollback_mf) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s);
                                        """

        sql_update_funfams = f"""
            UPDATE {self.phdCachePrefix}PFAM.pred_! SET related_sfamilies = %s, related_ffamilies = %s WHERE hashed_sequence = %s
        """

        search = Searches()
        pfamer = PFAMmer()

        m = blake2b(digest_size=blake2b.MAX_DIGEST_SIZE, key=b"JBPhD", salt=b"PhD is SALTY")

        m.update(sequence.encode())
        hashed_sequence = m.hexdigest()

        target_id, protein_description, searchable_sequence = self.step0_identify_sequence(sequence)

        from utilities.mysqlconnectionpool import MySQLConnectionPool
        cp = MySQLConnectionPool.get_instance()
        connection = None

        try:
            # with sqlite3.connect(cache_props['NAME'], check_same_thread=True, timeout=60.0) as cache_conn:

            connection = cp.get_connection(pool_name="CACHE")

            with connection.cursor() as cache_cursor:
                cache_cursor.execute(sql_retrieve_protein.replace("!", cafa_set), (hashed_sequence,))

                row = cache_cursor.fetchone()

                # Step 1: Get the necessary data from Cache or Search
                self.logger.debug(f"Step 1 - Processing Set {cafa_set}")
                if row is not None:
                    self.logger.debug("Step 1[DB]")
                    protein_id = row["protein_id"]
                    protein_name = row["protein_name"]
                    superfamilies = pickle.loads(row['related_sfamilies'])
                    funfams = pickle.loads(row['related_ffamilies'])
                    if self.use_rollback:
                        expected_go_terms = pickle.loads(row['expected_go_terms_rollback'])
                        bp_terms = pickle.loads(row['expected_go_terms_rollback_bp'])
                        cc_terms = pickle.loads(row['expected_go_terms_rollback_cc'])
                        mf_terms = pickle.loads(row['expected_go_terms_rollback_mf'])
                    else:
                        expected_go_terms = pickle.loads(row['expected_go_terms'])
                        bp_terms = pickle.loads(row['expected_go_terms_bp'])
                        cc_terms = pickle.loads(row['expected_go_terms_cc'])
                        mf_terms = pickle.loads(row['expected_go_terms_mf'])

                    self.logger.debug(
                        "Processing {} (Blasted protein id {}) from cache".format(protein_name, protein_id))
                    self.logger.debug(f"Superfamilies: <{superfamilies}>")
                    self.logger.debug(f"Funfams: <{funfams}>")
                    self.logger.debug(f"Expected Go Terms: <{expected_go_terms}>")

                    # Sanity Checks

                    if superfamilies is None or funfams is None or len(superfamilies.union(funfams)) == 0:
                        self.logger.warn(
                            "No PFAM family records found returned None for: {} - {}".format(protein_id,
                                                                                             protein_name))
                        self.logger.warn("Faulty sequence is: \n{}".format(sequence))

                        superfamilies, funfams = self.step1B_pfam_funfamilies(pfamer, protein_description,
                                                                              searchable_sequence, output)

                        if superfamilies is None or funfams is None or len(superfamilies.union(funfams)) == 0:
                            self.logger.warn("Retry to get PFAM Funfams failed. Aborting this prediction.")
                            return
                        else:
                            values = (superfamilies, funfams, hashed_sequence)
                            cache_cursor.execute(sql_update_funfams.replace("!", cafa_set), values)

                    if expected_go_terms is None or len(expected_go_terms) == 0:
                        self.logger.warn(
                            "Finding GO Terms for protein {}({}) failed, hence retrying".format(protein_id,
                                                                                                protein_name))

                        # Try repeating Step 1C - and update the record appropriately

                        expected_go_terms, bp_terms, cc_terms, mf_terms = self.get_expected_go_terms(protein_id,
                                                                                                     str(
                                                                                                         searchable_sequence),
                                                                                                     search)

                        if expected_go_terms is None or len(expected_go_terms) == 0:  # If it is still 0
                            self.logger.warn(
                                "No experimental GO Terms found for protein {}({}), hence skipped".format(protein_id,
                                                                                                          protein_name))
                            return
                else:  # Nothing found in cache
                    # Step 1A: Get updated protein_id and protein_name using BLAST (on SwissProt)
                    protein_id, protein_name = self.step1A_protein_name_and_id(searchable_sequence, protein_description)

                    # Step 1B: Pfamer the sequence
                    superfamilies, funfams = self.step1B_pfam_funfamilies(pfamer, protein_description,
                                                                          searchable_sequence,
                                                                          output)

                    # Step 1C and 1D: Get the expected Go Terms (True Positive Set) and split them into respective
                    # ontology
                    expected_go_terms, bp_terms, cc_terms, mf_terms = self.get_expected_go_terms(protein_id, connection,
                                                                                                 searchable_sequence)

                    # Step 1E: Save the data into the cache
                    self.logger.debug('Step 1E')
                    values = (hashed_sequence, protein_description, protein_id, protein_name,
                              pickle.dumps(superfamilies), pickle.dumps(funfams), pickle.dumps(expected_go_terms),
                              pickle.dumps(bp_terms), pickle.dumps(cc_terms), pickle.dumps(mf_terms))
                    if self.use_rollback:
                        cache_cursor.execute(sql_insert_protein_rollback.replace("!", cafa_set), values)
                    else:
                        cache_cursor.execute(sql_insert_protein.replace("!", cafa_set), values)

                    # Step 1F: Sanity Checks
                    self.logger.debug('Step 1F')

                    if funfams is None or superfamilies is None or len(funfams) == 0 or len(superfamilies) == 0:
                        self.logger.error(
                            "ERROR: Pfamer returned None for: {} - {}".format(protein_id, protein_name))
                        self.logger.error("ERROR: Faulty sequence is: \n{}".format(searchable_sequence))

                        return

                    if expected_go_terms is None:
                        self.logger.error(
                            "ERROR: Finding GO Terms for protein {}({}) failed, hence skipped".format(protein_id,
                                                                                                      protein_name))
                        return

                    if len(expected_go_terms) == 0:
                        self.logger.error(
                            "ERROR: No experimental GO Terms found for protein {}({}), hence skipped".format(
                                protein_id,
                                protein_name))
                        return

            # Step 2: Perform the prediction
            self.logger.debug("Step 2")

            with connection.cursor() as db_cursor:
                if "dcgo" in self.indexes_to_use:
                    term_file_funfams = term_file_dict[("dcgo", "FF")]
                    pred_file_funfams = file_dict[("dcgo", "FF")]

                    if self.separate_ontologies:
                        for ontology in ontologies_list:
                            try:
                                if ontology == 'BP':
                                    expected_go_terms = set(bp_terms)
                                elif ontology == 'MF':
                                    expected_go_terms = set(mf_terms)
                                else:
                                    expected_go_terms = set(cc_terms)

                                if len(expected_go_terms) == 0:  # There might not be terms in the ontology
                                    continue

                                results = self.perform_prediction_dcgo(db_cursor, expected_go_terms,
                                                                       # str(funfams).replace("{", "(")
                                                                       # .replace("}", ")"),
                                                                       list(funfams),
                                                                       ontology, term_file_funfams,
                                                                       protein_name, protein_id)

                                if results is None:
                                    continue

                                # save to result file

                                pred_file_funfams.write("{},{},{},".format(protein_name, protein_id, ontology))
                                pred_file_funfams.write(",".join(str(v) for v in results))
                                pred_file_funfams.write("\n")

                                term_file_funfams.write("\n")

                                self.logger.info("Prediction Results: {}".format(results))

                            except Exception as err:
                                self.logger.error(err)
                                self.logger.error(traceback.format_exc())
                                from utilities.emailmanager import EmailManager
                                EmailManager.send_message('joseph.bonello@um.edu.mt', 'Prediction Error',
                                                          "\r\n".join(["In predictions_pfam.py", traceback.format_exc(),
                                                                       "", "Current sequence is", sequence]))
                                continue
                    else:
                        try:
                            term_file_funfams = term_file_dict[("dcgo", "FF")]

                            results = self.perform_prediction_dcgo(db_cursor, expected_go_terms,
                                                                   # str(funfams).replace("{", "(").replace("}", ")"),
                                                                   list(funfams),
                                                                   None, term_file_funfams,
                                                                   protein_name, protein_id)

                            if results is None:
                                return

                            # save to result file

                            pred_file_funfams.write("{},{},{},".format(protein_name, protein_id, None))
                            pred_file_funfams.write(",".join(str(v) for v in results))
                            pred_file_funfams.write("\n")

                            term_file_funfams.write("\n")

                            self.logger.info("Prediction Results: {}".format(results))

                        except Exception as err:
                            self.logger.error(err)
                            self.logger.error(traceback.format_exc())
                            from utilities.emailmanager import EmailManager
                            EmailManager.send_message('joseph.bonello@um.edu.mt', 'Prediction Error',
                                                      "\r\n".join(["In predictions_pfam.py", traceback.format_exc(),
                                                                   "", "Current sequence is", sequence]))
                            return

                for index in self.indexes_to_use:
                    if index == "dcgo":
                        continue  # Skip dcgo in this part

                    if self.separate_ontologies:
                        for ontology in ontologies_list:
                            try:
                                if ontology == 'BP':
                                    expected_go_terms = set(bp_terms)
                                elif ontology == 'MF':
                                    expected_go_terms = set(mf_terms)
                                else:
                                    expected_go_terms = set(cc_terms)

                                if len(expected_go_terms) == 0:  # There might not be terms in the ontology
                                    continue

                                term_file_superfams = term_file_dict[(index, "SF")]
                                term_file_funfams = term_file_dict[(index, "FF")]

                                results = self.perform_prediction_set_based(db_cursor, index,
                                                                            expected_go_terms,
                                                                            list(superfamilies), list(funfams),
                                                                            '_samefunfam_similarity', ontology,
                                                                            term_file_superfams, term_file_funfams,
                                                                            protein_name, protein_id)

                                if results is None:
                                    continue

                                # save to result file

                                results_superfam = results[0]
                                results_funfams = results[1]

                                pred_file_superfams = file_dict[(index, "SF")]
                                pred_file_superfams.write("{},{},{},".format(protein_name, protein_id, ontology))
                                pred_file_superfams.write(",".join(str(v) for v in results_superfam))
                                pred_file_superfams.write("\n")

                                pred_file_funfams = file_dict[(index, "FF")]
                                pred_file_funfams.write("{},{},{},".format(protein_name, protein_id, ontology))
                                pred_file_funfams.write(",".join(str(v) for v in results_funfams))
                                pred_file_funfams.write("\n")

                                term_file_superfams.write("\n")
                                term_file_funfams.write("\n")

                                self.logger.info("Prediction Results: {}".format(results))

                            except Exception as err:
                                self.logger.error(err)
                                self.logger.error(traceback.format_exc())
                                from utilities.emailmanager import EmailManager
                                EmailManager.send_message('joseph.bonello@um.edu.mt', 'Prediction Error',
                                                          "\r\n".join(["In predictions_pfam.py", traceback.format_exc(),
                                                                       "", "Current sequence is", sequence]))
                                continue
                    else:
                        try:
                            term_file_superfams = term_file_dict[(index, "SF")]
                            term_file_funfams = term_file_dict[(index, "FF")]

                            results = self.perform_prediction_set_based(db_cursor, index,
                                                                        expected_go_terms,
                                                                        list(superfamilies), list(funfams),
                                                                        '_samefunfam_similarity', None,
                                                                        term_file_superfams, term_file_funfams,
                                                                        protein_name, protein_id)

                            if results is None:
                                continue

                            # save to result file

                            results_superfam = results[0]
                            results_funfams = results[1]

                            pred_file_superfams = file_dict[(index, "SF")]
                            pred_file_superfams.write("{},{},{},".format(protein_name, protein_id, ontology))
                            pred_file_superfams.write(",".join(str(v) for v in results_superfam))
                            pred_file_superfams.write("\n")

                            pred_file_funfams = file_dict[(index, "FF")]
                            pred_file_funfams.write("{},{},{},".format(protein_name, protein_id, ontology))
                            pred_file_funfams.write(",".join(str(v) for v in results_funfams))
                            pred_file_funfams.write("\n")

                            term_file_superfams.write("\n")
                            term_file_funfams.write("\n")

                            self.logger.info("Prediction Results: {}".format(results))

                        except Exception as err:
                            self.logger.error(err)
                            self.logger.error(traceback.format_exc())
                            from utilities.emailmanager import EmailManager
                            EmailManager.send_message('joseph.bonello@um.edu.mt', 'Prediction Error',
                                                      "\r\n".join(["In predictions_pfam.py", traceback.format_exc(),
                                                                   "", "Current sequence is", sequence]))
                            continue

        except Exception as err:
            self.logger.error(err)
            self.logger.error(traceback.format_exc())
            from utilities.emailmanager import EmailManager
            EmailManager.send_message('joseph.bonello@um.edu.mt', 'Prediction Error',
                                      "\r\n".join(["In predictions_pfam.py", traceback.format_exc(),
                                                   "", "Current sequence is", sequence]))

        finally:
            if db_cursor is not None:
                db_cursor.close()

            if connection is not None:
                cp.close_connection(connection, pool_name="CACHE")

        for index in file_dict:
            text_file = file_dict[index]
            text_file.flush()

        for index in term_file_dict:
            term_file = term_file_dict[index]
            term_file.flush()

    def ontology_database_filters(self, ontology: str, separate_ontologies: bool, table_alias: str):
        ontology_column = ''
        ontology_filter = ''
        if separate_ontologies:
            if ontology is None:
                raise Exception("Ontology cannot be None if separate_ontologies option is selected.")
            else:
                ontology_column = f', {table_alias}.ontology'
                ontology_filter = f"and {table_alias}.ontology = '{ontology}'"
        return ontology_column, ontology_filter

    def predict_from_funfams_dcgo(self, ontology: str, separate_ontologies: bool, funfams: list,
                                  cursor_go_terms, enriched: bool = True):
        ontology_column, ontology_filter = self.ontology_database_filters(ontology, separate_ontologies, 't')

        # get_go_term_from_funfams = f"""
        #             SELECT dcgff.go_term AS goterm, dcgff.hscore, cffp.proportion as confidence {ontology_column}
        #             FROM {self.phdPrefix}.dcgoFunFams dcgff
        #             JOIN {self.phdPrefix}.PFAMFFProportions cffp ON dcgff.go_term = cffp.go_term AND dcgff.pfamfunfamilyfull_id = cffp.pfamfunfamilyfull_id
        #             JOIN {self.phdPrefix}.GoGraph_Terms t on cffp.go_term = t.go_term
        #             WHERE dcgff.pfamfunfamilyfull_id IN {str(funfams).replace('[', '(').replace(']', ')')}
        #             AND dcgff.hscore > 0 {ontology_filter}
        #             GROUP BY dcgff.go_term
        #         """

        get_go_term_from_funfams = f"""
            select dff.go_term AS goterm, sum(hscore) as sum {ontology_column}
      from {self.phdPrefix}.dcgoFunFams dff
      join {self.phdPrefix}.GoGraph_Terms t on dff.go_term = t.go_term
      where dff.pfamfunfamilyfull_id in {str(funfams).replace('[', '(').replace(']', ')')}
        AND dff.hscore > 0 {ontology_filter}
      group by dff.go_term;
        """

        get_go_terms_in_funfams = f"""
                    SELECT DISTINCT dcgff.go_term
                    FROM {self.phdPrefix}.dcgoFunFams dcgff
                    WHERE dcgff.pfamfunfamilyfull_id IN {str(funfams).replace('[', '(').replace(']', ')')}
                """

        cursor_go_terms.execute(get_go_term_from_funfams)
        scored_funfams_go_terms = cursor_go_terms.fetchall()

        cursor_go_terms.execute(get_go_terms_in_funfams)
        all_gterms = cursor_go_terms.fetchall()

        all_funfams_go_terms = set(s['go_term'] for s in all_gterms)

        hscores = {s['goterm']: s['sum'] for s in scored_funfams_go_terms}
        minhscore = min(hscores.values()) if len(hscores) > 0 else 0  # ensure there are predictions
        maxhscore = max(hscores.values()) if len(hscores) > 0 else 0
        denominator = maxhscore - minhscore
        pscores = {s['goterm']: ((s['sum'] - minhscore) / denominator) if denominator > 0 else 0 for s in
                   scored_funfams_go_terms}

        pscores_kept = {s: pscores[s] for s in pscores.keys() if pscores[s] > 0.5}  # dcGO limit

        scored_funfams_go_terms = list()
        for s in pscores_kept:
            r = {'goterm': s}
            r['confidence'] = pscores_kept[s]
            if separate_ontologies:
                r['ontology'] = ontology
            scored_funfams_go_terms.append(r)

        funfams_predicted_go_terms = set(s['goterm'] for s in scored_funfams_go_terms)

        final_predicted_terms = funfams_predicted_go_terms

        parental_go_terms = None

        if len(funfams_predicted_go_terms) > 0:
            if enriched:
                parental_go_terms = self.go_tool.get_parental_terms_for_list(list(funfams_predicted_go_terms))

                pgt = set()
                for pTerm in parental_go_terms:
                    if pTerm not in ['all', 'GO:0005575', 'GO:0008150', 'GO:0003674']:
                        pgt.add(pTerm)

                final_predicted_terms = funfams_predicted_go_terms.union(pgt)
            else:
                final_predicted_terms = funfams_predicted_go_terms

        return scored_funfams_go_terms, final_predicted_terms, parental_go_terms, all_funfams_go_terms

    def perform_prediction_dcgo(self, cursor_go_terms, expected_go_terms, funfams, ontology,
                                term_file_funfams, protein_name, protein_id):
        result = list()

        # Preliminary part: Get Expected Terms

        for expTerm in expected_go_terms:
            outstr = f'{protein_name},{protein_id},{ontology if ontology is not None else ""},'
            outstr += f"{expTerm},1.0,Y,N\n"
            term_file_funfams.write(outstr)

        if self.inherit_go_terms:
            gold_standard_parental_go_terms = self.go_tool.get_parental_terms_for_list(list(expected_go_terms))

            gspgt = set()
            for pTerm in gold_standard_parental_go_terms:
                if pTerm not in ['all', 'GO:0005575', 'GO:0008150', 'GO:0003674']:  # Remove GO root terms
                    gspgt.add(pTerm)
                    outstr = f'{protein_name},{protein_id},{ontology if ontology is not None else ""},'
                    outstr += f"{pTerm},1.0,Y,Y\n"
                    term_file_funfams.write(outstr)

            final_expected_terms = set(expected_go_terms.union(gspgt))
        else:
            final_expected_terms = expected_go_terms

        result.append(len(final_expected_terms))

        # Get Funfam terms and scores

        predictions = self.predict_from_funfams_dcgo(ontology, self.separate_ontologies, funfams,
                                                     cursor_go_terms, self.inherit_go_terms)

        scored_funfams_go_terms = predictions[0]
        final_predicted_terms = predictions[1]
        parental_go_terms = predictions[2]
        all_funfams_go_terms = predictions[3]

        # output.append("Funfams GO Terms: {}".format(set(s[0] for s in scored_funfams_go_terms)))

        for pTerm in scored_funfams_go_terms:
            outstr = f'{protein_name},{protein_id},{ontology if ontology is not None else ""},'
            outstr += f"{pTerm['goterm']},{str(pTerm['confidence'])},N,N\n"
            term_file_funfams.write(outstr)

        if parental_go_terms is not None:
            for pTerm in parental_go_terms:
                outstr = f'{protein_name},{protein_id},{ontology if ontology is not None else ""},'

                cursor_go_terms.execute(self.sql_term_probability, pTerm)

                score = cursor_go_terms.fetchone()
                score = score['probability'] if score is not None else 0

                outstr += f"{pTerm},{score},N,Y\n"
                term_file_funfams.write(outstr)

        result.append(len(final_predicted_terms) if final_predicted_terms is not None else 0)

        if final_predicted_terms is not None and len(final_predicted_terms) > 0:

            intersection_expected_predicted = set(final_expected_terms).intersection(final_predicted_terms)
            funfams_tp = len(intersection_expected_predicted)  # true positive
            funfams_fp = len(set(final_predicted_terms).difference(final_expected_terms))  # false positive
            # superfamily_fn = len(set(expected_go_terms).difference(superfamily_predicted_go_terms))  # false negative
            funfams_fn = abs(len(set(final_expected_terms)) - funfams_tp)
            funfams_tn = len(all_funfams_go_terms.difference(final_expected_terms.union(final_predicted_terms)))

            # CAFA Style
            funfams_precision = len(intersection_expected_predicted) / len(final_predicted_terms)
            funfams_recall = len(intersection_expected_predicted) / len(final_expected_terms)

            # TODO in results, check when fp + tn == 0
            funfams_specificity = 0 if (funfams_fp + funfams_tn) == 0 else funfams_tn / (
                    funfams_fp + funfams_tn)  # True Negative Rate
            if funfams_tp == 0:
                funfams_f1_score = 0
            else:
                funfams_f1_score = ((2 * funfams_precision * funfams_recall) / (funfams_precision + funfams_recall))
            funfams_false_negative_rate = 0 if (funfams_tp + funfams_fn) == 0 else funfams_fn / (
                    funfams_fn + funfams_tp)  # Miss Rate
            funfams_false_positive_rate = 0 if (funfams_fp + funfams_tn) == 0 else funfams_fp / (
                    funfams_fp + funfams_tn)
            funfams_false_discovery_rate = 0 if (funfams_fp + funfams_tp) == 0 else funfams_fp / (
                    funfams_fp + funfams_tp)
            funfams_fowlkes_mallows_index = sqrt(funfams_precision * funfams_recall)
            funfams_accuracy = (funfams_tp + funfams_tn) / (funfams_tp + funfams_tn + funfams_fp + funfams_fn)

            result.append(funfams_tp)
            result.append(funfams_fp)
            result.append(funfams_fn)
            result.append(funfams_tn)
            result.append(funfams_precision)
            result.append(funfams_recall)
            result.append(funfams_specificity)
            result.append(funfams_f1_score)
            result.append(funfams_false_negative_rate)
            result.append(funfams_false_positive_rate)
            result.append(funfams_false_discovery_rate)
            result.append(funfams_fowlkes_mallows_index)
            result.append(funfams_accuracy)

        return result

    def predict_from_superfamilies_set_based(self, ontology: str, separate_ontologies: bool, superfamilies: list,
                                             index: str, cursor_go_terms, enriched: bool = True):
        ontology_column, ontology_filter = self.ontology_database_filters(ontology, separate_ontologies, 't')

        get_go_terms_from_superfamilies = f"""
                        select cfs.goterm{ontology_column},
                          cfs.{index} as score
                        from {self.phdPrefix}.pfamfamiliessimilarity cfs
                        join {self.phdPrefix}.GoGraph_Terms t on cfs.goterm = t.go_term
                        where cfs.pfamfam in {str(superfamilies).replace('[', '(').replace(']', ')')}
                        and cfs.{index} >= ({self.threshold} * 0.01) {ontology_filter}
                        order by cfs.{index};
                    """  # Arbitrarily set due to high dispersion in superfamilies

        get_go_terms_in_superfamily = f"""
                        SELECT distinct cfs1.goterm
                        FROM  {self.phdPrefix}.pfamfamiliessimilarity cfs1
                        WHERE cfs1.pfamfam in {str(superfamilies).replace('[', '(').replace(']', ')')};
                    """

        final_predicted_terms = None

        cursor_go_terms.execute(get_go_terms_from_superfamilies)
        scored_superfamily_go_terms = cursor_go_terms.fetchall()

        cursor_go_terms.execute(get_go_terms_in_superfamily)
        all_superfamily_go_terms = set(s['goterm'] for s in cursor_go_terms.fetchall())

        superfamily_predicted_go_terms = set(s['goterm'] for s in scored_superfamily_go_terms)

        parental_go_terms = None
        final_predicted_terms = superfamily_predicted_go_terms

        if len(superfamily_predicted_go_terms) > 0:
            if enriched:
                parental_go_terms = self.go_tool.get_parental_terms_for_list(list(superfamily_predicted_go_terms))

                pgt = set()
                for pTerm in parental_go_terms:
                    if pTerm not in ['all', 'GO:0005575', 'GO:0008150', 'GO:0003674']:  # Remove GO root terms
                        pgt.add(pTerm)

                final_predicted_terms = superfamily_predicted_go_terms.union(pgt)
            else:
                final_predicted_terms = superfamily_predicted_go_terms

        return scored_superfamily_go_terms, final_predicted_terms, parental_go_terms, all_superfamily_go_terms

    def predict_from_funfams_set_based(self, tableExtension: str, ontology: str, separate_ontologies: bool,
                                       funfams: list, index: str, cursor_go_terms, enriched: bool = True):
        ontology_column, ontology_filter = self.ontology_database_filters(ontology, separate_ontologies, 't')

        get_go_term_from_funfams = f"""
                        select DISTINCT
                            cffs.goterm{ontology_column},
                            cffs.{index} as score
                        from
                            {self.phdPrefix}.pfamfunfamiliessimilarity{tableExtension} cffs
                        join {self.phdPrefix}.GoGraph_Terms t on cffs.goterm = t.go_term
                        where cffs.pfamfunfam in {str(funfams).replace('[', '(').replace(']', ')')}
                        {ontology_filter}
                        and cffs.{index} >= {self.threshold}
                        order by cffs.{index}, cffs.goterm;
                    """  # default threshold 0.09

        get_go_terms_in_funfams = f"""
                        SELECT distinct cffs1.goterm
                        FROM  {self.phdPrefix}.pfamfunfamiliessimilarity{tableExtension} cffs1
                        WHERE cffs1.pfamfunfam in {str(funfams).replace('[', '(').replace(']', ')')}
                    """

        cursor_go_terms.execute(get_go_term_from_funfams)
        scored_funfams_go_terms = cursor_go_terms.fetchall()

        cursor_go_terms.execute(get_go_terms_in_funfams)
        all_funfams_go_terms = set(s['goterm'] for s in cursor_go_terms.fetchall())

        funfams_predicted_go_terms = set(s['goterm'] for s in scored_funfams_go_terms)

        final_predicted_terms = funfams_predicted_go_terms
        parental_go_terms = None

        if len(funfams_predicted_go_terms) > 0:
            if enriched:
                parental_go_terms = self.go_tool.get_parental_terms_for_list(list(funfams_predicted_go_terms))

                pgt = set()
                for pTerm in parental_go_terms:
                    if pTerm not in ['all', 'GO:0005575', 'GO:0008150', 'GO:0003674']:
                        pgt.add(pTerm)

                final_predicted_terms = funfams_predicted_go_terms.union(pgt)
            else:
                final_predicted_terms = funfams_predicted_go_terms

        return scored_funfams_go_terms, final_predicted_terms, parental_go_terms, all_funfams_go_terms

    def perform_prediction_set_based(self, cursor_go_terms, index, expected_go_terms, superfamilies,
                                     funfams, tableExtension, ontology, term_file_superfams, term_file_funfams,
                                     protein_name, protein_id):
        # Results for superfamilies (0) and funfams (1)
        results = (list(), list())

        ontology_column, ontology_filter = self.ontology_database_filters(ontology, self.separate_ontologies, 't')

        # Preliminary part: Get Expected Terms

        final_expected_terms = None

        for expTerm in expected_go_terms:
            outstr = f'{protein_name},{protein_id},{ontology if ontology is not None else ""},'
            outstr += f"{expTerm},1.0,Y,N\n"

            term_file_superfams.write(outstr)
            term_file_funfams.write(outstr)

        if self.inherit_go_terms:
            gold_standard_parental_go_terms = self.go_tool.get_parental_terms_for_list(list(expected_go_terms))

            gspgt = set()
            for pTerm in gold_standard_parental_go_terms:
                if pTerm not in ['all', 'GO:0005575', 'GO:0008150', 'GO:0003674']:  # Remove GO root terms
                    gspgt.add(pTerm)
                    outstr = f'{protein_name},{protein_id},{ontology if ontology is not None else ""},'
                    outstr += f"{pTerm},1.0,Y,Y\n"

                    term_file_superfams.write(outstr)
                    term_file_funfams.write(outstr)

            final_expected_terms = set(expected_go_terms.union(gspgt))
        else:
            final_expected_terms = expected_go_terms

        results[0].append(len(final_expected_terms))
        results[1].append(len(final_expected_terms))

        # Part 1; Get Go Terms from SuperFamilies
        result = results[0]

        predictions = self.predict_from_superfamilies_set_based(ontology, self.separate_ontologies, superfamilies,
                                                                index,
                                                                cursor_go_terms, self.inherit_go_terms)
        scored_superfamily_go_terms = predictions[0]
        final_predicted_terms = predictions[1]
        parental_go_terms = predictions[2]
        all_superfamily_go_terms = predictions[3]

        result.append(len(final_predicted_terms))

        # output.append("Superfamily GO Terms: {}".format(set(s[0] for s in scored_superfamily_go_terms)))

        for pTerm in scored_superfamily_go_terms:
            outstr = f'{protein_name},{protein_id},{ontology if ontology is not None else ""},'
            outstr += f"{pTerm['goterm']},{str(pTerm['score'])},N,N\n"
            term_file_superfams.write(outstr)

        if parental_go_terms is not None:
            for pTerm in parental_go_terms:
                if pTerm not in ['all', 'GO:0005575', 'GO:0008150', 'GO:0003674']:  # Remove GO root terms
                    cursor_go_terms.execute(self.sql_term_probability, pTerm)

                    score = cursor_go_terms.fetchone()
                    score = score['probability'] if score is not None else 0

                    outstr = f'{protein_name},{protein_id},{ontology if ontology is not None else ""},'
                    outstr += f"{pTerm},{score},N,Y\n"
                    term_file_superfams.write(outstr)

        if final_predicted_terms is not None and len(final_predicted_terms) > 0:

            intersection_expected_predicted = set(final_expected_terms).intersection(final_predicted_terms)
            superfamily_tp = len(intersection_expected_predicted)  # true positive
            superfamily_fp = len(set(final_predicted_terms).difference(final_expected_terms))  # false positive
            # superfamily_fn = len(set(expected_go_terms).difference(superfamily_predicted_go_terms))  # false negative
            superfamily_fn = abs(len(set(final_expected_terms)) - superfamily_tp)
            superfamily_tn = len(all_superfamily_go_terms.difference(final_expected_terms.union(final_predicted_terms)))

            # CAFA Style
            superfamily_precision = len(intersection_expected_predicted) / len(final_predicted_terms)
            superfamily_recall = len(intersection_expected_predicted) / len(final_expected_terms)

            superfamily_specificity = 0 if (superfamily_fp + superfamily_tn) == 0 else superfamily_tn / (
                    superfamily_fp + superfamily_tn)  # True Negative Rate
            if superfamily_precision == 0:
                superfamily_f1_score = 0
            else:
                superfamily_f1_score = (
                        (2 * superfamily_precision * superfamily_recall) / (superfamily_precision + superfamily_recall))
            superfamily_false_negative_rate = 0 if (superfamily_tp + superfamily_fn) == 0 else superfamily_fn / (
                    superfamily_fn + superfamily_tp)  # Miss Rate
            superfamily_false_positive_rate = 0 if (superfamily_fp + superfamily_tn) == 0 else superfamily_fp / (
                    superfamily_fp + superfamily_tn)
            superfamily_false_discovery_rate = 0 if (superfamily_fp + superfamily_tp) == 0 else superfamily_fp / (
                    superfamily_fp + superfamily_tp)
            superfamily_fowlkes_mallows_index = sqrt(superfamily_precision * superfamily_recall)
            superfamily_accuracy = (superfamily_tp + superfamily_tn) / (
                    superfamily_tp + superfamily_tn + superfamily_fp + superfamily_fn)

            result.append(superfamily_tp)
            result.append(superfamily_fp)
            result.append(superfamily_fn)
            result.append(superfamily_tn)
            result.append(superfamily_precision)
            result.append(superfamily_recall)
            result.append(superfamily_specificity)
            result.append(superfamily_f1_score)
            result.append(superfamily_false_negative_rate)
            result.append(superfamily_false_positive_rate)
            result.append(superfamily_false_discovery_rate)
            result.append(superfamily_fowlkes_mallows_index)
            result.append(superfamily_accuracy)
        else:
            return None

        # Part 2: Get Go Terms from Funfams

        result = results[1]
        intersection_expected_predicted = None

        predictions = self.predict_from_funfams_set_based(tableExtension, ontology, self.separate_ontologies, funfams,
                                                          index, self.inherit_go_terms)

        scored_funfams_go_terms = predictions[0]
        final_predicted_terms = predictions[1]
        parental_go_terms = predictions[2]
        all_funfams_go_terms = predictions[3]

        # output.append("Funfams GO Terms: {}".format(set(s[0] for s in scored_funfams_go_terms)))

        for pTerm in scored_funfams_go_terms:
            outstr = f'{protein_name},{protein_id},{ontology if ontology is not None else ""},'
            outstr += f"{pTerm['goterm']},{pTerm['score']},N,N\n"
            term_file_funfams.write(outstr)

        if parental_go_terms is not None:
            for pTerm in parental_go_terms:
                if pTerm not in ['all', 'GO:0005575', 'GO:0008150', 'GO:0003674']:  # Remove GO root terms
                    outstr = f'{protein_name},{protein_id},{ontology if ontology is not None else ""},'

                    cursor_go_terms.execute(self.sql_term_probability, pTerm)

                    score = cursor_go_terms.fetchone()
                    score = score['probability'] if score is not None else 0

                    outstr += f"{pTerm},{score},N,Y\n"

                    term_file_funfams.write(outstr)

        result.append(len(final_predicted_terms) if final_predicted_terms is not None else 0)

        if final_predicted_terms is not None and len(final_predicted_terms) > 0:
            intersection_expected_predicted = set(final_expected_terms).intersection(final_predicted_terms)
            funfams_tp = len(intersection_expected_predicted)  # true positive
            funfams_fp = len(set(final_predicted_terms).difference(final_expected_terms))  # false positive
            # superfamily_fn = len(set(expected_go_terms).difference(superfamily_predicted_go_terms))  # false negative
            funfams_fn = abs(len(set(final_expected_terms)) - funfams_tp)
            funfams_tn = len(all_funfams_go_terms.difference(final_expected_terms.union(final_predicted_terms)))

            # CAFA Style
            funfams_precision = len(intersection_expected_predicted) / len(final_predicted_terms)
            funfams_recall = len(intersection_expected_predicted) / len(final_expected_terms)

            # TODO in results, check when fp + tn == 0
            funfams_specificity = 0 if (funfams_fp + funfams_tn) == 0 else funfams_tn / (
                    funfams_fp + funfams_tn)  # True Negative Rate
            if funfams_tp == 0:
                funfams_f1_score = 0
            else:
                funfams_f1_score = ((2 * funfams_precision * funfams_recall) / (funfams_precision + funfams_recall))
            funfams_false_negative_rate = 0 if (funfams_tp + funfams_fn) == 0 else funfams_fn / (
                    funfams_fn + funfams_tp)  # Miss Rate
            funfams_false_positive_rate = 0 if (funfams_fp + funfams_tn) == 0 else funfams_fp / (
                    funfams_fp + funfams_tn)
            funfams_false_discovery_rate = 0 if (funfams_fp + funfams_tp) == 0 else funfams_fp / (
                    funfams_fp + funfams_tp)
            funfams_fowlkes_mallows_index = sqrt(funfams_precision * funfams_recall)
            funfams_accuracy = (funfams_tp + funfams_tn) / (funfams_tp + funfams_tn + funfams_fp + funfams_fn)

            result.append(funfams_tp)
            result.append(funfams_fp)
            result.append(funfams_fn)
            result.append(funfams_tn)
            result.append(funfams_precision)
            result.append(funfams_recall)
            result.append(funfams_specificity)
            result.append(funfams_f1_score)
            result.append(funfams_false_negative_rate)
            result.append(funfams_false_positive_rate)
            result.append(funfams_false_discovery_rate)
            result.append(funfams_fowlkes_mallows_index)
            result.append(funfams_accuracy)
        else:
            return None

        return results

    def predict_from_path(self, path, output, parallel_backend='threading'):
        import os
        for file in os.listdir(path):
            print(file)

            sequences = set()
            sequence = ''

            with open(os.path.join(path, file), 'r') as file_handle:
                for line in file_handle:
                    if line.strip() == '' or line.find('>') > -1:  # Handle new sequences; either through a line break
                        # or through the discovery of a new sequence right after this one ends
                        if sequence.strip() != '':
                            sequences.add(sequence)  # Checks if empty
                            sequence = ''
                        if line.find('>') > -1:
                            if sequence != '':
                                sequences.add(sequence)
                            sequence = line
                    else:
                        sequence += line.strip()

            if sequence != '':
                sequences.add(sequence)  # set will handle doubles; avoids mossing last sequence in the file

            # self.logger.debugsequences)

            file_dict = dict()
            term_file_dict = dict()

            cafa_set = file[0:file.find('.')]

            # cache_conn = sqlite3.connect(self.cache_properties['NAME'])
            from utilities.mysqlconnectionpool import MySQLConnectionPool
            connection = MySQLConnectionPool.get_instance().get_instance(pool_name="CACHE")

            try:

                for index in self.indexes_to_use:
                    import os

                    if not os.path.exists(os.path.join(tempfile.gettempdir(), "predictions", self.ts)):
                        os.makedirs(os.path.join(tempfile.gettempdir(), "predictions", self.ts))

                    if (index != 'dcgo'):
                        pred_file_superfams = open(
                            os.path.join(tempfile.gettempdir(), "predictions", self.ts,
                                         "prediction-{0}-{1}-t{2}-SF.csv".format(index, cafa_set,
                                                                                 str(self.threshold).replace('.',
                                                                                                             '_'))),
                            "w")
                    pred_file_funfams = open(
                        os.path.join(tempfile.gettempdir(), "predictions", self.ts,
                                     "prediction-{0}-{1}-t{2}-FF.csv".format(index, cafa_set,
                                                                             str(self.threshold).replace('.', '_'))),
                        "w")

                    if index != 'dcgo':
                        term_file_superfams = open(
                            os.path.join(tempfile.gettempdir(), "predictions", self.ts,
                                         "terms-{0}-{1}-t{2}-SF.csv".format(index, cafa_set,
                                                                            str(self.threshold).replace('.', '_'))),
                            "w")
                    term_file_funfams = open(
                        os.path.join(tempfile.gettempdir(), "predictions", self.ts,
                                     "terms-{0}-{1}-t{2}-FF.csv".format(index, cafa_set,
                                                                        str(self.threshold).replace('.', '_'))), "w")

                    if index != 'dcgo':
                        file_dict[(index, 'SF')] = pred_file_superfams
                    file_dict[(index, 'FF')] = pred_file_funfams

                    if index != 'dcgo':
                        term_file_dict[(index, 'SF')] = term_file_superfams
                    term_file_dict[(index, 'FF')] = term_file_funfams

                    if index != 'dcgo':
                        pred_file_superfams.write(
                            "Protein Name,Protein Unique Id,{}expected_count,superfamily_predicted_count,superfamily_tp,"
                            "superfamily_fp,superfamily_fn,superfamily_tn,"
                            "superfamily_precision,superfamily_recall,superfamily_specificity,superfamily_f1_score,"
                            "superfamily_false_negative_rate,superfamily_false_positive_rate,superfamily_false_discovery_rate,"
                            "superfamily_fowlkes_mallows_index,superfamily_accuracy"
                                .format("Ontology, " if self.separate_ontologies else ""))

                    pred_file_funfams.write(
                        "Protein Name,Protein Unique Id,{}expected_count,funfams_predicted_count,funfams_tp,"
                        "funfams_fp,funfams_fn,funfams_tn,funfams_precision,funfams_recall,"
                        "funfams_specificity,funfams_f1_score,funfams_false_negative_rate,funfams_false_positive_rate,"
                        "funfams_false_discovery_rate, funfams_fowlkes_mallows_index,funfams_accuracy"
                            .format("Ontology, " if self.separate_ontologies else ""))

                    if index != 'dcgo':
                        term_file_superfams.write("Protein Name,Protein Unique Id,{}GO Term,Score,Is From TPS,"
                                                  "Is Inherited From Parents"
                                                  .format("Ontology, " if self.separate_ontologies else ""))

                    term_file_funfams.write("Protein Name,Protein Unique Id,{}GO Term,Score,Is From TPS,"
                                            "Is Inherited From Parents"
                                            .format("Ontology, " if self.separate_ontologies else ""))

                    # if self.inherit_go_terms:
                    #    text_file.write(", Inherited GO Terms, Count Correct Inherited")

                    pred_file_superfams.write("\n")
                    pred_file_superfams.flush()

                    pred_file_funfams.write("\n")
                    pred_file_funfams.flush()

                    term_file_superfams.write("\n")
                    term_file_superfams.flush()

                    term_file_funfams.write("\n")
                    term_file_funfams.flush()

                with connection.cursor() as cache_cursor:

                    # Create if the cache table does not exist
                    sql_check_table = f"""
                        CREATE TABLE IF NOT EXISTS {self.phdCachePrefix}PFAM.pred_!
                        (
                            id INT PRIMARY KEY NOT NULL AUTO_INCREMENT,
                            hashed_sequence VARCHAR(128),
                            protein_description TEXT,
                            protein_id TEXT,
                            protein_name TEXT,
                            related_sfamilies LONGBLOB,
                            related_ffamilies LONGBLOB,
                            expected_go_terms LONGBLOB,
                            expected_go_terms_bp LONGBLOB,
                            expected_go_terms_cc LONGBLOB,
                            expected_go_terms_mf LONGBLOB
                        );
                    """.replace("!", cafa_set)

                    sql_rollback_expected_results = f"""
                        ALTER TABLE {self.phdCachePrefix}PFAM.pred_!
                        ADD COLUMN IF NOT EXISTS expected_go_terms_rollback LONGBLOB,
                        ADD COLUMN IF NOT EXISTS expected_go_terms_rollback_bp LONGBLOB,
                        ADD COLUMN IF NOT EXISTS expected_go_terms_rollback_cc LONGBLOB,
                        ADD COLUMN IF NOT EXISTS expected_go_terms_rollback_mf LONGBLOB;
                    """.replace("!", cafa_set)

                    sql_constraint = f"""
                        CREATE UNIQUE INDEX IF NOT EXISTS pred_{0}_hashed_sequence_uindex ON {self.phdCachePrefix}PFAM.pred_{0} (hashed_sequence);
                    """.replace("!", cafa_set)

                    cache_cursor.execute(sql_check_table)
                    cache_cursor.execute(sql_rollback_expected_results)
                    cache_cursor.execute(sql_constraint)

                import os
                max_threads = int(os.getenv('NUM_THREADS'))

                sequences = list(sequences)

                Parallel(n_jobs=max_threads, verbose=5, backend=parallel_backend)(delayed(unwrap_self_predict)
                                                                                  (sequence, output=output,
                                                                                   cafa_set=cafa_set,
                                                                                   file_dict=file_dict,
                                                                                   term_file_dict=term_file_dict)
                                                                                  for sequence in
                                                                                  zip([self] * len(sequences),
                                                                                      sequences))
            except Exception as err:
                self.logger.error(err)
                self.logger.error(traceback.format_exc())
                from utilities.emailmanager import EmailManager
                EmailManager.send_message('joseph.bonello@um.edu.mt', 'Prediction Error',
                                          "\r\n".join(["In predictions_pfam.py", traceback.format_exc()]))

            finally:
                self.logger.warn("In finally ...")

                if connection is not None:
                    MySQLConnectionPool.close_connection(connection, pool_name="CACHE")

                for index in file_dict:
                    text_file = file_dict[index]
                    text_file.flush()
                    text_file.close()

                for index in term_file_dict:
                    term_file = term_file_dict[index]
                    term_file.flush()
                    term_file.close()
