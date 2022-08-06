import io
import traceback
from hashlib import blake2b
from math import sqrt
from typing import List, Set, TextIO, Tuple

import pymysql
from Bio import SeqIO
from Bio.Seq import Seq
from joblib import Parallel, delayed

from GoGraph.classes.blast import Blast
from GoGraph.classes.funfhmmer import Funfhmmer
from GoGraph.classes.geneontology import GeneOntologyUtility
from GoGraph.classes.searches import Searches


def unwrap_self_predict(arg, **kwarg):
    return Predictions_CATH.predict(*arg, **kwarg)


# noinspection SqlResolve
class Predictions_CATH:
    __slots__ = ['threshold', 'inherit_go_terms', 'separate_ontologies', 'indexes_to_use', 'ts', 'evidenceCodes', 'go_tool', 'logger', 'pool', 'sql_term_probability', 'funfhmmer',
                 'use_rollback', 'sql_update_funfams', 'b2b', 'sql_update_hashed_sequence', 'sql_update_expected_terms', 'apply_cache_update_to_db', 'phdPrefix', 'phdCachePrefix']

    def __init__(self, threshold: float, inherit_go_terms: bool = False, separate_ontologies: bool = False, indexes_to_use: list = [], use_rollback: bool = False,
                 apply_cache_update_to_db: bool = False, log_to_file: bool = True):

        # evidenceCodes = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'TAS', 'IEP', 'IC']
        self.evidenceCodes = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'TAS', 'IEP', 'IC', 'ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'IBA', 'IBD', 'IKR', 'IRD', 'RCA']  # , 'IEA']

        self.go_tool = GeneOntologyUtility()

        self.inherit_go_terms = inherit_go_terms
        self.separate_ontologies = separate_ontologies
        self.indexes_to_use = indexes_to_use
        self.threshold = threshold

        self.funfhmmer = Funfhmmer()

        import datetime
        self.ts = datetime.datetime.now().strftime("%Y%m%d%H%M%S")

        import socket
        import os
        tmp_path = os.getenv('PREDICTION_SAVE_PATH_' + socket.gethostname().upper())
        self.phdPrefix = os.getenv("PHD_DB_PREFIX")
        self.phdCachePrefix = os.getenv("PHDCACHE_DB_PREFIX")

        from GoGraph.classes.decorators import create_logger
        import logging
        if log_to_file:
            self.logger = create_logger(os.path.join(tmp_path, 'predictions', 'phd.log'), loglevel=logging.DEBUG)
        else:
            self.logger = create_logger(loglevel=logging.DEBUG)

        self.sql_term_probability = f"""
                SELECT format(ic.probability,2) as "probability"
                FROM {self.phdPrefix}.infocontentcafa ic
                WHERE ic.goterm = %s;
            """

        # self.sql_term_probability = f"""
        #     select round(sin( i.infocontent / a.maxinfocontent),2) as "probability"
        #     from {self.phdPrefix}.infocontentcafa i,
        #     (
        #     select max(ic.infocontent) as maxinfocontent
        #     from {self.phdPrefix}.infocontentcafa ic
        #     ) a
        #     where i.goterm = %s;
        # """

        self.use_rollback = use_rollback

        from utilities.mysqlconnectionpool import MySQLConnectionPool
        self.pool = MySQLConnectionPool.get_instance()

        self.sql_update_funfams = f"""
            UPDATE {self.phdCachePrefix}.pred_! SET related_sfamilies = %s, related_ffamilies = %s WHERE protein_description = %s;
        """

        self.sql_update_hashed_sequence = f"""
            UPDATE {self.phdCachePrefix}.pred_! SET hashed_sequence = %s WHERE protein_description = %s
        """

        self.sql_update_expected_terms = f"""
                        UPDATE {self.phdCachePrefix}.pred_! SET  expected_go_terms = %s, expected_go_terms_bp = %s, expected_go_terms_cc = %s, expected_go_terms_mf = %s WHERE 
                        protein_description = %s;
                    """

        self.b2b = blake2b(digest_size=blake2b.MAX_DIGEST_SIZE, key=b"JBPhD", salt=b"PhD is SALTY")

        self.apply_cache_update_to_db = apply_cache_update_to_db

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.logger.debug("Closing MySQL Connections")

    def get_data_from_cache(self, cafa_set: str, hashed_sequence: str, sequence: str, funfhmmer: Funfhmmer, protein_description: str, searchable_sequence: Seq,
                            output: List, search: Searches, skip_sanity_checks: bool = False):
        sql_retrieve_protein = f"""
            SELECT * FROM {self.phdCachePrefix}.pred_! WHERE protein_description = %s
        """

        protein_id = None
        protein_name = None
        expected_go_terms = None
        bp_terms = None
        cc_terms = None
        mf_terms = None
        superfamilies = None
        funfams = None

        cache_cursor = None
        cache_connection = self.pool.get_connection(pool_name='CACHE')

        try:
            # with sqlite3.connect(cache_props['NAME'], check_same_thread=True, timeout=60.0) as cache_conn:

            with cache_connection.cursor() as cache_cursor:
                # self.logger.debug(f"SQL: {sql_retrieve_protein.replace("!", cafa_set)}; Protein Description: {protein_description}")
                cache_cursor.execute(sql_retrieve_protein.replace("!", cafa_set), (protein_description,))

                row = cache_cursor.fetchone()

                # Step 1: Get the necessary data from Cache or Search
                self.logger.debug("Step 1 - Processing Set {}".format(cafa_set))
                if row is not None:
                    if row['hashed_sequence'] != hashed_sequence and self.apply_cache_update_to_db:
                        cache_cursor.execute(self.sql_update_hashed_sequence.replace("!", cafa_set), (hashed_sequence, protein_description,))

                    self.logger.debug("Step 1[DB]")
                    protein_id = row["protein_id"]
                    protein_name = row["protein_name"]

                    try:
                        superfamilies = set(row['related_sfamilies'].split(','))
                    except:
                        superfamilies = None

                    try:
                        funfams = set(row['related_ffamilies'].split(','))
                    except:
                        funfams = None

                    if self.use_rollback:
                        expected_go_terms = set(row['expected_go_terms_rollback'].split(','))
                        bp_terms = set(row['expected_go_terms_rollback_bp'].split(','))
                        cc_terms = set(row['expected_go_terms_rollback_cc'].split(','))
                        mf_terms = set(row['expected_go_terms_rollback_mf'].split(','))
                    else:
                        expected_go_terms = set(str(row['expected_go_terms']).split(',')) if str(row['expected_go_terms']) is not None else None
                        bp_terms = set(str(row['expected_go_terms_bp']).split(',')) if str(row['expected_go_terms_bp']) is not None else None
                        cc_terms = set(str(row['expected_go_terms_cc']).split(',')) if str(row['expected_go_terms_cc']) is not None else None
                        mf_terms = set(str(row['expected_go_terms_mf']).split(',')) if str(row['expected_go_terms_mf']) is not None else None

                    self.logger.debug(f"Processing {protein_name} (Blasted protein id {protein_id}) - Target {protein_description} from cache")
                    self.logger.debug(f"Superfamilies: <{superfamilies}>")
                    self.logger.debug(f"Funfams: <{funfams}>")
                    self.logger.debug(f"Expected Go Terms: <{expected_go_terms}>")

            # Sanity Checks
            if not skip_sanity_checks:
                sanity_check_results = self.sanity_checks(cafa_set, protein_id, protein_name, protein_description, sequence, searchable_sequence, funfhmmer,
                                                          superfamilies, funfams, expected_go_terms, bp_terms, cc_terms, mf_terms, search, output)
                if not sanity_check_results[0]:
                    return protein_id, protein_name, None, None, None, None, None, None
                else:
                    superfamilies = sanity_check_results[1]
                    funfams = sanity_check_results[2]
                    expected_go_terms = sanity_check_results[3]
                    bp_terms = sanity_check_results[4]
                    cc_terms = sanity_check_results[5]
                    mf_terms = sanity_check_results[6]
        except Exception as err:
            self.logger.error("An Exception in CACHE lookup occurred")
            self.logger.error(err)
            self.logger.error(traceback.format_exc())
            from utilities.emailmanager import EmailManager
            EmailManager.send_message('joseph.bonello@um.edu.mt', 'Prediction Error',
                                      "\r\n".join(["In predictions_cath.py", traceback.format_exc(), "", "Current sequence is", sequence]))
            return protein_id, protein_name, None, None, None, None, None, None

        finally:
            if cache_cursor is not None:
                cache_cursor.close()

            self.pool.close_connection(cache_connection, pool_name='CACHE')

        if protein_id is None or protein_name is None:
            if protein_description is not None:
                protein_id = protein_description.split(' ')[0]
                protein_name = protein_description.split(' ')[1]

        return protein_id, protein_name, expected_go_terms, bp_terms, cc_terms, mf_terms, superfamilies, funfams

    def calculate_protein_prediction_metrics(self, final_expected_terms, final_predicted_terms, all_family_go_terms=None):
        # Some references on re-implementation without TN
        # https://stats.stackexchange.com/questions/61829/given-true-positive-false-negative-rates-can-you-calculate-false-positive-tru

        result = list()

        intersection_expected_predicted = set(final_expected_terms).intersection(set(final_predicted_terms))
        family_tp = len(intersection_expected_predicted)  # true positive
        family_fp = len(set(final_predicted_terms).difference(final_expected_terms))  # false positive
        # family_fn = len(set(expected_go_terms).difference(family_predicted_go_terms))  # false negative
        family_fn = abs(len(set(final_expected_terms)) - family_tp)
        family_tn = len(all_family_go_terms.difference(final_expected_terms.union(final_predicted_terms)))

        # CAFA Style
        family_precision = 0 if len(final_predicted_terms) == 0 else len(intersection_expected_predicted) / len(final_predicted_terms)
        family_recall = 0 if len(final_expected_terms) == 0 else len(intersection_expected_predicted) / len(final_expected_terms)

        # Alternative
        family_precision_alternative = 0 if (family_tp + family_fp) == 0 else family_tp / (family_tp + family_fp)
        family_recall_alternative = 0 if (family_tp + family_fn) == 0 else family_tp / (family_tp + family_fn)

        # family_specificity = 0 if (family_fp + family_tn) == 0 else family_tn / (family_fp + family_tn)  # True Negative Rate

        if family_precision == 0 and family_recall == 0:
            family_f1_score = 0
        else:
            family_f1_score = ((2 * family_precision * family_recall) / (family_precision + family_recall))

        family_false_negative_rate = 0 if (family_tp + family_fn) == 0 else family_fn / (family_fn + family_tp)  # Miss Rate
        # family_false_positive_rate = 0 if (family_fp + family_tn) == 0 else family_fp / (family_fp + family_tn)
        family_false_positive_rate = 0 if (family_fp + family_tp) == 0 else family_fp / (family_tp + family_fp)

        # family_specificity = 0 if (family_fp + family_tn) == 0 else family_tn / (family_fp + family_tn)  # True Negative Rate
        family_specificity = 1 - family_false_positive_rate

        family_false_discovery_rate = 0 if (family_fp + family_tp) == 0 else family_fp / (family_fp + family_tp)
        family_fowlkes_mallows_index = sqrt(family_precision * family_recall)

        # family_accuracy = (family_tp + family_tn) / (family_tp + family_tn + family_fp + family_fn)
        family_accuracy = 0 if (family_tp + family_fp + family_fn) == 0 else family_tp / (family_tp + family_fp + family_fn)

        result.append(family_tp)
        result.append(family_fp)
        result.append(family_fn)
        result.append(family_tn)
        result.append(family_precision)
        result.append(family_recall)
        result.append(family_specificity)
        result.append(family_f1_score)
        result.append(family_false_negative_rate)
        result.append(family_false_positive_rate)
        result.append(family_false_discovery_rate)
        result.append(family_fowlkes_mallows_index)
        result.append(family_accuracy)

        return result

    def search_expected_terms_local(self, protein_id: str, only_admit_exp_evidence_codes: bool):
        expected_go_terms = set()

        terms_sql = f"select group_concat(DISTINCT upg.go_term) as go_terms from uniprotJun16.uniprot_goa_20160607 upg where upg.protein_name = '{protein_id}'"

        connection = self.pool.get_connection()
        cursor = None

        try:
            with connection.cursor() as cursor:
                cursor.execute(terms_sql)

                for row in cursor:
                    if row is None:
                        expected_go_terms = set()
                    else:
                        expected_go_terms = set(row['go_terms'].split(',')) if row['go_terms'] is not None else set()
        except Exception as err:
            self.logger.error("An Exception in CACHE lookup occurred")
            self.logger.error(err)
            self.logger.error(traceback.format_exc())
            from utilities.emailmanager import EmailManager
            EmailManager.send_message('joseph.bonello@um.edu.mt', 'Prediction Error',
                                      "\r\n".join(["In predictions_cath.py", traceback.format_exc(), "", "Current protein id is ", protein_id]))
            return None, None, None, None, None, None, None, None

        finally:
            if cursor is not None:
                cursor.close()

            self.pool.close_connection(connection)

        return expected_go_terms

    def search_expected_terms_quickgo(self, protein_id: str, searchable_sequence: str, search: Searches, only_admit_exp_evidence_codes=bool):
        if only_admit_exp_evidence_codes:
            expected_go_terms = search.quickgo_terms_by_protein_id(protein_id=protein_id, only_admit_evidence_codes=self.evidenceCodes)
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
                    if only_admit_exp_evidence_codes:
                        expected_go_terms = search.quickgo_terms_by_protein_id(results[0], self.evidenceCodes)
                    else:
                        expected_go_terms = search.quickgo_terms_by_protein_id(results[0])
                else:
                    expected_go_terms = None
            else:
                self.logger.info("Expected GO Terms for protein: {}".format(expected_go_terms))

        return expected_go_terms

    def search_expected_terms_from_rollback(self, protein_id: str, only_admit_evidence_codes=bool):
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

        connection = self.pool.get_connection()
        cursor = None

        try:
            with connection.cursor() as cursor:
                cursor.execute(terms_sql)

                for row in cursor:
                    expected_go_terms.append(row['go_term'])
        except Exception as err:
            self.logger.error("An Exception in CACHE lookup occurred")
            self.logger.error(err)
            self.logger.error(traceback.format_exc())
            from utilities.emailmanager import EmailManager
            EmailManager.send_message('joseph.bonello@um.edu.mt', 'Prediction Error',
                                      "\r\n".join(["In predictions_cath.py", traceback.format_exc(), "", "Current protein id is ", protein_id]))
            return None, None, None, None, None, None, None, None

        finally:
            if cursor is not None:
                cursor.close()

            self.pool.close_connection(connection)

        return expected_go_terms

    def get_expected_go_terms(self, protein_id: str, searchable_sequence: Seq, search: Searches):
        expected_go_terms = set()
        bp_terms = set()
        mf_terms = set()
        cc_terms = set()

        # Step 1C: Get the expected Go Terms
        self.logger.debug("Step 1C")
        # expected_go_terms = set(search.go_terms_by_protein_id(protein_id))
        # expected_go_terms = expected_go_terms.union(set(search.go_terms_by_protein_id(protein_name)))
        if protein_id is not None:
            # if self.use_rollback:
            #     expected_go_terms = self.search_expected_terms_from_rollback(protein_id, only_admit_evidence_codes=self.evidenceCodes)
            # else:
            #     # expected_go_terms = self.search_expected_terms_quickgo(protein_id, str(searchable_sequence), search, only_admit_exp_evidence_codes=self.evidenceCodes)
            #     expected_go_terms = self.search_expected_terms_local(protein_id, only_admit_exp_evidence_codes=self.evidenceCodes)

            if expected_go_terms is None or len(expected_go_terms) == 0:
                expected_go_terms = self.search_expected_terms_quickgo(protein_id, searchable_sequence, search, only_admit_exp_evidence_codes=True)

            if expected_go_terms is not None:
                full_parental_terms = set(expected_go_terms)
                for term in expected_go_terms:
                    for pterm in self.go_tool.get_parental_terms(term):
                        full_parental_terms.add(pterm)

                expected_go_terms = full_parental_terms

            # Step 1D: Since we have expected terms, check what ontology they belong to
            self.logger.debug("Step 1D")

            if expected_go_terms is not None and len(expected_go_terms) > 0:
                term_ontologies = self.go_tool.get_ontology_for_term_list(expected_go_terms)

                for term in term_ontologies:
                    ttype = term_ontologies[term]
                    if ttype == 'biological_process':
                        bp_terms.add(term)
                    if ttype == 'molecular_function':
                        mf_terms.add(term)
                    if ttype == 'cellular_component':
                        cc_terms.add(term)

        return expected_go_terms, bp_terms, cc_terms, mf_terms

    def step0_identify_sequence(self, sequence: str) -> Tuple[str, str, Seq]:
        self.logger.debug("Step 0")

        fasta_handle = io.StringIO(sequence)

        target_id = None
        protein_description = None
        searchable_sequence = None

        for seq_record in SeqIO.parse(fasta_handle, "fasta"):
            target_id = seq_record.id if '|' not in seq_record.id else seq_record.id.split('|')[1]
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
            protein_id = protein_description.split(' ')[0].strip() if '|' not in protein_description else protein_description.split(' ')[0].split('|')[1]
            if '(' in protein_description:
                protein_name = protein_description.split(' ')[1].strip(')').strip('(').strip()
            else:
                protein_name = protein_description[protein_description.find(' '):].strip() if '|' not in protein_description else protein_description.split(' ')[0].split('|')[2]

            self.logger.warn("Warning: Cannot get details from Uniprot. Reverted to provided details.")

        # protein_gene_access = search.gene_access_from_protein_id(protein_id)

        self.logger.debug("Processing {}".format(protein_description))
        self.logger.debug("Processing {} (Blasted protein id {})".format(protein_description, protein_id))

        return protein_id, protein_name

    def step1B_funfhmmer(self, funfhmmer: Funfhmmer, protein_description: str, searchable_sequence: Seq, output: List):
        self.logger.debug("Step 1B")
        related_families = funfhmmer.fhmmer_search(protein_description, sequence=str(searchable_sequence), output=output)
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
            self.logger.error(f"ERROR: Funfhmmer returned None for: {protein_description}")
            self.logger.error(f"ERROR: Faulty sequence is: \n{str(searchable_sequence)}")

        return superfamilies, funfams

    def get_hashed_sequence(self, sequence: str):
        self.b2b.update(sequence.encode())
        return self.b2b.hexdigest()

    def target_in_cafa(self, target: str):
        sql_query = f"""SELECT Target FROM {self.phdPrefix}.CAFA3Targets WHERE Target=%s"""
        connection = self.pool.get_connection()
        cursor = None
        try:
            with connection.cursor() as cursor:
                cursor.execute(sql_query, (target,))

                row = cursor.fetchone()

                if row is not None:
                    return target in row.values()
        except Exception as err:
            self.logger.error("An Exception in CACHE lookup occurred")
            self.logger.error(err)
            self.logger.error(traceback.format_exc())
            from utilities.emailmanager import EmailManager
            EmailManager.send_message('joseph.bonello@um.edu.mt', 'Prediction Error',
                                      "\r\n".join(["In predictions_cath.py", traceback.format_exc(), "", "Current target is: ", target]))
            return None, None, None, None, None, None, None, None

        finally:
            if cursor is not None:
                cursor.close()

            self.pool.close_connection(connection)

        return False

    def predict(self, sequence: str, output: list, cafa_set: str, file_dict: dict, term_file_dict: dict):
        ontologies_list = ['BP', 'CC', 'MF']

        sql_insert_protein = f"""
                                REPLACE INTO {self.phdCachePrefix}.pred_!(hashed_sequence, protein_description, protein_id, protein_name, related_sfamilies,
                                related_ffamilies, expected_go_terms, expected_go_terms_bp, expected_go_terms_cc,
                                expected_go_terms_mf, cafa_target, Target) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s);
                            """

        # noinspection SQL error handling
        sql_insert_protein_rollback = f"""
                                            REPLACE INTO {self.phdCachePrefix}.pred_!(hashed_sequence, protein_description, protein_id, protein_name, related_sfamilies,
                                            related_ffamilies, expected_go_terms_rollback, expected_go_terms_rollback_bp, expected_go_terms_rollback_cc,
                                            expected_go_terms_rollback_mf,cafa_target,Target) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s);
                                        """

        search = Searches()
        funfhmmer = self.funfhmmer

        hashed_sequence = self.get_hashed_sequence(sequence)

        target_id, protein_description, searchable_sequence = self.step0_identify_sequence(sequence)

        # if not self.cath_predictor.target_in_cafa(target_id):
        #    return

        cache_connection = self.pool.get_connection(pool_name='CACHE')
        cache_cursor = None

        try:
            protein_id, protein_name, expected_go_terms, bp_terms, cc_terms, mf_terms, superfamilies, funfams = self.get_data_from_cache(cafa_set, hashed_sequence, sequence,
                                                                                                                                         funfhmmer, protein_description,
                                                                                                                                         searchable_sequence, output, search,
                                                                                                                                         skip_sanity_checks=False)

            if protein_id is not None or protein_name is not None:
                self.logger.debug("Step 1[DB]")

                self.logger.debug(f"Processing {protein_name} (Blasted protein id {protein_id}, Description {protein_description}) from cache")
                self.logger.debug(f"Superfamilies: <{superfamilies}>")
                self.logger.debug(f"Funfams: <{funfams}>")
                self.logger.debug(f"Expected Go Terms: <{expected_go_terms}>")
            else:  # Nothing found in cache
                self.logger.debug(f"Protein Target {protein_description} was not found in cache.")

                # Step 1A: Get updated protein_id and protein_name using BLAST (on SwissProt)
                protein_id, protein_name = self.step1A_protein_name_and_id(searchable_sequence, protein_description)

                # Step 1B: Funfhmmer the sequence
                superfamilies, funfams = self.step1B_funfhmmer(funfhmmer, protein_description, searchable_sequence, output)

                # Step 1C and 1D: Get the expected Go Terms (True Positive Set) and split them into respective
                # ontology
                expected_go_terms, bp_terms, cc_terms, mf_terms = self.get_expected_go_terms(protein_id, searchable_sequence, search)

                # Step 1E: Save the data into the cache
                self.logger.debug('Step 1E')
                cafa_target = protein_description.split(' ')[0] if 'sp' not in protein_description else protein_id
                protein_description = protein_description if 'sp' not in protein_description else f"{protein_id} {protein_description}"
                if expected_go_terms is None:
                    expected_go_terms = set()
                values = (
                    hashed_sequence, protein_description, protein_id, protein_name, str(','.join(superfamilies)), str(','.join(funfams)), str(','.join(expected_go_terms)),
                    str(','.join(bp_terms)), str(','.join(cc_terms)), str(','.join(mf_terms)), cafa_target, cafa_target)

                with cache_connection.cursor() as cache_cursor:
                    if self.use_rollback:
                        cache_cursor.execute(sql_insert_protein_rollback.replace("!", cafa_set), values)
                    else:
                        cache_cursor.execute(sql_insert_protein.replace("!", cafa_set), values)

            # Step 1F: Sanity Checks
            self.logger.debug('Step 1F')

            if superfamilies is None or funfams is None or len(superfamilies) == 0 or len(funfams) == 0:
                self.logger.warn("Cannot do predictions if superfamilies or funfams are None")
                return  # Cannot do predictions

            # Step 2: Perform the prediction
            self.logger.debug("Step 2")

            target_id = protein_description.split(' ')[0] if 'sp' not in protein_description else protein_id

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

                            # if len(expected_go_terms) == 0:  # There might not be terms in the ontology
                            #     continue

                            results = self.perform_prediction_dcgo(expected_go_terms,  # str(funfams).replace("{", "(")
                                                                   # .replace("}", ")"),
                                                                   list(funfams), ontology, term_file_funfams, protein_name, protein_id, target_id)

                            if results is None:
                                continue

                            # save to result file

                            if len(results) <= 2:
                                self.logger.debug(f"No results obtained for protein {protein_description}")
                            else:
                                self.logger.info(f"Prediction Results(dcGO,{ontology}): {','.join(str(v) for v in results)}")
                                self.write_to_pred_file(pred_file_funfams, results, target_id, protein_name, protein_id, ontology)

                        except Exception as err:
                            self.logger.error(err)
                            self.logger.error(traceback.format_exc())
                            from utilities.emailmanager import EmailManager
                            EmailManager.send_message('joseph.bonello@um.edu.mt', 'Prediction Error',
                                                      "\r\n".join(["In predictions_cath.py", traceback.format_exc(), "", "Current sequence is", sequence]))
                            continue
                else:
                    try:
                        term_file_funfams = term_file_dict[("dcgo", "FF")]

                        results = self.perform_prediction_dcgo(expected_go_terms, # str(funfams).replace("{", "(").replace("}", ")"),
                                                               list(funfams), None, term_file_funfams, protein_name, protein_id, target_id)

                        if results is None:
                            return

                        # save to result file

                        if len(results) <= 2:
                            self.logger.debug(f"No results obtained for protein {protein_description}")
                        else:
                            self.logger.info(f"Prediction Results(dcGO,all_ontologies): {','.join(str(v) for v in results)}")
                            self.write_to_pred_file(pred_file_funfams, results, target_id, protein_name, protein_id, '')

                    except Exception as err:
                        self.logger.error(err)
                        self.logger.error(traceback.format_exc())
                        from utilities.emailmanager import EmailManager
                        EmailManager.send_message('joseph.bonello@um.edu.mt', 'Prediction Error',
                                                  "\r\n".join(["In predictions_cath.py", traceback.format_exc(), "", "Current sequence is", sequence]))
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

                            # if len(expected_go_terms) == 0:  # There might not be terms in the ontology
                            #     continue

                            term_file_superfams = term_file_dict[(index, "SF")]
                            term_file_funfams = term_file_dict[(index, "FF")]

                            results = self.perform_prediction_set_based(index, expected_go_terms, list(superfamilies), list(funfams), '_samefunfam_similarity_allEC', ontology,
                                                                        term_file_superfams, term_file_funfams, protein_name, protein_id, target_id)

                            if results is None:
                                continue

                            # save to result file

                            results_superfam = results[0]
                            results_funfams = results[1]

                            pred_file_superfams = file_dict[(index, "SF")]

                            pred_file_funfams = file_dict[(index, "FF")]

                            if len(results_superfam) <= 2:
                                self.logger.debug(f"No results obtained for protein {protein_description} with Superfamilies")
                            else:
                                self.logger.info(f"Prediction Results({index},{ontology},SF): {','.join(str(v) for v in results_superfam)}")
                                self.write_to_pred_file(pred_file_superfams, results_superfam, target_id, protein_name, protein_id, ontology)

                            if len(results_funfams) <= 2:
                                self.logger.debug(f"No results obtained for protein {protein_description} with Funfams")
                            else:
                                self.logger.info(f"Prediction Results({index},{ontology},FF): {','.join(str(v) for v in results_funfams)}")
                                self.write_to_pred_file(pred_file_funfams, results_funfams, target_id, protein_name, protein_id, ontology)

                        except Exception as err:
                            self.logger.error(err)
                            self.logger.error(traceback.format_exc())
                            from utilities.emailmanager import EmailManager
                            EmailManager.send_message('joseph.bonello@um.edu.mt', 'Prediction Error',
                                                      "\r\n".join(["In predictions_cath.py", traceback.format_exc(), "", "Current sequence is", sequence]))
                            continue
                else:
                    try:
                        term_file_superfams = term_file_dict[(index, "SF")]
                        term_file_funfams = term_file_dict[(index, "FF")]

                        results = self.perform_prediction_set_based(index, expected_go_terms, list(superfamilies), list(funfams), '_samefunfam_similarity_allEC',
                                                                    None, term_file_superfams, term_file_funfams, protein_name, protein_id, target_id)

                        if results is None:
                            continue

                        # save to result file

                        results_superfam = results[0]
                        results_funfams = results[1]

                        pred_file_superfams = file_dict[(index, "SF")]

                        pred_file_funfams = file_dict[(index, "FF")]

                        if len(results_superfam) <= 2:
                            self.logger.debug(f"No results obtained for protein {protein_description} with Superfamilies")
                        else:
                            self.logger.info(f"Prediction Results({index},all_ontologies,SF): {','.join(str(v) for v in results_superfam)}")
                            self.write_to_pred_file(pred_file_superfams, results_superfam, target_id, protein_name, protein_id, '')

                        if len(results_funfams) <= 2:
                            self.logger.debug(f"No results obtained for protein {protein_description} with Funfams")
                        else:
                            self.logger.info(f"Prediction Results({index},all_ontologies,FF): {','.join(str(v) for v in results_funfams)}")
                            self.write_to_pred_file(pred_file_funfams, results_funfams, target_id, protein_name, protein_id, '')

                    except Exception as err:
                        self.logger.error(err)
                        self.logger.error(traceback.format_exc())
                        from utilities.emailmanager import EmailManager
                        EmailManager.send_message('joseph.bonello@um.edu.mt', 'Prediction Error',
                                                  "\r\n".join(["In predictions_cath.py", traceback.format_exc(), "", "Current sequence is", sequence]))
                        continue

        except Exception as err:
            self.logger.error(err)
            self.logger.error(traceback.format_exc())
            from utilities.emailmanager import EmailManager
            EmailManager.send_message('joseph.bonello@um.edu.mt', 'Prediction Error',
                                      "\r\n".join(["In predictions_cath.py", traceback.format_exc(), "", "Current sequence is", sequence]))

        finally:
            if cache_cursor is not None:
                cache_cursor.close()

            if cache_connection is not None:
                cache_connection.commit()
                self.pool.close_connection(cache_connection, pool_name="CACHE")

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

    def get_funfam_go_terms(self, ontology: str, separate_ontologies: bool, funfams: List, index: str, table_extension: str = None) -> (Set, List):
        connection = self.pool.get_connection()
        cursor_go_terms = None

        try:
            with connection.cursor() as cursor_go_terms:
                ontology_column, ontology_filter = self.ontology_database_filters(ontology, separate_ontologies, 't')

                # get_go_term_from_funfams = f"""
                #             SELECT dcgff.go_term AS goterm, dcgff.hscore, cffp.proportion as confidence {ontology_column}
                #             FROM {self.phdPrefix}.dcgoFunFams dcgff
                #             JOIN {self.phdPrefix}.CATHFFProportions cffp ON dcgff.go_term = cffp.go_term AND dcgff.cathfunfamilyfull_id = cffp.cathfunfamilyfull_id
                #             JOIN {self.phdPrefix}.GoGraph_Terms t on cffp.go_term = t.go_term
                #             WHERE dcgff.cathfunfamilyfull_id IN {str(funfams).replace('[', '(').replace(']', ')')}
                #             AND dcgff.hscore > 0 {ontology_filter}
                #             GROUP BY dcgff.go_term
                #         """

                if index == 'dcgo':
                    get_go_term_from_funfams = f"""
                                select dff.go_term AS goterm, sum(hscore) as sum {ontology_column}
                          from {self.phdPrefix}.dcgoFunFams dff
                          join {self.phdPrefix}.GoGraph_Terms t on dff.go_term = t.go_term
                          where dff.cathfunfamilyfull_id in {str(funfams).replace('[', '(').replace(']', ')')}
                            AND dff.hscore > 0 {ontology_filter}
                          group by dff.go_term;
                            """
                else:
                    get_go_term_from_funfams = f"""
                                            select DISTINCT
                                                cffs.goterm{ontology_column},
                                                cffs.{index} as score
                                            from
                                                {self.phdPrefix}.cathfunfamiliessimilarity{table_extension} cffs
                                            join {self.phdPrefix}.GoGraph_Terms t on cffs.goterm = t.go_term
                                            where cffs.cathfunfam in {str(funfams).replace('[', '(').replace(']', ')')}
                                            {ontology_filter}
                                            and cffs.{index} >= {self.threshold} -- default threshold 0.09
                                            order by cffs.{index}, cffs.goterm;
                                        """

                if index == 'dcgo':
                    get_go_terms_in_funfams = f"""
                                        SELECT DISTINCT dcgff.go_term as "go_term"
                                        FROM {self.phdPrefix}.dcgoFunFams dcgff
                                        WHERE dcgff.cathfunfamilyfull_id IN {str(funfams).replace('[', '(').replace(']', ')')}
                                    """
                else:
                    get_go_terms_in_funfams = f"""
                                            SELECT distinct cffs1.goterm as "go_term"
                                            FROM  {self.phdPrefix}.cathfunfamiliessimilarity{table_extension} cffs1
                                            WHERE cffs1.cathfunfam in {str(funfams).replace('[', '(').replace(']', ')')}
                                        """

                cursor_go_terms.execute(get_go_term_from_funfams)
                scored_funfams_go_terms = cursor_go_terms.fetchall()

                cursor_go_terms.execute(get_go_terms_in_funfams)
                all_gterms = cursor_go_terms.fetchall()

                all_funfams_go_terms = set(s['go_term'] for s in all_gterms)

                connection.commit()
        except Exception as err:
            self.logger.error(err)
            self.logger.error(traceback.format_exc())
            from utilities.emailmanager import EmailManager
            EmailManager.send_message('joseph.bonello@um.edu.mt', 'Prediction Error', "\r\n".join(["In predictions_cath.py", traceback.format_exc()]))

        finally:
            if cursor_go_terms is not None:
                cursor_go_terms.close()

            self.pool.close_connection(connection)

        return all_funfams_go_terms, scored_funfams_go_terms

    def get_superfamily_go_terms(self, ontology: str, separate_ontologies: List, superfamilies: List, index: str):
        cursor_go_terms = None
        connection = self.pool.get_connection()

        try:
            with connection.cursor() as cursor_go_terms:
                ontology_column, ontology_filter = self.ontology_database_filters(ontology, separate_ontologies, 't')

                superfamilies_formatted = str(superfamilies).replace('[', '(').replace(']', ')') if superfamilies is not None else '()'

                get_go_terms_from_superfamilies = f"""
                                        select cfs.goterm{ontology_column},
                                          cfs.{index} as score
                                        from {self.phdPrefix}.cathfamiliessimilarity cfs
                                        join {self.phdPrefix}.GoGraph_Terms t on cfs.goterm = t.go_term
                                        where cfs.cathfam in {superfamilies_formatted}
                                        and cfs.{index} >= ({self.threshold} * 0.01) {ontology_filter}
                                        order by cfs.{index};
                                    """  # Arbitrarily set due to high dispersion in superfamilies

                get_go_terms_in_superfamily = f"""
                                        SELECT distinct cfs1.goterm
                                        FROM  {self.phdPrefix}.cathfamiliessimilarity cfs1
                                        WHERE cfs1.cathfam in {str(superfamilies).replace('[', '(').replace(']', ')')};
                                    """

                cursor_go_terms.execute(get_go_terms_from_superfamilies)
                scored_superfamily_go_terms = cursor_go_terms.fetchall()

                cursor_go_terms.execute(get_go_terms_in_superfamily)
                all_superfamily_go_terms = set(s['goterm'] for s in cursor_go_terms.fetchall())

                connection.commit()
        except Exception as err:
            self.logger.error(err)
            self.logger.error(traceback.format_exc())
            from utilities.emailmanager import EmailManager
            EmailManager.send_message('joseph.bonello@um.edu.mt', 'Prediction Error', "\r\n".join(["In predictions_cath.py", traceback.format_exc()]))

        finally:
            if cursor_go_terms is not None:
                cursor_go_terms.close()

            self.pool.close_connection(connection)

            return all_superfamily_go_terms, scored_superfamily_go_terms

    def predict_from_funfams_dcgo(self, ontology: str, separate_ontologies: bool, funfams: list, enriched: bool = True):
        all_funfams_go_terms, scored_funfams_go_terms = self.get_funfam_go_terms(ontology, separate_ontologies, funfams, 'dcgo')

        if enriched:
            all_funfams_go_terms = self.enrich_go_terms(all_funfams_go_terms)

        hscores = {s['goterm']: s['sum'] for s in scored_funfams_go_terms}

        minhscore = min(hscores.values()) if len(hscores) > 0 else 0  # ensure there are predictions
        maxhscore = max(hscores.values()) if len(hscores) > 0 else 0
        denominator = maxhscore - minhscore
        pscores = {s['goterm']: abs((s['sum'] - minhscore) / denominator) if denominator > 0 else minhscore / denominator if denominator > 0 and s['sum'] == minhscore else 0 for s
                   in scored_funfams_go_terms}  # Ensures that the minhscore is not discarded

        pscores_kept = {s: pscores[s] for s in pscores.keys() if pscores[s] > 0.5}  # dcGO limit

        scored_funfams_go_terms = list()
        for s in pscores_kept:
            r = {'goterm': s}
            r['score'] = pscores_kept[s]
            if separate_ontologies:
                r['ontology'] = ontology
            scored_funfams_go_terms.append(r)

        funfams_predicted_go_terms = set(s['goterm'] for s in scored_funfams_go_terms)

        final_predicted_terms = funfams_predicted_go_terms
        parental_go_terms = None

        if len(funfams_predicted_go_terms) > 0:
            if enriched:
                final_predicted_terms = self.enrich_go_terms(funfams_predicted_go_terms)
                parental_go_terms = final_predicted_terms.difference(funfams_predicted_go_terms)

        return scored_funfams_go_terms, final_predicted_terms, parental_go_terms, all_funfams_go_terms

    def perform_prediction_dcgo(self, expected_go_terms: Set, funfams: Set, ontology: str, term_file_funfams: TextIO, protein_name: str, protein_id: str, target_id: str):
        result = list()

        # Preliminary part: Get Expected Terms
        if self.inherit_go_terms:
            final_expected_terms = self.enrich_go_terms(expected_go_terms)
        else:
            final_expected_terms = expected_go_terms

        result.append(len(final_expected_terms))

        # Get Funfam terms and scores

        predictions = self.predict_from_funfams_dcgo(ontology, self.separate_ontologies, funfams, self.inherit_go_terms)

        scored_funfams_go_terms = predictions[0]
        final_predicted_terms = predictions[1]
        parental_go_terms = predictions[2]
        all_funfams_go_terms = predictions[3]

        result.append(len(final_predicted_terms) if final_predicted_terms is not None else 0)

        if final_predicted_terms is not None and len(final_predicted_terms) > 0:
            for res in self.calculate_protein_prediction_metrics(final_expected_terms, final_predicted_terms, all_funfams_go_terms):
                result.append(res)

            self.write_to_term_file(term_file_funfams, expected_go_terms, final_expected_terms, scored_funfams_go_terms, final_predicted_terms, target_id, protein_name, protein_id,
                                    ontology, parental_go_terms=parental_go_terms, is_scored_pred_terms=True)
        else:
            result = None  # No predictions could be made

        return result

    def predict_from_superfamilies_set_based(self, ontology: str, separate_ontologies: bool, superfamilies: list, index: str, enriched: bool = True):
        all_superfamily_go_terms, scored_superfamily_go_terms = self.get_superfamily_go_terms(ontology, separate_ontologies, superfamilies, index)

        if enriched:
            all_superfamily_go_terms = self.enrich_go_terms(all_superfamily_go_terms)

        superfamily_predicted_go_terms = set(s['goterm'] for s in scored_superfamily_go_terms)

        parental_go_terms = None
        final_predicted_terms = superfamily_predicted_go_terms

        if len(superfamily_predicted_go_terms) > 0:
            if enriched:
                final_predicted_terms = self.enrich_go_terms(superfamily_predicted_go_terms)
                parental_go_terms = final_predicted_terms.difference(superfamily_predicted_go_terms)

        return scored_superfamily_go_terms, final_predicted_terms, parental_go_terms, all_superfamily_go_terms

    def predict_from_funfams_set_based(self, table_extension: str, ontology: str, separate_ontologies: bool, funfams: list, index: str, enriched: bool = True):
        all_funfams_go_terms, scored_funfams_go_terms = self.get_funfam_go_terms(ontology, separate_ontologies, funfams, index, table_extension)

        if enriched:
            all_funfams_go_terms = self.enrich_go_terms(all_funfams_go_terms)

        funfams_predicted_go_terms = set(s['goterm'] for s in scored_funfams_go_terms)

        final_predicted_terms = funfams_predicted_go_terms
        parental_go_terms = None

        if len(funfams_predicted_go_terms) > 0:
            if enriched:
                final_predicted_terms = self.enrich_go_terms(funfams_predicted_go_terms)
                parental_go_terms = final_predicted_terms.difference(funfams_predicted_go_terms)

        return scored_funfams_go_terms, final_predicted_terms, parental_go_terms, all_funfams_go_terms

    def perform_prediction_set_based(self, index: str, expected_go_terms: Set, superfamilies: Set, funfams: Set, tableExtension: str, ontology: str, term_file_superfams: TextIO,
                                     term_file_funfams: TextIO, protein_name: str, protein_id: str, target_id: str, include_superfamily_predictions: bool = False):
        # Results for superfamilies (0) and funfams (1)
        results = (list(), list())

        # Preliminary part: Get Expected Terms
        if self.inherit_go_terms:
            final_expected_terms = self.enrich_go_terms(expected_go_terms)
        else:
            final_expected_terms = expected_go_terms

        results[0].append(len(final_expected_terms))
        results[1].append(len(final_expected_terms))

        # Part 1; Get Go Terms from SuperFamilies
        result = results[0]

        if include_superfamily_predictions:
            predictions = self.predict_from_superfamilies_set_based(ontology, self.separate_ontologies, superfamilies, index, self.inherit_go_terms)
            scored_superfamily_go_terms = predictions[0]
            final_predicted_terms = predictions[1]
            parental_go_terms = predictions[2]
            all_superfamily_go_terms = predictions[3]

            result.append(len(final_predicted_terms))

            if final_predicted_terms is not None and len(final_predicted_terms) > 0:
                for res in self.calculate_protein_prediction_metrics(final_expected_terms, final_predicted_terms, all_superfamily_go_terms):
                    result.append(res)

                self.write_to_term_file(term_file_superfams, expected_go_terms, final_expected_terms, scored_superfamily_go_terms, final_predicted_terms, target_id, protein_name,
                                        protein_id, ontology, parental_go_terms=parental_go_terms, is_scored_pred_terms=True)
            else:
                result.append(list())
        else:
            result.append(0)  # No predicted terms
            result.append(list())  # No calculations possible

        # Part 2: Get Go Terms from Funfams

        result = results[1]

        predictions = self.predict_from_funfams_set_based(tableExtension, ontology, self.separate_ontologies, funfams, index, self.inherit_go_terms)

        scored_funfams_go_terms = predictions[0]
        final_predicted_terms = predictions[1]
        parental_go_terms = predictions[2]
        all_funfams_go_terms = predictions[3]

        result.append(len(final_predicted_terms) if final_predicted_terms is not None else 0)

        if final_predicted_terms is not None and len(final_predicted_terms) > 0:
            for res in self.calculate_protein_prediction_metrics(final_expected_terms, final_predicted_terms, all_funfams_go_terms):
                result.append(res)

            self.write_to_term_file(term_file_funfams, expected_go_terms, final_expected_terms, scored_funfams_go_terms, final_predicted_terms, target_id, protein_name, protein_id,
                                    ontology, parental_go_terms=parental_go_terms, is_scored_pred_terms=True)
        else:
            result.append(list())

        return results

    def predict_from_path(self, path, output, parallel_backend='threading'):
        import os
        for file in os.listdir(path):
            if not os.path.isdir(os.path.join(path, file)):
                print(file)

                sequences = set()
                sequence = ''

                import re
                regex = re.compile('>sp\|([A-Z0-9]+)\|([A-Z0-9_]+).*')

                with open(os.path.join(path, file), 'r') as file_handle:
                    for line in file_handle:
                        if line.strip() == '' or line.find('>') > -1:  # Handle new sequences; either through a line break
                            if '>' in line:
                                if '|' in line:
                                    line = regex.sub('> \1 \2', line.strip())
                            # or through the discovery of a new sequence right after this one ends
                            if sequence.strip() != '':
                                sequences.add(sequence)  # Checks if empty
                                sequence = ''
                            if line.find('>') > -1:
                                if sequence != '':
                                    sequences.add(sequence)
                                if '>' in line:
                                    if '|' in line:
                                        line = regex.sub('> \1 \2', line.strip())
                                sequence = line
                        else:
                            sequence += line.strip()

                if sequence != '':
                    sequences.add(sequence)  # set will handle doubles; avoids mossing last sequence in the file

                # self.logger.debugsequences)

                file_dict = dict()
                term_file_dict = dict()

                cafa_edition = os.getenv('CAFA_EDITION')
                if os.getenv('PREDICTION_NAMING_SCHEME') == 'CAFA2':
                    cafa_set = file[0:file.find('.')]
                else:
                    cafa_set = f"cafa{cafa_edition}_{file[0:file.rfind('.')].replace('.', '_').lower()}"

                import os
                import socket
                tmp_path = os.getenv('PREDICTION_SAVE_PATH_' + socket.gethostname().upper())

                cache_cursor = None
                cache_connection = self.pool.get_connection(pool_name='CACHE')

                try:

                    for index in self.indexes_to_use:
                        if not os.path.exists(os.path.join(tmp_path, "predictions", self.ts)):
                            os.makedirs(os.path.join(tmp_path, "predictions", self.ts))

                        if index != 'dcgo':
                            pred_file_superfams = open(
                                os.path.join(tmp_path, "predictions", self.ts, "prediction-{0}-{1}-t{2}-SF.csv".format(index, cafa_set, str(self.threshold).replace('.', '_'))), "w")
                        pred_file_funfams = open(
                            os.path.join(tmp_path, "predictions", self.ts, "prediction-{0}-{1}-t{2}-FF.csv".format(index, cafa_set, str(self.threshold).replace('.', '_'))), "w")

                        if index != 'dcgo':
                            term_file_superfams = open(
                                os.path.join(tmp_path, "predictions", self.ts, "terms-{0}-{1}-t{2}-SF.csv".format(index, cafa_set, str(self.threshold).replace('.', '_'))), "w")
                        term_file_funfams = open(
                            os.path.join(tmp_path, "predictions", self.ts, "terms-{0}-{1}-t{2}-FF.csv".format(index, cafa_set, str(self.threshold).replace('.', '_'))), "w")

                        if index != 'dcgo':
                            file_dict[(index, 'SF')] = pred_file_superfams
                        file_dict[(index, 'FF')] = pred_file_funfams

                        if index != 'dcgo':
                            term_file_dict[(index, 'SF')] = term_file_superfams
                        term_file_dict[(index, 'FF')] = term_file_funfams

                        if index != 'dcgo':
                            pred_file_superfams.write("CAFA Target,Protein Name,Protein Unique Id,{}expected_count,superfamily_predicted_count,superfamily_tp,"
                                                      "superfamily_fp,superfamily_fn,superfamily_tn,"
                                                      "superfamily_precision,superfamily_recall,superfamily_specificity,superfamily_f1_score,"
                                                      "superfamily_false_negative_rate,superfamily_false_positive_rate,superfamily_false_discovery_rate,"
                                                      "superfamily_fowlkes_mallows_index,superfamily_accuracy".format("Ontology, " if self.separate_ontologies else ""))

                        pred_file_funfams.write("CAFA Target,Protein Name,Protein Unique Id,{}expected_count,funfams_predicted_count,funfams_tp,"
                                                "funfams_fp,funfams_fn,funfams_tn,funfams_precision,funfams_recall,"
                                                "funfams_specificity,funfams_f1_score,funfams_false_negative_rate,funfams_false_positive_rate,"
                                                "funfams_false_discovery_rate, funfams_fowlkes_mallows_index,funfams_accuracy".format("Ontology, " if self.separate_ontologies else ""))

                        if index != 'dcgo':
                            term_file_superfams.write("CAFA Target,Protein Name,Protein Unique Id,{}GO Term,Score,Is From TPS,"
                                                      "Is Inherited From Parents".format("Ontology, " if self.separate_ontologies else ""))

                        term_file_funfams.write("CAFA Target,Protein Name,Protein Unique Id,{}GO Term,Score,Is From TPS,"
                                                "Is Inherited From Parents".format("Ontology, " if self.separate_ontologies else ""))

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

                    self.create_cache_tables(cache_connection, cafa_set)

                    import os
                    max_threads = int(os.getenv('NUM_THREADS'))

                    sequences = sorted(list(sequences))

                    Parallel(n_jobs=max_threads, verbose=5, backend=parallel_backend)(
                        delayed(unwrap_self_predict)(sequence, output=output, cafa_set=cafa_set, file_dict=file_dict, term_file_dict=term_file_dict) for sequence in
                        zip([self] * len(sequences), sequences))
                except Exception as err:
                    self.logger.error(err)
                    self.logger.error(traceback.format_exc())
                    from utilities.emailmanager import EmailManager
                    EmailManager.send_message('joseph.bonello@um.edu.mt', 'Prediction Error', "\r\n".join(["In predictions_cath.py", traceback.format_exc()]))

                finally:
                    self.logger.info("In finally ...")

                    for index in file_dict:
                        text_file = file_dict[index]
                        text_file.flush()
                        text_file.close()

                    for index in term_file_dict:
                        term_file = term_file_dict[index]
                        term_file.flush()
                        term_file.close()

                    if cache_cursor is not None:
                        cache_cursor.close()

                    self.pool.close_connection(cache_connection, pool_name='CACHE')

    def sanity_checks(self, cafa_set: str, protein_id: str, protein_name: str, protein_description: str, sequence: str, searchable_sequence: Seq, funfhmmer: Funfhmmer,
                      superfamilies: Set, funfams: Set, expected_go_terms: Set, bp_terms: Set, cc_terms: Set, mf_terms: Set, search: Searches, output: List,
                      do_not_retry_fetches: bool = False):

        cache_cursor = None
        cache_connection = self.pool.get_connection(pool_name='CACHE')
        values = tuple()

        try:
            if protein_id is None or protein_name is None:
                protein_id = protein_description.split(' ')[0]
                protein_name = protein_description.split(' ')[1]

            # Sanity Checks
            with cache_connection.cursor() as cache_cursor:
                if superfamilies is None or funfams is None or len(superfamilies.union(funfams)) == 0:
                    self.logger.warn(f"No CATH family records found returned None for: {protein_id} - {protein_name} - {protein_description}")
                    self.logger.warn("Faulty sequence is: \n{}".format(sequence))

                    if not do_not_retry_fetches:
                        superfamilies, funfams = self.step1B_funfhmmer(funfhmmer, protein_description, searchable_sequence, output)

                    if superfamilies is None or funfams is None or len(superfamilies.union(funfams)) == 0:
                        self.logger.warn("Retry to get CATH Funfams failed. Aborting this prediction.")
                        return False,
                    else:
                        values = (','.join(superfamilies), ','.join(funfams), protein_description)
                        cache_cursor.execute(self.sql_update_funfams.replace('!', cafa_set), values)

                expected_terms_needs_db_update = False
                if expected_go_terms is None or len(expected_go_terms) == 0:
                    self.logger.warn("Finding GO Terms for protein {}({}) failed, hence retrying".format(protein_id, protein_name))

                    # Try repeating Step 1C - and update the record appropriately
                    if not do_not_retry_fetches:
                        expected_go_terms, bp_terms, cc_terms, mf_terms = self.get_expected_go_terms(protein_id, searchable_sequence, search)

                    if expected_go_terms is None or len(expected_go_terms) == 0:  # If it is still 0
                        self.logger.warn("No experimental GO Terms found for protein {}({}), hence skipped".format(protein_id, protein_name))
                        # return False, These are not needed
                    else:
                        expected_terms_needs_db_update = True

                new_bp_terms = set()
                new_cc_terms = set()
                new_mf_terms = set()

                bp_terms = set() if bp_terms is None else bp_terms
                cc_terms = set() if cc_terms is None else cc_terms
                mf_terms = set() if mf_terms is None else mf_terms

                if bp_terms is None or cc_terms is None or mf_terms is None or len(bp_terms) == 0 or \
                        len(cc_terms) == 0 or len(mf_terms) == 0 or (len(bp_terms) + len(cc_terms) + len(mf_terms)) != len(expected_go_terms):
                    term_ontologies = self.go_tool.get_ontology_for_term_list(expected_go_terms)

                    for term in term_ontologies:
                        ttype = term_ontologies[term]
                        if ttype == 'biological_process':
                            new_bp_terms.add(term)
                        if ttype == 'molecular_function':
                            new_mf_terms.add(term)
                        if ttype == 'cellular_component':
                            new_cc_terms.add(term)

                    if len(new_bp_terms) != len(bp_terms) or len(new_cc_terms) != len(cc_terms) or len(new_mf_terms) != len(mf_terms):  # Only update if different
                        expected_terms_needs_db_update = True

                if expected_terms_needs_db_update:
                    # update experimental terms
                    e_terms = '' if expected_go_terms is None else ','.join(expected_go_terms)
                    b_terms = '' if bp_terms is None else ','.join(bp_terms)
                    c_terms = '' if cc_terms is None else ','.join(cc_terms)
                    m_terms = '' if mf_terms is None else ','.join(mf_terms)
                    p_desc = '' if protein_description is None else protein_description
                    protein_description = '' if protein_description is None else protein_description
                    values = (e_terms, b_terms, c_terms, m_terms, p_desc,)
                    cache_cursor.execute(self.sql_update_expected_terms.replace("!", cafa_set), values)
                    cache_connection.commit()

        except Exception as err:
            self.logger.error(err)
            self.logger.error(traceback.format_exc())
            from utilities.emailmanager import EmailManager
            e_terms = '' if expected_go_terms is None else ','.join(expected_go_terms)
            b_terms = '' if bp_terms is None else ','.join(bp_terms)
            c_terms = '' if cc_terms is None else ','.join(cc_terms)
            m_terms = '' if mf_terms is None else ','.join(mf_terms)
            p_desc = '' if protein_description is None else protein_description

            self.logger.error(f"Expected Terms = {e_terms}, BP Terms = {b_terms}, CC Terms = {c_terms}, MF Terms = {m_terms}, Protein Desc = {p_desc}")

            EmailManager.send_message('joseph.bonello@um.edu.mt', 'Prediction Error',
                                      "\r\n".join(["In predictions_cath.py", traceback.format_exc(), "-" * 20, e_terms, b_terms, c_terms, m_terms, p_desc, "-" * 20]))

        finally:
            if cache_cursor is not None:
                cache_cursor.close()

            self.pool.close_connection(cache_connection, pool_name='CACHE')

        sfams = set() if superfamilies is None else set(superfamilies)
        ffams = set() if funfams is None else set(funfams)
        exterms = set() if expected_go_terms is None else set(expected_go_terms)
        bterms = set() if bp_terms is None else set(bp_terms)
        cterms = set() if bp_terms is None else set(cc_terms)
        mterms = set() if bp_terms is None else set(mf_terms)
        return True, sfams, ffams, exterms, bterms, cterms, mterms

    def enrich_go_terms(self, go_term_set: Set):
        if self.inherit_go_terms:
            gold_standard_parental_go_terms = self.go_tool.get_parental_terms_for_list(list(go_term_set))

            gspgt = set()
            for pTerm in gold_standard_parental_go_terms:
                if pTerm not in ['all', 'GO:0005575', 'GO:0008150', 'GO:0003674']:  # Remove GO root terms
                    gspgt.add(pTerm)

            enriched_term_set = set(go_term_set.union(gspgt))
        else:
            enriched_term_set = go_term_set

        return enriched_term_set

    def calculate_probability(self, term: str):
        connection = self.pool.get_connection()
        cursor_go_terms = None
        score = 0.0

        try:
            with connection.cursor() as cursor_go_terms:
                cursor_go_terms.execute(self.sql_term_probability, term)

                score_res = cursor_go_terms.fetchone()
                score = score_res['probability'] if score_res is not None else 0
        except Exception as err:
            self.logger.error(err)
            self.logger.error(traceback.format_exc())
            from utilities.emailmanager import EmailManager
            EmailManager.send_message('joseph.bonello@um.edu.mt', 'Prediction Error', "\r\n".join(["In predictions_cath.py", traceback.format_exc()]))

        finally:
            if cursor_go_terms is not None:
                cursor_go_terms.close()

            self.pool.close_connection(connection)

        return score

    def write_to_term_file(self, term_file: TextIO, expected_go_terms: Set, final_expected_terms: Set, pred_terms: Set, final_pred_terms: Set, target_id: str, protein_name: str,
    protein_id: str, ontology: str, parental_go_terms: Set = None, is_scored_pred_terms: bool = False):
        if len(expected_go_terms) > 0:
            for expTerm in expected_go_terms:  # List of expected terms
                outstr = f'{target_id},{protein_name},{protein_id},{ontology},'
                outstr += f"{expTerm},1.0,Y,N\n"
                term_file.write(outstr)

            for pTerm in final_expected_terms.difference(expected_go_terms):  # List of parental expected terms
                outstr = f'{target_id},{protein_name},{protein_id},{ontology if ontology is not None else ""},'
                outstr += f"{pTerm},1.0,Y,Y\n"
                term_file.write(outstr)

        if len(pred_terms) > 0:
            if is_scored_pred_terms:
                for pTerm in pred_terms:  # Scored Predicted Terms
                    if float(pTerm['score']) > 0.0:
                        outstr = f"{target_id},{protein_name},{protein_id},{ontology if ontology is not None else ''},"
                        outstr += f"{pTerm['goterm']},{pTerm['score']},N,N\n"
                        term_file.write(outstr)

                if parental_go_terms is not None:
                    for pTerm in parental_go_terms:
                        if pTerm not in ['all', 'GO:0005575', 'GO:0008150', 'GO:0003674']:  # Remove GO root terms
                            score = self.calculate_probability(pTerm)

                            if float(score) > 0.0:
                                outstr = f"{target_id},{protein_name},{protein_id},{ontology if ontology is not None else ''},"
                                outstr += f"{pTerm},{score},N,Y\n"

                                term_file.write(outstr)
            else:
                for pred_term in pred_terms:  # Predicted Terms
                    score = self.cath_predictor.calculate_probability(pred_term)
                    if float(score) > 0.0:
                        outstr = f'{target_id},{protein_name},{protein_id},{ontology},'
                        outstr += f"{pred_term},{score},N,N\n"
                        term_file.write(outstr)

                for pTerm in final_pred_terms.difference(pred_terms): # Parental GO Terms for the predicted
                    score = self.cath_predictor.calculate_probability(pTerm)
                    if float(score) > 0.0:
                        outstr = f'{target_id},{protein_name},{protein_id},{ontology if ontology is not None else ""},'
                        outstr += f"{pTerm},{score},N,Y\n"
                        term_file.write(outstr)

            term_file.flush()

    def write_to_pred_file(self, pred_file_funfams: TextIO, results: List, target_id: str, protein_name: str, protein_id: str, ontology: str):
        if (len(results) > 0):
            pred_file_funfams.write(f"{target_id},{protein_name},{protein_id},{ontology},")
            pred_file_funfams.write(",".join(str(v) for v in results))
            pred_file_funfams.write("\n")

            pred_file_funfams.flush()

    def create_cache_tables(self, cache_connection: pymysql.Connection, cafa_set: str):
        with cache_connection.cursor() as cache_cursor:
            # Create if the cache table does not exist
            sql_check_table = f"""
                CREATE TABLE IF NOT EXISTS {self.phdCachePrefix}.pred_!
                (
                    id INT PRIMARY KEY NOT NULL AUTO_INCREMENT,
                    hashed_sequence VARCHAR(128),
                    protein_description VARCHAR(350),
                    protein_id VARCHAR(50),
                    protein_name VARCHAR(50),
                    related_sfamilies MEDIUMTEXT,
                    related_ffamilies MEDIUMTEXT,
                    expected_go_terms MEDIUMTEXT,
                    expected_go_terms_bp MEDIUMTEXT,
                    expected_go_terms_cc MEDIUMTEXT,
                    expected_go_terms_mf MEDIUMTEXT
                ) ROW_FORMAT = COMPRESSED;
            """.replace("!", cafa_set)

            sql_rollback_expected_results = f"""
                ALTER TABLE {self.phdCachePrefix}.pred_!
                ADD COLUMN IF NOT EXISTS expected_go_terms_rollback MEDIUMTEXT,
                ADD COLUMN IF NOT EXISTS expected_go_terms_rollback_bp MEDIUMTEXT,
                ADD COLUMN IF NOT EXISTS expected_go_terms_rollback_cc MEDIUMTEXT,
                ADD COLUMN IF NOT EXISTS expected_go_terms_rollback_mf MEDIUMTEXT,
                ADD COLUMN IF NOT EXISTS cafa_target VARCHAR(50),
                ADD COLUMN IF NOT EXISTS Target VARCHAR(50);
            """.replace("!", cafa_set)

            # sql_constraint = f"""
            #     CREATE UNIQUE INDEX IF NOT EXISTS pred_!_hashed_sequence_uindex ON {self.phdCachePrefix}.pred_! (hashed_sequence);
            # """.replace("!", cafa_set)

            # sql_constraint_target = f"""CREATE UNIQUE INDEX IF NOT EXISTS pred_!_hashed_Target_uindex ON {self.phdCachePrefix}.pred_! (Target);""".replace("!", cafa_set)

            cache_cursor.execute(sql_check_table)
            cache_cursor.execute(sql_rollback_expected_results)
            # cache_cursor.execute(sql_constraint)
            # cache_cursor.execute(sql_constraint_target)
