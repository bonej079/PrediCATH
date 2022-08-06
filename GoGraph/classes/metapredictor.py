import io
import os
import traceback
from typing import List, Set, Tuple

from joblib import Parallel, delayed

from GoGraph.classes.decorators import *
from GoGraph.classes.geneontology import GeneOntologyUtility
from GoGraph.classes.searches import Searches
from utilities.mysqlconnectionpool import MySQLConnectionPool


def unwrap_self_predict(arg, **kwarg):
    return MetaPredictor.predict_a_sequence(*arg, **kwarg)


# noinspection SqlResolve,DuplicatedCode
class MetaPredictor:
    # indexes_to_use = []

    __slots__ = ['threshold', 'inherit_go_terms', 'separate_ontologies', 'indexes_to_use', 'ts', 'evidenceCodes', 'go_tool', 'logger', 'db_connection', 'sql_term_probability',
                 'use_rollback', 'ontologies', 'cath_predictor', 'funfhmmer', 'skip_not_in_cache', 'target_logger', 'notfound_logger', 'phdPrefix', 'use_database_for_output',
                 'pred_file_dict', 'term_file_dict', 'cafa_file_dict', 'aggressive_debug', 'aggressive_debug_dict', 'aggressive_debug_pred_dict', 'cafa_edition', 'cafa_set',
                 'cafa_species', 'cafa_group']

    def __init__(self, indexes_to_use: list, ontologies: list, threshold: float, inherit_go_terms: bool, skip_not_in_cache: bool = True, log_to_file: bool = True,
                 aggressive_debug: bool = True):
        self.indexes_to_use = indexes_to_use
        self.ontologies = ontologies
        self.logger = create_logger(loglevel=logging.DEBUG)

        if log_to_file:
            self.target_logger = create_custom_logger(loglevel=logging.DEBUG, name="sequence_logger",
                                                      path=os.path.join(os.getenv('PREDICTION_SAVE_PATH_UOM-1A26'), 'predictions', 'sequence.log'))
            self.notfound_logger = create_custom_logger(loglevel=logging.DEBUG, name="notfound_logger",
                                                        path=os.path.join(os.getenv('PREDICTION_SAVE_PATH_UOM-1A26'), 'predictions', 'not_found.log'))
            self.use_database_for_output = False
        else:
            self.use_database_for_output = True
            self.target_logger = create_custom_logger(loglevel=logging.DEBUG, name="sequence_logger", path=None)
            self.notfound_logger = create_custom_logger(loglevel=logging.DEBUG, name="notfound_logger", path=None)

        from GoGraph.classes.predictions_cath import Predictions_CATH
        self.cath_predictor = Predictions_CATH(threshold, inherit_go_terms=inherit_go_terms, separate_ontologies=True, indexes_to_use=indexes_to_use, log_to_file=log_to_file)

        from GoGraph.classes.funfhmmer import Funfhmmer
        self.funfhmmer = Funfhmmer()

        self.db_connection = MySQLConnectionPool.get_instance().get_connection()

        import datetime
        self.ts = datetime.datetime.now().strftime("%Y%m%d%H%M%S")

        self.threshold = threshold

        self.separate_ontologies = True

        self.go_tool = GeneOntologyUtility()

        self.inherit_go_terms = inherit_go_terms

        self.skip_not_in_cache = skip_not_in_cache

        self.phdPrefix = os.getenv("PHD_DB_PREFIX")

        self.cafa_group = os.getenv('CAFA_TEAM_NAME')

        self.pred_file_dict = None
        self.term_file_dict = None
        self.cafa_file_dict = None
        self.aggressive_debug_dict = None
        self.aggressive_debug_pred_dict = None

        self.aggressive_debug = aggressive_debug

        self.cafa_edition = None
        self.cafa_set = None
        self.cafa_species = None

    def __exit__(self, exc_type, exc_val, exc_tb):
        MySQLConnectionPool.get_instance().close_connection(self.db_connection)

    @exception(create_logger())
    def predict_by_funfhmmer(self, target: str):
        ret_val_bp = set()
        ret_val_cc = set()
        ret_val_mf = set()
        scoredFunFamerPredictions = list()

        db_cursor = None
        oracle_db_connection = None

        local_logger = self.logger

        try:
            oracle_db_connection = MySQLConnectionPool.get_instance().get_connection(pool_name='CACHE')
            with oracle_db_connection.cursor() as db_cursor:

                db_cursor.callproc(f"{self.phdPrefix}.funfhmmer_predictions", (target.strip(),))
                rows = db_cursor.fetchall()

                for row in rows:
                    if row['term_type'] in ['BP', 'biological_process']:
                        ret_val_bp.add(row['go_term'])
                    if row['term_type'] in ['CC', 'cellular_component']:
                        ret_val_cc.add(row['go_term'])
                    if row['term_type'] in ['MF', 'molecular_function']:
                        ret_val_mf.add(row['go_term'])
                    scoredFunFamerPredictions.append({'goterm': row['go_term'], 'score': row['confidence']})
        except Exception as ex:
            local_logger.error(ex)
            self.logger.error(traceback.format_exc())
            from utilities.emailmanager import EmailManager
            EmailManager.send_message('joseph.bonello@um.edu.mt', 'Prediction Error in MetaPredictor', "\r\n".join(["In metapredictor.py", traceback.format_exc()]))
        finally:
            if db_cursor is not None:
                db_cursor.close()

            if oracle_db_connection is not None:
                MySQLConnectionPool.get_instance().close_connection(oracle_db_connection, pool_name='CACHE')

        return ret_val_bp, ret_val_cc, ret_val_mf, scoredFunFamerPredictions

    def write_to_term_file(self, term_file_key: Tuple, expected_go_terms: Set, final_expected_terms: Set, pred_terms: Set, final_pred_terms: Set, target_id: str, protein_name: str,
                           protein_id: str, ontology: str):

        predictions_schema = os.getenv('PREDICTIONS_SCHEMA')

        sql_write_to_terms_table = f"""
                                                        INSERT INTO `{predictions_schema}`.`terms-{term_file_key}-cafa{self.cafa_edition}_sp_species_{self.cafa_species}-t
                                                        {str(self.threshold).replace('.', '_')}-FF` (
                                                            `CAFA_Target`,
                                                            `Protein_Name`,
                                                            `Protein_Unique_Id`,
                                                            `Ontology`,
                                                            `GO_Term`,
                                                            `Score`,
                                                            `Is_From_TPS`,
                                                            `Is_Inherited_From_Parents`)
                                                        VALUES(%s, %s, %s, %s, %s, %s, %s, %s);
                                                    """

        db_connection = None

        if len(expected_go_terms) > 0:
            if self.use_database_for_output:
                try:
                    db_connection = MySQLConnectionPool.get_instance().get_connection(pool_name='CACHE')
                    with db_connection.cursor() as db_cursor:
                        for expTerm in expected_go_terms:  # List of expected terms
                            db_cursor.execute(sql_write_to_terms_table, (target_id, protein_name, protein_id, ontology, expTerm, 1.0, 1, 0,))

                        for pTerm in final_expected_terms.difference(expected_go_terms):  # List of parental expected terms
                            db_cursor.execute(sql_write_to_terms_table, (target_id, protein_name, protein_id, ontology if ontology is not None else "", pTerm, 1.0, 1, 1))
                finally:
                    if db_connection is not None:
                        MySQLConnectionPool.get_instance().close_connection(db_connection, pool_name='CACHE')
            else:
                term_file = self.term_file_dict[term_file_key] if term_file_key in self.term_file_dict else self.aggressive_debug_dict[term_file_key]
                if term_file is None:
                    logger.error(f"Term file is None. Term file dict is {self.term_file_dict.keys()}; Aggressive Debug Dict is {self.aggressive_debug_dict.keys()}")
                for expTerm in expected_go_terms:  # List of expected terms
                    outstr = f'{target_id},{protein_name},{protein_id},{ontology},'
                    outstr += f"{expTerm},1.0,Y,N\n"
                    term_file.write(outstr)

                for pTerm in final_expected_terms.difference(expected_go_terms):  # List of parental expected terms
                    outstr = f'{target_id},{protein_name},{protein_id},{ontology if ontology is not None else ""},'
                    outstr += f"{pTerm},1.0,Y,Y\n"
                    term_file.write(outstr)

        if len(pred_terms) > 0:
            if self.use_database_for_output:
                try:
                    db_connection = MySQLConnectionPool.get_instance().get_connection(pool_name='CACHE')
                    with db_connection.cursor() as db_cursor:
                        for pred_term in pred_terms:  # Predicted Terms
                            score = self.cath_predictor.calculate_probability(pred_term)
                            if float(score) > 0.0:
                                db_cursor.execute(sql_write_to_terms_table, (target_id, protein_name, protein_id, ontology, pred_term, score, 0, 0,))

                        for pTerm in final_pred_terms.difference(pred_terms):
                            score = self.cath_predictor.calculate_probability(pTerm)
                            if float(score) > 0.0:
                                db_cursor.execute(sql_write_to_terms_table,
                                                  (target_id, protein_name, protein_id, ontology if ontology is not None else "", pTerm, score, 0, 1,))
                finally:
                    if db_connection is not None:
                        MySQLConnectionPool.get_instance().close_connection(db_connection, pool_name='CACHE')

            else:
                term_file = self.term_file_dict[term_file_key] if term_file_key in self.term_file_dict else self.aggressive_debug_dict[term_file_key]
                for pred_term in pred_terms:  # Predicted Terms
                    score = self.cath_predictor.calculate_probability(pred_term)
                    if float(score) > 0.0:
                        outstr = f'{target_id},{protein_name},{protein_id},{ontology},'
                        outstr += f"{pred_term},{score},N,N\n"
                        term_file.write(outstr)

                for pTerm in final_pred_terms.difference(pred_terms):
                    score = self.cath_predictor.calculate_probability(pTerm)
                    if float(score) > 0.0:
                        outstr = f'{target_id},{protein_name},{protein_id},{ontology if ontology is not None else ""},'
                        outstr += f"{pTerm},{score},N,Y\n"
                        term_file.write(outstr)

                term_file.flush()

    def write_to_cafa_file(self, cafa_file_key: Tuple, pred_terms: Set, final_pred_terms: Set, target_id: str, scoredResults: List):
        if len(pred_terms) > 0:
            if len(scoredResults) > 0:
                scored_terms = [x['goterm'] for x in scoredResults]
            else:
                scored_terms = list()

            if self.use_database_for_output:
                # fl_name = cafa_file.name
                # tbl_name = os.path.splitext(os.path.basename(fl_name))[0]
                # tbl_name = cafa_file_key[0]

                predictions_schema = os.getenv('PREDICTIONS_SCHEMA')

                sql_write_to_cafa_table = f"""INSERT INTO `{predictions_schema}`.`{self.cafa_group}_{cafa_file_key}_{self.cafa_species}_go` (`CAFA_Target`, `GO_Term`, 
                `Score`) VALUES(%s, %s, %s);"""

                db_connection = None
                try:
                    db_connection = MySQLConnectionPool.get_instance().get_connection(pool_name='CACHE')
                    with db_connection.cursor() as db_cursor:
                        final_pred_terms_calc = set()
                        actual_parent_scores = dict()
                        for pred_term in pred_terms:  # Predicted Terms
                            if pred_term in scored_terms:
                                score = float(scoredResults[scored_terms.index(pred_term)]['score'])
                            else:
                                score = float(self.cath_predictor.calculate_probability(pred_term))
                            if score > 0.0:
                                db_cursor.execute(sql_write_to_cafa_table, (target_id, pred_term, score,))

                            parent_terms = self.go_tool.get_parental_terms_with_level(pred_term)
                            levels = len(set([x for x in parent_terms.values()]))
                            score_diff = (1 - float(score)) / levels
                            final_pred_terms_calc.add(pred_term)
                            for pTerm in parent_terms.keys():
                                tmp_score = score + (score_diff * int(parent_terms[pTerm]))
                                tmp_score = tmp_score if tmp_score < 1.0 and pTerm not in ['GO:0008150', 'GO:0003674', 'GO:0005575'] else 1.0
                                final_pred_terms_calc.add(pTerm)
                                if pTerm in actual_parent_scores.keys():
                                    # actual_parent_scores[pTerm] = max(actual_parent_scores[pTerm], tmp_score)
                                    actual_parent_scores[pTerm] = ((actual_parent_scores[pTerm] + tmp_score) / 2)  # Use average
                                else:
                                    actual_parent_scores[pTerm] = tmp_score

                            for pTerm in actual_parent_scores.keys():
                                db_cursor.execute(sql_write_to_cafa_table, (target_id, pTerm, actual_parent_scores[pTerm],))

                        for pTerm in final_pred_terms.difference(final_pred_terms_calc):
                            if pTerm in scored_terms:
                                score = float(scoredResults[scored_terms.index(pTerm)]['score'])
                            else:
                                score = float(self.cath_predictor.calculate_probability(pTerm))
                            if score > 0.0:
                                db_cursor.execute(sql_write_to_cafa_table, (target_id, pTerm, score,))
                finally:
                    if db_connection is not None:
                        MySQLConnectionPool.get_instance().close_connection(db_connection, pool_name='CACHE')
            else:
                cafa_file = self.cafa_file_dict[cafa_file_key]
                final_pred_terms_calc = set()
                actual_parent_scores = dict()
                to_print_terms = dict()
                for pred_term in pred_terms:  # Predicted Terms
                    # score = self.cath_predictor.calculate_probability(pred_term)
                    if pred_term in scored_terms:
                        score = float(scoredResults[scored_terms.index(pred_term)]['score'])
                    else:
                        score = float(self.cath_predictor.calculate_probability(pred_term))
                    if score > 0.0:
                        if pred_term not in to_print_terms.keys():
                            to_print_terms[pred_term] = score
                        else:
                            # to_print_terms[pred_term] = max(score, to_print_terms[pred_term])  # Using max
                            to_print_terms[pred_term] = ((score + to_print_terms[pred_term]) / 2)  # Using average
                        # outstr = f'{target_id}\t{pred_term}\t{float(score):.2f}\n'
                        # cafa_file.write(outstr)

                    parent_terms = self.go_tool.get_parental_terms_with_level(pred_term)
                    levels = len(set([x for x in parent_terms.values()]))
                    score_diff = (1 - float(score)) / levels if levels > 0 else 0
                    final_pred_terms_calc.add(pred_term)
                    for pTerm in parent_terms.keys():
                        tmp_score = score + (score_diff * int(parent_terms[pTerm]))
                        tmp_score = tmp_score if tmp_score < 1.0 and pTerm not in ['GO:0008150', 'GO:0003674', 'GO:0005575'] else 1.0
                        final_pred_terms_calc.add(pTerm)
                        if pTerm in actual_parent_scores.keys():
                            # actual_parent_scores[pTerm] = max(actual_parent_scores[pTerm], tmp_score)
                            actual_parent_scores[pTerm] = ((actual_parent_scores[pTerm] + tmp_score) / 2)  # Use average
                        else:
                            actual_parent_scores[pTerm] = tmp_score

                    for pTerm in actual_parent_scores.keys():
                        # outstr = f'{target_id}\t{pTerm}\t{actual_parent_scores[pTerm]:.2f}\n'
                        # cafa_file.write(outstr)
                        if pTerm not in to_print_terms.keys():
                            to_print_terms[pTerm] = actual_parent_scores[pTerm]
                        else:
                            # to_print_terms[pTerm] = max(actual_parent_scores[pTerm], to_print_terms[pTerm])
                            to_print_terms[pTerm] = ((actual_parent_scores[pTerm] + to_print_terms[pTerm]) / 2)  # Use average

                for pTerm in final_pred_terms.difference(final_pred_terms_calc):
                    if pTerm in scored_terms:
                        score = float(scoredResults[pTerm.index(pTerm)]['score'])
                    else:
                        score = float(self.cath_predictor.calculate_probability(pTerm))

                    # score = self.cath_predictor.calculate_probability(pTerm)
                    if score > 0.0:
                        # outstr = f'{target_id}\t{pTerm}\t{float(score):.2f}\n'
                        # cafa_file.write(outstr)
                        if pTerm not in to_print_terms.keys():
                            to_print_terms[pTerm] = score
                        else:
                            # to_print_terms[pTerm] = max(score, to_print_terms[pTerm]) # use Max
                            to_print_terms[pTerm] = ((score + to_print_terms[pTerm]) / 2)  # Use average

                for term in to_print_terms.keys():
                    outstr = f'{target_id}\t{term}\t{to_print_terms[term]:.2f}\n'
                    cafa_file.write(outstr)

                cafa_file.flush()

    def write_to_pred_file(self, pred_file_key: Tuple, results: List, target_id: str, protein_name: str, protein_id: str, ontology: str, expected_go_terms: int,
                           predicted_terms: int):

        if predicted_terms > 0:
            if self.use_database_for_output:
                # fl_name = pred_file_funfams.name
                # tbl_name = os.path.splitext(os.path.basename(fl_name))[0]
                # tbl_name = pred_file_key[0]

                predictions_schema = os.getenv('PREDICTIONS_SCHEMA')

                sql_write_to_pred_file_table = f"""INSERT INTO `{predictions_schema}`.`prediction-{pred_file_key}-cafa{self.cafa_edition}_sp_species_{self.cafa_species}-t
                {str(self.threshold).replace('.', '_')}-FF` (
                                    `CAFA_Target`,  
                                    `Protein_Name`,  
                                    `Protein_Unique_Id`,  
                                    `Ontology`,  
                                    `expected_count`,  
                                    `funfams_predicted_count`,  
                                    `funfams_tp`,  
                                    `funfams_fp`,  
                                    `funfams_fn`,  
                                    `funfams_tn`,  
                                    `funfams_precision`,  
                                    `funfams_recall`,  
                                    `funfams_specificity`,  
                                    `funfams_f1_score`,  
                                    `funfams_false_negative_rate`,  
                                    `funfams_false_positive_rate`,  
                                    `funfams_false_discovery_rate`,  
                                    `funfams_fowlkes_mallows_index`,  
                                    `funfams_accuracy`)
                                    VALUES(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s);
                            """

                db_connection = None
                try:
                    db_connection = MySQLConnectionPool.get_instance().get_connection(pool_name='CACHE')
                    with db_connection.cursor() as db_cursor:
                        values = [target_id, protein_name, protein_id, ontology, expected_go_terms, predicted_terms]
                        values += results

                        db_cursor.execute(sql_write_to_pred_file_table, tuple(values))
                finally:
                    if db_connection is not None:
                        MySQLConnectionPool.get_instance().close_connection(db_connection, pool_name='CACHE')
            else:
                self.logger.debug(
                    f"Writing to pred file {pred_file_key}. In pred_file_dict = {pred_file_key in self.pred_file_dict}; In aggressive_debug_pred_dict = "
                    f"{pred_file_key in self.aggressive_debug_pred_dict}")
                pred_file_funfams = self.pred_file_dict[pred_file_key] if pred_file_key in self.pred_file_dict else self.aggressive_debug_pred_dict[pred_file_key]
                if pred_file_funfams is None:
                    logger.error(f"Pred file is None. Pred file dict is {self.pred_file_dict.keys()}; Aggressive Pred Debug Dict is {self.aggressive_debug_pred_dict.keys()}")

                pred_file_funfams.write(f"{target_id},{protein_name},{protein_id},{ontology},{expected_go_terms},{predicted_terms},")
                if results is None or results == []:
                    self.logger.debug("In write to pred file: Results are empty")
                pred_file_funfams.write(",".join(str(v) for v in results))
                pred_file_funfams.write("\n")

                pred_file_funfams.flush()

    @exception(create_logger())
    def predict_a_sequence(self, sequence: str):
        cafa_set = self.cafa_set
        indexes_to_use = self.indexes_to_use
        term_file_dict = self.term_file_dict
        pred_file_dict = self.pred_file_dict
        cafa_file_dict = self.cafa_file_dict
        aggressive_debug = self.aggressive_debug

        ontologies_list = ['BP', 'CC', 'MF']
        searches = Searches()

        local_logger = self.logger
        target_logger = self.target_logger

        local_logger.debug(sequence)
        fasta_handle = io.StringIO(sequence)
        hashed_sequence = self.cath_predictor.get_hashed_sequence(sequence)

        # Step 1A: Identify the sequence
        # for seq_record in SeqIO.parse(fasta_handle, "fasta"):
        target_id, protein_description, searchable_sequence = self.cath_predictor.step0_identify_sequence(sequence)

        # if not self.cath_predictor.target_in_cafa(target_id):
        #     return

        #  Get predictions from funfhmmer
        ff_preds_BP, ff_preds_CC, ff_preds_MF, scoredFunFamerPredictions = self.predict_by_funfhmmer(target_id)

        if (ff_preds_BP is None and ff_preds_CC is None and ff_preds_MF is None) or (len(ff_preds_BP) == 0 and len(ff_preds_CC) == 0 and len(ff_preds_MF) == 0):
            self.notfound_logger.debug(f"Target {target_id} not found by Funfhmmer predictor")

        # Attempt to get data from cache
        protein_id, protein_name, expected_go_terms, bp_terms, cc_terms, mf_terms, superfamilies, funfams = \
            self.cath_predictor.get_data_from_cache(cafa_set=cafa_set, hashed_sequence=hashed_sequence, sequence=sequence, funfhmmer=self.funfhmmer,
                                                    protein_description=protein_description, searchable_sequence=searchable_sequence, search=searches,
                                                    skip_sanity_checks=self.skip_not_in_cache, output=[])

        if target_id is not None and target_id.strip() != '':
            if protein_id is None or protein_name is None:
                if protein_description is not None:
                    protein_id = protein_description.split(' ')[0]
                    protein_name = protein_description.split(' ')[1]
            self.logger.debug(f"Processing {protein_name} (Blasted protein id {protein_id}, Description {protein_description}), target_id {target_id} from cache")
            self.logger.debug(f"Superfamilies: <{superfamilies}>")
            self.logger.debug(f"Funfams: <{funfams}>")
            self.logger.debug(f"Expected Go Terms: <{expected_go_terms}>")
        else:
            if self.skip_not_in_cache:
                self.logger.debug("Data from Cache is not available. Skipping due to skip_not_in_cache setting.")  # return
            else:
                self.logger.debug("Data from Cache is not available. Trying to find the data again.")

                protein_id, protein_name = self.cath_predictor.step1A_protein_name_and_id(searchable_sequence, protein_description)

                expected_go_terms, bp_terms, cc_terms, mf_terms = self.cath_predictor.get_expected_go_terms(protein_id, searchable_sequence, searches)

                funfams = self.cath_predictor.step1B_funfhmmer(self.funfhmmer, protein_description, searchable_sequence, None)[1]

                local_logger.debug("Related families for protein: {}".format(funfams))
                funfams = set()
                if funfams is not None:
                    for family in funfams:
                        # superfamily = family[0:family.find('.FF')]
                        # superfamilies.add(superfamily)
                        funfams.add(family)

        # Sanity Checks
        sanity_check_results = self.cath_predictor.sanity_checks(cafa_set, protein_id, protein_name, protein_description, sequence, searchable_sequence, self.funfhmmer,
                                                                 superfamilies, funfams, expected_go_terms, bp_terms, cc_terms, mf_terms, searches, None,
                                                                 do_not_retry_fetches=self.skip_not_in_cache)
        sanity_check_results_done = True
        if not sanity_check_results[0]:
            self.logger.warn(f"Sanity checks failed. Cannot do predictions for {protein_description}")
            sanity_check_results_done = False  # return  # Sanity checks failed. Exit
        else:
            superfamilies = sanity_check_results[1]
            funfams = sanity_check_results[2]
            # expected_go_terms = set(sanity_check_results[3])
            bp_terms = set(sanity_check_results[4])
            cc_terms = set(sanity_check_results[5])
            mf_terms = set(sanity_check_results[6])

        if superfamilies is None or funfams is None or len(superfamilies) == 0 or len(funfams) == 0:
            self.logger.warn("Cannot do predictions if superfamilies and funfams are None")
            sanity_check_results_done = False  # return  # Cannot do predictions

        resultsList = {}
        all_family_go_terms_ontology = {}
        scoredResults = {}
        scoresResultsList = {}
        for ontology in self.ontologies:
            if funfams is not None and len(funfams) > 0:
                for index in indexes_to_use:
                    if index != 'dcgo':
                        results = self.cath_predictor.predict_from_funfams_set_based('_samefunfam_similarity_allEC', ontology, True, list(funfams), index, enriched=True)
                    else:
                        results = self.cath_predictor.predict_from_funfams_dcgo(ontology, True, list(funfams), enriched=True)

                    all_family_go_terms_ontology[ontology] = results[3]

                    if index in resultsList.keys():
                        res = resultsList[index]
                        scoredResults = scoresResultsList[index]
                        if results[1] is not None:
                            res[ontology] = set(results[1])
                            scoredResults[ontology] = results[0]
                        else:
                            res[ontology] = set()
                            scoredResults[ontology] = list()
                    else:
                        res = dict()
                        scoredResults = dict()
                        if results[1] is not None:
                            res[ontology] = set(results[1])
                            scoredResults[ontology] = results[0]
                        else:
                            res[ontology] = set()
                            scoredResults[ontology] = list()
                        resultsList[index] = res
                        scoresResultsList[index] = scoredResults

                # resultsList[index + '-' + ontology] = results

        local_logger.debug(resultsList)
        local_logger.debug(scoresResultsList)

        # Combine results and return
        bp_predicted_intersection = set()
        bp_predicted_union = set()
        cc_predicted_intersection = set()
        cc_predicted_union = set()
        mf_predicted_intersection = set()
        mf_predicted_union = set()

        for index in indexes_to_use:
            local_logger.debug(f'Currently doing index: {index}')

            results = resultsList[index] if index in resultsList else None
            scoredResults = scoresResultsList[index] if index in scoresResultsList else list()
            ex_go_terms = set()
            ex_go_terms_enriched = set()
            if results is not None:
                for ontology in self.ontologies:
                    if ontology == 'BP':
                        ex_go_terms = bp_terms
                        if len(bp_predicted_intersection) == 0:
                            bp_predicted_intersection = results[ontology]
                        else:
                            if len(results[ontology]) > 0:  # Ensure that there are elements in the set as otherwise you end up with an equally empty set - no predictions
                                bp_predicted_intersection = bp_predicted_intersection.intersection(set(results[ontology]))

                        if len(bp_predicted_union) == 0:
                            bp_predicted_union = results[ontology]
                        else:
                            bp_predicted_union = bp_predicted_union.union(results[ontology])

                    if ontology == 'CC':
                        ex_go_terms = cc_terms
                        if len(cc_predicted_intersection) == 0:
                            cc_predicted_intersection = set(results[ontology])
                        else:
                            if len(results[ontology]) > 0:  # Same reason as BP
                                cc_predicted_intersection = cc_predicted_intersection.intersection(set(results[ontology]))

                        if len(cc_predicted_union) == 0:
                            cc_predicted_union = set(results[ontology])
                        else:
                            cc_predicted_union = cc_predicted_union.union(set(results[ontology]))

                    if ontology == 'MF':
                        ex_go_terms = mf_terms
                        if len(mf_predicted_intersection) == 0:
                            mf_predicted_intersection = set(results[ontology])
                        else:
                            if len(results[ontology]) > 0:  # Same reason as BP
                                mf_predicted_intersection = mf_predicted_intersection.intersection(set(results[ontology]))

                        if len(mf_predicted_union) == 0:
                            mf_predicted_union = set(results[ontology])
                        else:
                            mf_predicted_union = mf_predicted_union.union(set(results[ontology]))

                    # noinspection PyTypeChecker
                    if ex_go_terms is not None or len(ex_go_terms) > 0:
                        ex_go_terms_enriched = self.cath_predictor.enrich_go_terms(set(ex_go_terms))

                    if aggressive_debug and results is not None:
                        # term_file = aggressive_debug_dict[index.lower()]
                        # pred_file = aggressive_debug_pred_dict[index.lower()]
                        # cafa_file = cafa_file_dict[index.lower()]

                        if results[ontology] is not None and len(results[ontology]) > 0:  # Only do if there are predictions!
                            pred_terms_enriched = self.cath_predictor.enrich_go_terms(results[ontology])

                            self.write_to_term_file((index.lower(),), ex_go_terms, ex_go_terms_enriched, results[ontology], pred_terms_enriched, target_id, protein_name,
                                                    protein_id, ontology)

                            self.write_to_cafa_file((index.lower(),), results[ontology], pred_terms_enriched, target_id, scoredResults[ontology])

                            if self.inherit_go_terms:
                                metrics_results = self.cath_predictor.calculate_protein_prediction_metrics(ex_go_terms_enriched, pred_terms_enriched, all_family_go_terms_ontology[
                                    ontology]) if ex_go_terms_enriched is not None else set()

                                self.write_to_pred_file((index.lower(),), metrics_results, target_id, protein_name, protein_id, ontology, len(ex_go_terms_enriched),
                                                        len(pred_terms_enriched))
                            else:
                                metrics_results = self.cath_predictor.calculate_protein_prediction_metrics(ex_go_terms_enriched, results[ontology], all_family_go_terms_ontology[
                                    ontology]) if ex_go_terms_enriched is not None else set()

                                self.write_to_pred_file((index.lower(),), metrics_results, target_id, protein_name, protein_id, ontology, len(ex_go_terms),
                                                        len(results[ontology]))

        local_logger.debug("Intersected only")
        local_logger.debug(f"BP JB_Intersected (only) count: {len(bp_predicted_intersection)}")
        local_logger.debug(f"CC JB_Intersected (only) count: {len(cc_predicted_intersection)}")
        local_logger.debug(f"MF JB_Intersected (only) count: {len(mf_predicted_intersection)}")

        local_logger.debug("Intersected and intersected with funfhmmer")
        # noinspection GrazieInspection
        local_logger.debug(f"BP JB_Intersected intersected with FF count: {len(bp_predicted_intersection.intersection(ff_preds_BP))}")
        # noinspection GrazieInspection
        local_logger.debug(f"CC JB_Intersected intersected with FF count: {len(cc_predicted_intersection.intersection(ff_preds_CC))}")
        # noinspection GrazieInspection
        local_logger.debug(f"MF JB_Intersected intersected with FF count: {len(mf_predicted_intersection.intersection(ff_preds_MF))}")

        local_logger.debug("Intersected and unioned with funfhmmer")
        local_logger.debug(f"BP JB_Intersected unioned with FF count: {len(bp_predicted_intersection.union(ff_preds_BP))}")
        local_logger.debug(f"CC JB_Intersected unioned with FF count: {len(cc_predicted_intersection.union(ff_preds_CC))}")
        local_logger.debug(f"MF JB_Intersected unioned with FF count: {len(mf_predicted_intersection.union(ff_preds_MF))}")

        local_logger.debug("Unioned only")
        local_logger.debug(f"BP JB_Unioned (only) count: {len(bp_predicted_union)}")
        local_logger.debug(f"CC JB_Unioned (only) count: {len(cc_predicted_union)}")
        local_logger.debug(f"MF JB_Unioned (only) count: {len(mf_predicted_union)}")

        local_logger.debug("Unioned and intersected with funfhmmer")
        local_logger.debug(f"BP JB_Unioned intersected with FF count: {len(bp_predicted_union.intersection(ff_preds_BP))}")
        local_logger.debug(f"CC JB_Unioned intersected with FF count: {len(cc_predicted_union.intersection(ff_preds_CC))}")
        local_logger.debug(f"MF JB_Unioned intersected with FF count: {len(mf_predicted_union.intersection(ff_preds_MF))}")

        local_logger.debug("Unioned and unioned with funfhmmer")
        # noinspection GrazieInspection
        local_logger.debug(f"BP JB_Unioned unioned with FF count: {len(bp_predicted_union.union(ff_preds_BP))}")
        # noinspection GrazieInspection
        local_logger.debug(f"CC JB_Unioned unioned with FF count: {len(cc_predicted_union.union(ff_preds_CC))}")
        # noinspection GrazieInspection
        local_logger.debug(f"MF JB_Unioned unioned with FF count: {len(mf_predicted_union.union(ff_preds_MF))}")

        local_logger.debug("Funfhmmer only data")
        local_logger.debug(f"BP FF count: {len(ff_preds_BP)}")
        local_logger.debug(f"CC FF count: {len(ff_preds_CC)}")
        local_logger.debug(f"MF FF count: {len(ff_preds_MF)}")

        if term_file_dict is None or pred_file_dict is None or cafa_file_dict is None:
            raise BaseException("The Term File Dictionary or the Prediction File Dictionary or the Cafa File Dictionary is/are None")

        # Output to the Term Files
        prediction_files = ['jb', 'ff', 'jb-i-ff', 'jb-u-ff']
        variants = ['int', 'un']

        for pf in prediction_files:
            if pf != 'ff':
                for variant in variants:
                    file_index = (f'{pf}-{variant}', 'FF')
                    # term_file = term_file_dict[(f'{pf}-{variant}', 'FF')]
                    # pred_file = pred_file_dict[(f'{pf}-{variant}', 'FF')]
                    # cafa_file = cafa_file_dict[(f'{pf}-{variant}', 'FF')]

                    # 1. Output the expected terms
                    pred_terms = None
                    ex_go_terms = None
                    ex_go_terms_enriched = None
                    for ontology in ontologies_list:
                        if ontology == 'BP':
                            ex_go_terms = set(bp_terms) if bp_terms is not None else set()
                            ex_go_terms_enriched = self.cath_predictor.enrich_go_terms(ex_go_terms) if self.inherit_go_terms else ex_go_terms
                            if variant == 'int':
                                pred_terms = bp_predicted_intersection
                            elif variant == 'un':
                                pred_terms = bp_predicted_union

                            if pf == 'jb-i-ff':
                                pred_terms = pred_terms.intersection(ff_preds_BP)
                            elif pf == 'jb-u-ff':
                                pred_terms = pred_terms.union(ff_preds_BP)

                        elif ontology == 'CC':
                            ex_go_terms = set(cc_terms) if cc_terms is not None else set()
                            ex_go_terms_enriched = self.cath_predictor.enrich_go_terms(ex_go_terms) if self.inherit_go_terms else ex_go_terms
                            if variant == 'int':
                                pred_terms = cc_predicted_intersection
                            elif variant == 'un':
                                pred_terms = cc_predicted_union

                            if pf == 'jb-i-ff':
                                pred_terms = pred_terms.intersection(ff_preds_CC)
                            elif pf == 'jb-u-ff':
                                pred_terms = pred_terms.union(ff_preds_CC)

                        elif ontology == 'MF':
                            ex_go_terms = set(mf_terms) if mf_terms is not None else set()
                            ex_go_terms_enriched = self.cath_predictor.enrich_go_terms(ex_go_terms) if self.inherit_go_terms else ex_go_terms
                            if variant == 'int':
                                pred_terms = mf_predicted_intersection
                            elif variant == 'un':
                                pred_terms = mf_predicted_union

                            if pf == 'jb-i-ff':
                                pred_terms = pred_terms.intersection(ff_preds_MF)
                            elif pf == 'jb-u-ff':
                                pred_terms = pred_terms.union(ff_preds_MF)

                        pred_terms_enriched = self.cath_predictor.enrich_go_terms(pred_terms)

                        if len(pred_terms) == 0:  # If there are no predicted terms for the ontology, do not write anything
                            continue

                        self.write_to_term_file(file_index, ex_go_terms, ex_go_terms_enriched, pred_terms, pred_terms_enriched, target_id, protein_name, protein_id, ontology)

                        scoredResultsToCafaFile = []
                        if 'Sorensen' in scoresResultsList:
                            scoredResultsToCafaFile.extend(scoresResultsList['Sorensen'][ontology])  # This gives the most reliable results, higher than Jaccard, less than Overlap
                        else:
                            scoredResultsToCafaFile.extend(scoredFunFamerPredictions)  # This is when methods do not return predictions, so use Funfamer

                        self.write_to_cafa_file(file_index, pred_terms, pred_terms_enriched, target_id, scoredResultsToCafaFile)

                        if sanity_check_results_done:
                            if self.inherit_go_terms:
                                results = self.cath_predictor.calculate_protein_prediction_metrics(ex_go_terms_enriched, pred_terms_enriched,
                                                                                                   all_family_go_terms_ontology[ontology]) if len(ex_go_terms) > 0 else set()
                                self.write_to_pred_file(file_index, results, target_id, protein_name, protein_id, ontology, len(ex_go_terms_enriched), len(pred_terms_enriched))
                            else:
                                results = self.cath_predictor.calculate_protein_prediction_metrics(ex_go_terms, pred_terms, all_family_go_terms_ontology[ontology]) if len(
                                    ex_go_terms) > 0 else set()

                                self.write_to_pred_file(file_index, results, target_id, protein_name, protein_id, ontology, len(ex_go_terms), len(pred_terms))

            else:
                # print('Outputting Funfhmmer data')
                file_index = ('ff', 'FF')
                # term_file = term_file_dict[('ff', 'FF')]
                # pred_file = pred_file_dict[('ff', 'FF')]
                # cafa_file = cafa_file_dict[('ff', 'FF')]

                for ontology in ontologies_list:
                    pred_terms = set()
                    ex_go_terms = set()
                    ex_go_terms_enriched = set()
                    pred_terms_enriched = set()

                    if ontology == 'BP':
                        ex_go_terms = set(bp_terms) if bp_terms is not None else set()
                        ex_go_terms_enriched = self.cath_predictor.enrich_go_terms(ex_go_terms) if self.inherit_go_terms else ex_go_terms
                        pred_terms = ff_preds_BP
                        pred_terms_enriched = self.cath_predictor.enrich_go_terms(pred_terms) if self.inherit_go_terms else pred_terms

                    elif ontology == 'CC':
                        ex_go_terms = set(cc_terms) if cc_terms is not None else set()
                        ex_go_terms_enriched = self.cath_predictor.enrich_go_terms(ex_go_terms) if self.inherit_go_terms else ex_go_terms
                        pred_terms = ff_preds_CC
                        pred_terms_enriched = self.cath_predictor.enrich_go_terms(pred_terms) if self.inherit_go_terms else pred_terms

                    elif ontology == 'MF':
                        ex_go_terms = set(mf_terms) if mf_terms is not None else set()
                        ex_go_terms_enriched = self.cath_predictor.enrich_go_terms(ex_go_terms) if self.inherit_go_terms else ex_go_terms
                        pred_terms = ff_preds_MF
                        pred_terms_enriched = self.cath_predictor.enrich_go_terms(pred_terms) if self.inherit_go_terms else pred_terms

                    if len(pred_terms) == 0:  # If there are no expected terms for the ontology, do not write anything
                        continue

                    self.write_to_term_file(file_index, ex_go_terms, ex_go_terms_enriched, pred_terms, pred_terms_enriched, target_id, protein_name, protein_id, ontology)

                    self.write_to_cafa_file(file_index, pred_terms, pred_terms_enriched, target_id, list())

                    if sanity_check_results_done and len(ex_go_terms) > 0:
                        results = self.cath_predictor.calculate_protein_prediction_metrics(ex_go_terms_enriched, pred_terms, all_family_go_terms_ontology[ontology])
                        self.write_to_pred_file(file_index, results, target_id, protein_name, protein_id, ontology, len(ex_go_terms), len(pred_terms))

                    target_logger.debug(target_id)

    @exception(create_logger())
    def predict_list_of_sequences(self, sequences_to_annotate: List[str]):
        for sequence in sequences_to_annotate:
            self.predict_a_sequence(sequence)

    def create_tables(self, predictions_schema: str, filename: str):
        sql_create_pred_table = f"""
                                CREATE OR REPLACE TABLE `{predictions_schema}`.`prediction-{filename}-{self.cafa_set}-t{str(self.threshold).replace('.', '_')}-FF` (
                                    `id` MEDIUMINT NOT NULL AUTO_INCREMENT,
                                    `CAFA_Target` VARCHAR(20) NOT NULL,  
                                    `Protein_Name` VARCHAR(15),  
                                    `Protein_Unique_Id` VARCHAR(30),  
                                    `Ontology` VARCHAR(2),  
                                    `expected_count` INT,  
                                    `funfams_predicted_count` INT,  
                                    `funfams_tp` INT,  
                                    `funfams_fp` INT,  
                                    `funfams_fn` INT,  
                                    `funfams_tn` INT,  
                                    `funfams_precision` FLOAT,  
                                    `funfams_recall` FLOAT,  
                                    `funfams_specificity` FLOAT,  
                                    `funfams_f1_score` FLOAT,  
                                    `funfams_false_negative_rate` FLOAT,  
                                    `funfams_false_positive_rate` FLOAT,  
                                    `funfams_false_discovery_rate` FLOAT,  
                                    `funfams_fowlkes_mallows_index` FLOAT,  
                                    `funfams_accuracy` FLOAT,  
                                    KEY `idx_protein_name` (`Protein_Name`) USING BTREE,  
                                    KEY `idx_protein_id` (`Protein_Unique_Id`) USING BTREE,  
                                    KEY `idx_ontology` (`Ontology`) USING BTREE,  
                                        PRIMARY KEY (`id`)
                                    ) ENGINE=InnoDB;
                            """

        sql_create_terms_table = f"""
                                CREATE OR REPLACE TABLE `{predictions_schema}`.`terms-{filename}-{self.cafa_set}-t{str(self.threshold).replace('.', '_')}-FF` (
                                    `id` MEDIUMINT NOT NULL AUTO_INCREMENT,
                                    `CAFA_Target` VARCHAR(20) NOT NULL,
                                    `Protein_Name` VARCHAR(15),
                                    `Protein_Unique_Id` VARCHAR(30),
                                    `Ontology` VARCHAR(2),
                                    `GO_Term` VARCHAR(15),
                                    `Score` FLOAT,
                                    `Is_From_TPS` BIT,
                                    `Is_Inherited_From_Parents` BIT,
                                    KEY `idx_protein_name` (`Protein_Name`) USING BTREE,
                                KEY `idx_protein_id` (`Protein_Unique_Id`) USING BTREE,
                                KEY `idx_ontology` (`Ontology`) USING BTREE,
                                    PRIMARY KEY (`id`)
                                ) ENGINE=InnoDB;
                            """

        sql_create_cafa_table = f"""
                                CREATE OR REPLACE TABLE `{predictions_schema}`.`{self.cafa_group}_{filename}_{self.cafa_species}_go` (
                                    `id` MEDIUMINT NOT NULL AUTO_INCREMENT,
                                    `CAFA_Target` VARCHAR(20) NOT NULL,
                                    `GO_Term` VARCHAR(15),
                                    `Score` FLOAT,
                                    PRIMARY KEY (`id`)
                                ) ENGINE=InnoDB;
                            """

        db_connection = MySQLConnectionPool.get_instance().get_connection(pool_name='CACHE')
        with db_connection.cursor() as db_cursor:
            db_cursor._defer_warnings = True

            db_cursor.execute(sql_create_pred_table)
            db_cursor.execute(sql_create_terms_table)
            db_cursor.execute(sql_create_cafa_table)

    def setup_files(self):
        try:
            prediction_files = ['jb', 'ff', 'jb-i-ff', 'jb-u-ff']

            variants = ['int', 'un']

            if self.use_database_for_output:
                if self.aggressive_debug:
                    prediction_files += ['jaccard', 'overlap', 'sorensen', 'dcgo']

                predictions_schema = os.getenv('PREDICTIONS_SCHEMA')

                for pf in prediction_files:
                    if pf not in ['jaccard', 'overlap', 'sorensen', 'dcgo', 'ff']:
                        for variant in variants:
                            self.create_tables(predictions_schema, f'{pf}-{variant}')
                    elif pf == 'ff':
                        self.create_tables(predictions_schema, 'ff')
                    else:
                        self.create_tables(predictions_schema, pf)
            else:
                import socket
                tmp_path = os.getenv('PREDICTION_SAVE_PATH_' + socket.gethostname().upper())

                cath_version = os.getenv('CATH_VERSION_' + socket.gethostname().upper())
                pathname = ''
                if 'dcgo' in self.indexes_to_use:
                    pathname = os.path.join(tmp_path, "predictions", f"{cath_version}_IndexMethods+dcgo_t{str(self.threshold).replace('.', '_')}")
                else:
                    pathname = os.path.join(tmp_path, "predictions", f"{cath_version}_IndexMethods_t{str(self.threshold).replace('.', '_')}")  # self.ts))

                os.makedirs(pathname, exist_ok=True)

                for pf in prediction_files:
                    if pf not in ['jaccard', 'overlap', 'sorensen', 'dcgo', 'ff']:
                        for variant in variants:
                            pred_file_method_variant = open(os.path.join(pathname,
                                                                         "prediction-{0}-{1}-t{2}-FF.csv".format(f'{pf}-{variant}',
                                                                                                                 self.cafa_set,
                                                                                                                 str(self.threshold).replace('.', '_'))), "w")

                            term_file_method_variant = open(os.path.join(pathname,
                                                                         "terms-{0}-{1}-t{2}-FF.csv".format(f'{pf}-{variant}',
                                                                                                            self.cafa_set,
                                                                                                            str(self.threshold).replace('.', '_'))), "w")

                            cafa_file_method_variant = open(
                                os.path.join(pathname, "{0}_{1}_{2}_go.txt".format(self.cafa_group, f'{pf}-{variant}',
                                                                                   self.cafa_species)), "w")

                            self.pred_file_dict[(f'{pf}-{variant}', 'FF')] = pred_file_method_variant
                            self.term_file_dict[(f'{pf}-{variant}', 'FF')] = term_file_method_variant
                            self.cafa_file_dict[(f'{pf}-{variant}', 'FF')] = cafa_file_method_variant

                            pred_file_method_variant.write("CAFA Target, Protein Name,Protein Unique Id,{}expected_count,funfams_predicted_count,funfams_tp,"
                                                           "funfams_fp,funfams_fn,funfams_tn,funfams_precision,funfams_recall,"
                                                           "funfams_specificity,funfams_f1_score,funfams_false_negative_rate,funfams_false_positive_rate,"
                                                           "funfams_false_discovery_rate, funfams_fowlkes_mallows_index,funfams_accuracy".format(
                                "Ontology, " if self.separate_ontologies else ""))

                            term_file_method_variant.write("CAFA Target, Protein Name,Protein Unique Id,{}GO Term,Score,Is From TPS,"
                                                           "Is Inherited From Parents".format("Ontology, " if self.separate_ontologies else ""))

                            cafa_file_method_variant.write(f"AUTHOR\t{self.cafa_group}\n")
                            cafa_file_method_variant.write(f"MODEL\t1\n")
                            cafa_file_method_variant.write(f"KEYWORDS\tsequence-profile alignment, homolog, hidden Markov model.\n")

                            pred_file_method_variant.write("\n")
                            pred_file_method_variant.flush()

                            term_file_method_variant.write("\n")
                            term_file_method_variant.flush()

                            cafa_file_method_variant.flush()
                    else:
                        pred_file_ff = open(
                            os.path.join(pathname, "prediction-{0}-{1}-t{2}-FF.csv".format('ff', self.cafa_set, str(self.threshold).replace('.', '_'))),
                            "w")

                        term_file_ff = open(
                            os.path.join(pathname, "terms-{0}-{1}-t{2}-FF.csv".format('ff', self.cafa_set, str(self.threshold).replace('.', '_'))), "w")

                        cafa_file_ff = open(os.path.join(pathname, "{0}_ff_{1}_go.txt".format(self.cafa_group, self.cafa_species)), "w")

                        self.pred_file_dict[('ff', 'FF')] = pred_file_ff
                        self.term_file_dict[('ff', 'FF')] = term_file_ff
                        self.cafa_file_dict[('ff', 'FF')] = cafa_file_ff

                        pred_file_ff.write("CAFA Target, Protein Name,Protein Unique Id,{}expected_count,funfams_predicted_count,funfams_tp,"
                                           "funfams_fp,funfams_fn,funfams_tn,funfams_precision,funfams_recall,"
                                           "funfams_specificity,funfams_f1_score,funfams_false_negative_rate,funfams_false_positive_rate,"
                                           "funfams_false_discovery_rate, funfams_fowlkes_mallows_index,funfams_accuracy".format("Ontology, " if self.separate_ontologies else ""))

                        term_file_ff.write("CAFA Target, Protein Name,Protein Unique Id,{}GO Term,Score,Is From TPS,"
                                           "Is Inherited From Parents".format("Ontology, " if self.separate_ontologies else ""))

                        cafa_file_ff.write(f"AUTHOR\t{self.cafa_group}\n")
                        cafa_file_ff.write(f"MODEL\t1\n")
                        cafa_file_ff.write(f"KEYWORDS\tsequence-profile alignment, homolog, hidden Markov model.\n")

                        pred_file_ff.write("\n")
                        pred_file_ff.flush()

                        term_file_ff.write("\n")
                        term_file_ff.flush()

                        cafa_file_ff.flush()

                if self.aggressive_debug:
                    debug_index_files = ['jaccard', 'overlap', 'sorensen', 'dcgo']
                    for f in debug_index_files:
                        term_file_method = open(
                            os.path.join(pathname, "terms-{0}-{1}-t{2}-FF.csv".format(f, self.cafa_set, str(self.threshold).replace('.', '_'))), "w")

                        pred_file_method = open(
                            os.path.join(pathname, "prediction-{0}-{1}-t{2}-FF.csv".format(f, self.cafa_set, str(self.threshold).replace('.', '_'))), "w")

                        cafa_file_method = open(os.path.join(pathname, "{0}_{1}_{2}_go.txt".format(self.cafa_group, f, self.cafa_species)), "w")

                        term_file_method.write("CAFA Target, Protein Name,Protein Unique Id,{}GO Term,Score,Is From TPS,"
                                               "Is Inherited From Parents".format("Ontology, " if self.separate_ontologies else ""))

                        cafa_file_method.write(f"AUTHOR\t{self.cafa_group}\n")
                        cafa_file_method.write(f"MODEL\t1\n")
                        cafa_file_method.write(f"KEYWORDS\tsequence-profile alignment, homolog, hidden Markov model.\n")

                        pred_file_method.write("CAFA Target, Protein Name,Protein Unique Id,{}expected_count,funfams_predicted_count,funfams_tp,"
                                               "funfams_fp,funfams_fn,funfams_tn,funfams_precision,funfams_recall,"
                                               "funfams_specificity,funfams_f1_score,funfams_false_negative_rate,funfams_false_positive_rate,"
                                               "funfams_false_discovery_rate, funfams_fowlkes_mallows_index,funfams_accuracy".format(
                            "Ontology, " if self.separate_ontologies else ""))

                        term_file_method.write("\n")
                        term_file_method.flush()

                        cafa_file_method.flush()

                        ft = (f,)

                        self.cafa_file_dict[ft] = cafa_file_method

                        pred_file_method.write("\n")
                        pred_file_method.flush()

                        self.aggressive_debug_dict[ft] = term_file_method
                        self.aggressive_debug_pred_dict[ft] = pred_file_method
        except Exception as ex:
            raise BaseException(ex)

    def close_files(self):
        if not self.use_database_for_output:
            self.logger.info("Closing all files ...")

            for index in self.pred_file_dict:
                text_file = self.pred_file_dict[index]
                text_file.flush()
                text_file.close()

            for index in self.term_file_dict:
                term_file = self.term_file_dict[index]
                term_file.flush()
                term_file.close()

            for index in self.cafa_file_dict:
                cafa_file = self.cafa_file_dict[index]
                cafa_file.write('END')
                cafa_file.flush()
                cafa_file.close()

            if self.aggressive_debug:
                for index in self.aggressive_debug_dict:
                    debug_file = self.aggressive_debug_dict[index]
                    debug_file.flush()
                    debug_file.close()

                for index in self.aggressive_debug_pred_dict:
                    debug_file = self.aggressive_debug_pred_dict[index]
                    debug_file.flush()
                    debug_file.close()

    def setup_cafa_details(self, file: str):
        cafa_edition = os.getenv('CAFA_EDITION')
        filename = os.path.splitext(file)[0]

        cafa_species = ''
        if os.getenv('PREDICTION_NAMING_SCHEME') == 'CAFA2':
            cafa_set = filename[0:file.find('.')].replace('.', '_').lower()
        else:
            cafa_set = f"cafa{cafa_edition}_{filename.replace('.', '_').lower()}"
            cafa_species = filename.replace('sp_species.', '')

        return cafa_edition, cafa_set, cafa_species

    def cafa_predictor_preprocessing(self, path: str):
        print('Pre-Processing')

        if self.pred_file_dict is None:
            self.pred_file_dict = dict()

        if self.term_file_dict is None:
            self.term_file_dict = dict()

        if self.cafa_file_dict is None:
            self.cafa_file_dict = dict()

        if self.aggressive_debug_dict is None:
            self.aggressive_debug_dict = dict()

        if self.aggressive_debug_pred_dict is None:
            self.aggressive_debug_pred_dict = dict()

        (self.cafa_edition, self.cafa_set, self.cafa_species) = self.setup_cafa_details(path)

        # self.logger.debugsequences)

        self.setup_files()

        cache_connection = MySQLConnectionPool.get_instance().get_connection(pool_name='CACHE')
        self.cath_predictor.create_cache_tables(cache_connection, self.cafa_set)

    def predict_from_path(self, path, parallel_backend='threading', use_threads: bool = False):
        for file in os.listdir(path):
            if not os.path.isdir(os.path.join(path, file)):
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

                self.cafa_predictor_preprocessing(file)

                try:
                    max_threads = int(os.getenv('NUM_THREADS'))

                    sequences = sorted(list(sequences))

                    if use_threads:
                        # Parallel(n_jobs=max_threads, verbose=5, backend=parallel_backend)(
                        #     delayed(unwrap_self_predict)(sequence, indexes_to_use=self.indexes_to_use, cafa_set=self.cafa_set, pred_file_dict=self.pred_file_dict,
                        #                                  term_file_dict=self.term_file_dict, cafa_file_dict=self.cafa_file_dict, aggressive_debug=self.aggressive_debug,
                        #                                  aggressive_debug_dict=self.aggressive_debug_dict, aggressive_debug_pred_dict=self.aggressive_debug_pred_dict) for sequence
                        #     in zip([self] * len(sequences), sequences))

                        import sys
                        sys.setrecursionlimit(60000)

                        with Parallel(backend=parallel_backend, n_jobs=max_threads, verbose=55) as parallel:
                            parallel(delayed(unwrap_self_predict)(sequence) for sequence in zip([self] * len(sequences), sequences))
                    else:
                        self.predict_list_of_sequences(sequences)

                    from utilities.emailmanager import EmailManager
                    EmailManager.send_message('joseph.bonello@um.edu.mt', 'CAFA set finished', f"CAFA set {self.cafa_set} completed")
                except Exception as err:
                    self.logger.error(err)
                    self.logger.error(traceback.format_exc())
                    from utilities.emailmanager import EmailManager
                    EmailManager.send_message('joseph.bonello@um.edu.mt', 'Prediction Error', "\r\n".join(["In predictions_cath.py", traceback.format_exc()]))

                finally:
                    self.close_files()
