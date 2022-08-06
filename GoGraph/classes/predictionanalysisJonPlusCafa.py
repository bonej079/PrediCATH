import csv
import math
import os
import pickle as pkl
import re
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from goatools.obo_parser import GODag

from GoGraph.classes.benchmark import Benchmark
from GoGraph.evaluation.Ontology.IO import OboIO


class PredictionAnalysisCafa:
    information_contents = None
    ontology_graph = None

    datasets = dict()
    true_positive_set = dict()
    preds = defaultdict(lambda: defaultdict(dict))

    sets = set()
    ontologies = set()
    root_go_terms = set()

    cutoffs = list()
    output = list()

    basedir = ''
    data_path = ''
    output_path = ''
    currently_processing_index = ''
    currently_processing_set = ''
    currently_processing_ontology = ''

    def __init__(self, output: list, data_path: str, benchmark_directory: str):
        print("Prediction Analysis CAFA version 2.1 (based on code by Dr Jon Lees and CAFA Matlab Code (for reproducibility)")

        self.output = output
        self.data_path = data_path
        self.benchmark_directory = benchmark_directory

        # score cut-offs you want to test for
        self.cutoffs = ['%.2f' % elem for elem in np.arange(0.0, 1.01, 0.01)]  # [0.06, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        self.cutoffs = [float(i) for i in self.cutoffs]

        self.basedir = "./"

        self.root_go_terms = {"GO:0008150", "GO:0005575", "GO:0003674"}

        # Load pre-calculated Information contents of GO terms from pickled object
        # self.obo_path = os.path.join(self.basedir, "gene_ontology_edit_20160601.obo")  # CAFA 3!
        self.obo_path = os.path.join(self.basedir, "go_20191007.obo")  # CAFA 4!
        # self.obo_path = os.path.join(self.basedir, "go_20130615-termdb.obo")  # CAFA 2
        self.ontology_graph = GODag(obo_file=self.obo_path)

        file_location = os.path.join(self.basedir, "icd.pkl")
        self.information_contents = pkl.load(open(file_location, 'rb'))

        plt.switch_backend('agg')

        self.benchmark = dict()
        self.benchmark_set = set()

        self.benchmark_available_taxa_BP = ['ARATH', 'DANRE', 'DICDI', 'DROME', 'ECOLI', 'HUMAN', 'MOUSE', 'PSEAE', 'RAT', 'SCHPO', 'YEAST']
        self.benchmark_available_taxa_CC = ['ARATH', 'DROME', 'ECOLI', 'HUMAN', 'MOUSE', 'RAT', 'YEAST']
        self.benchmark_available_taxa_MF = ['ARATH', 'DROME', 'ECOLI', 'HUMAN', 'MOUSE', 'PSEAE', 'RAT', 'SCHPO', 'YEAST']

    def load_csv(self, path: str):
        file = os.path.basename(path)
        match_result = re.match("^terms-([A-Za-z-]+)-[A-Za-z0-9_]+_([a-z0-9]+)-t([0-9_]+).*", file)

        index = match_result.group(1)
        pred_set = match_result.group(2)
        threshold = match_result.group(3).replace("_", ".")

        self.output.append(f"Loading {file} into data frame")

        df = pd.read_csv(path, header=0, skip_blank_lines=True, error_bad_lines=False, warn_bad_lines=True, engine='c', low_memory=True, memory_map=False)
        df.rename(columns=lambda x: x.strip().replace(" ", "_").lower(), inplace=True)
        self.output.append(f"Loaded {df.size} records from file")

        df = df.dropna()
        df = df.drop_duplicates()
        self.output.append(f"After cleansing, {df.size} records remain")

        df = df.sort_values('cafa_target')  # 'protein_name')  # Sort by protein name just in case the names get confused by threading

        self.datasets[(index, pred_set)] = df
        self.sets.add(pred_set)

    def fill_meta_data(self):
        self.output.append("Processing meta-data")
        df = pd.DataFrame()
        iter = (x for x in self.datasets)
        while df.empty:
            df = self.datasets[next(iter)]

        for ont in df.ontology.drop_duplicates():
            self.ontologies.add(ont)

        self.output.append("Found ontologies: {0}".format(self.ontologies))
        self.output.append("")

    def data_index_set_ontology(self):
        self.output.append("Splitting the datasets by index, set and ontology for predictions")
        data_for_index_set_ontology = {}

        for index in self.datasets:
            df = self.datasets[index]
            for ontology in self.ontologies:
                dfont = df[df.ontology == ontology]

                # index[0] = Scoring Method
                # index[1] = CAFA Set
                data_for_index_set_ontology[(index[0], index[1], ontology)] = dfont

        self.output.append('')
        return data_for_index_set_ontology

    def load_predictions(self, data: pd.DataFrame, output_to_disk: bool = False, path_to_output: str = None):
        self.output.append("Loading the predictions from dataset")

        predicted_data = data[(data.is_from_tps == 'N')]  # & (data.is_inherited_from_parents == 'N')]

        count = 0
        for (index, row) in predicted_data.iterrows():
            sequence_id = row.cafa_target  # row.protein_name
            go_term = row.go_term
            score = float(row.score)

            # Ignore record if score is 0
            if score <= 0:
                continue

            self.preds[row.ontology][sequence_id][go_term] = 0.0 if score < 0.0 else score

            count += 1

        if output_to_disk:
            print("Outputting the predictions file for verification")

            if path_to_output is None or path_to_output.strip() == '':
                raise Exception("If output_to_disk is specified, then you must provide an output path for "
                                "the predictions file")

            file_location = os.path.join(path_to_output, f"{self.currently_processing_index}_{self.currently_processing_set}_"
                                                         f"{self.currently_processing_ontology}.pred_overide")
            predictions = csv.writer(open(file_location, "w"), dialect='excel', delimiter='\t', lineterminator='\n', quoting=csv.QUOTE_NONE)
            predictions.writerow([f"AUTHOR    FUNFAM-JB-{self.currently_processing_index}"])

            for ontology in self.preds:
                for sequence in self.preds[ontology]:
                    for term in self.preds[ontology][sequence]:
                        predictions.writerow([sequence, term, self.preds[ontology][sequence][term]])

        self.output.append(f"Finished loading {count} predictions for set {self.currently_processing_set} and {self.currently_processing_index}")
        self.output.append("")

    def go_ontology_split(self, ontology):
        """
        Split an GO obo file into three ontologies
        by Dr. Friedberg
        """
        MF_terms = set({})
        BP_terms = set({})
        CC_terms = set({})
        for node in ontology.get_ids():  # loop over node IDs and alt_id's
            if ontology.namespace[node] == "molecular_function":
                MF_terms.add(node)
            elif ontology.namespace[node] == "biological_process":
                BP_terms.add(node)
            elif ontology.namespace[node] == "cellular_component":
                CC_terms.add(node)
            else:
                raise (ValueError, "%s has no namespace" % node)
        return MF_terms, BP_terms, CC_terms

    def go_ontology_ancestors_split_write(self, obo_path):
        """
        Input: an OBO file
        Output: 3 files with ancestors
        by Dr. Friedberg
        Updated 20190906 by Ashley: ony write ancestor within the same ontology
        """

        obo_BP_out = open("%s_ancestors_BP.txt" % (os.path.splitext(obo_path)[0]), "w")
        obo_CC_out = open("%s_ancestors_CC.txt" % (os.path.splitext(obo_path)[0]), "w")
        obo_MF_out = open("%s_ancestors_MF.txt" % (os.path.splitext(obo_path)[0]), "w")
        print(os.path.abspath(obo_path))
        obo_parser = OboIO.OboReader(open(os.path.abspath(obo_path)))
        go = obo_parser.read()
        MF_terms, BP_terms, CC_terms = self.go_ontology_split(go)
        for term in MF_terms:
            ancestors = go.get_ancestors(term)
            ancestors_filter = set()
            if len(ancestors) > 0:
                for anc in ancestors:
                    anc_ont = go.get_namespace(anc)
                    if anc_ont == 'molecular_function':
                        ancestors_filter.add(anc)
            obo_MF_out.write("%s\t%s\n" % (term, ",".join(ancestors_filter)))
        for term in BP_terms:
            ancestors = go.get_ancestors(term)
            ancestors_filter = set()
            if len(ancestors) > 0:
                for anc in ancestors:
                    anc_ont = go.get_namespace(anc)
                    if anc_ont == 'biological_process':
                        ancestors_filter.add(anc)
            obo_BP_out.write("%s\t%s\n" % (term, ",".join(ancestors_filter)))
        for term in CC_terms:
            ancestors = go.get_ancestors(term)
            ancestors_filter = set()
            if len(ancestors) > 0:
                for anc in ancestors:
                    anc_ont = go.get_namespace(anc)
                    if anc_ont == 'cellular_component':
                        ancestors_filter.add(anc)
            obo_CC_out.write("%s\t%s\n" % (term, ",".join(ancestors_filter)))

        obo_MF_out.close()
        obo_BP_out.close()
        obo_CC_out.close()
        return [len(BP_terms), len(CC_terms), len(MF_terms)]

    def read_benchmark(self, namespace, species, types, fullbenchmarkfolder, obopath):
        """
        Read Benchmark.

        Input:
        namespace
        species
        types
        fullbenchmarkfolder
        obopath

        Output:
        bench
        """
        # Ancestor files here are precomputed
        # To get the ancestor files, use preprocess.py (go_ontology_ancestors_split_write) DOES NOT EXIST
        legal_types = ["type1", "type2", "typex"]
        legal_subtypes = ["easy", "hard"]
        legal_namespace = ["BP", "MF", "CC", "hpo"]
        legal_species = ['all', 'all-archaea', 'all-bacteria', 'all-eukarya', 'ARATH', 'BACSU', 'CANAX', 'DANRE', 'DICDI', 'DROME', 'ECOLI', 'eukarya', 'HALS3', 'HALVD', 'HELPY',
                         'HUMAN', 'IGNH4', 'METJA', 'MOUSE', 'MYCGE', 'NITMS', 'prokarya', 'PSEAE', 'PSEPK', 'PSESM', 'PYRFU', 'RAT', 'SALCH', 'SALTY', 'SCHPO', 'STRPN', 'SULSO',
                         'XENLA', 'YEAST', ]
        # fullbenchmarkfolder = './precrec/benchmark/'
        if namespace not in legal_namespace:
            print("Namespace not accepted, choose from 'BP', 'CC', 'MF' and 'hpo'\n")
        elif (species not in legal_species) and (species not in legal_subtypes):
            print('Species not accepted')
        elif types not in legal_types:
            print('Type not accepted, choose from "type1","type2" and "typex"\n')

        converted_namespace = 'bpo' if namespace == 'BP' else ('mfo' if namespace == 'MF' else 'cco')
        matchname = converted_namespace + '_' + species + '_' + types + '.txt'
        # generate ancestor files
        # obocounts = self.go_ontology_ancestors_split_write(obopath)
        # obocountDict = {'BP': obocounts[0], 'CC': obocounts[1], 'MF': obocounts[2]}
        # ontology-specific calculations
        if namespace == 'BP':
            full_benchmark_path = os.path.abspath(fullbenchmarkfolder + '/groundtruth/' + 'leafonly_BPO.txt')
            ancestor_path = os.path.splitext(obopath)[0] + "_ancestors_BP.txt"
        elif namespace == 'CC':
            full_benchmark_path = os.path.abspath(fullbenchmarkfolder + '/groundtruth/' + 'leafonly_CCO.txt')
            ancestor_path = os.path.splitext(obopath)[0] + "_ancestors_CC.txt"
        elif namespace == 'MF':
            full_benchmark_path = os.path.abspath(fullbenchmarkfolder + '/groundtruth/' + 'leafonly_MFO.txt')
            ancestor_path = os.path.splitext(obopath)[0] + "_ancestors_MF.txt"

        benchmarkListPath = os.path.abspath(fullbenchmarkfolder + '/lists/' + matchname)
        if not os.path.exists(benchmarkListPath):
            benchmarkListPath = os.path.abspath(fullbenchmarkfolder + '/lists/too_few/' + matchname)

        print(f"Now loading benchmark from {benchmarkListPath}")

        if os.path.isfile(benchmarkListPath) and os.path.getsize(benchmarkListPath) > 0:
            handle = open(benchmarkListPath, 'r')
            prots = set()
            for line in handle:
                prots.add(line.strip())
            handle.close()
            tempfilename = 'temp_%s_%s_%s.txt' % (namespace, species, types)
            tempfile = open(os.path.abspath(fullbenchmarkfolder + '/' + tempfilename), 'w')
            for line in open(os.path.abspath(full_benchmark_path), 'r'):
                prot = line.split('\t')[0]
                if prot in prots:
                    tempfile.write(line)
            tempfile.close()

            if not os.path.exists(ancestor_path):
                self.go_ontology_ancestors_split_write(obopath)

            bench = Benchmark(ancestor_path, tempfile.name)
            bench.propagate()
            os.remove(tempfile.name)
        else:
            print('Benchmark set is empty.\n')
            bench = None
        return bench

    def taxon_name_converter(self, taxonID):
        # convert from taxonomy ID to name (i.e. from 9606 to HUMAN）
        taxonTable = {'10090': 'MOUSE', '10116': 'RAT', '160488': 'PSEPK', '170187': 'STRPN', '186497': 'PYRFU', '208964': 'PSEAE', '223283': 'PSESM', '224308': 'BACSU',
                      '237561': 'CANAX', '243232': 'METJA', '243273': 'MYCGE', '273057': 'SULSO', '284812': 'SCHPO', '309800': 'HALVD', '321314': 'SALCH', '3702': 'ARATH',
                      '436308': 'NITMS', '44689': 'DICDI', '453591': 'IGNH4', '478009': 'HALS3', '559292': 'YEAST', '7227': 'DROME', '7955': 'DANRE', '83333': 'ECOLI',
                      '8355': 'XENLA', '85962': 'HELPY', '9606': 'HUMAN', '99287': 'SALTY', 'all': 'all', 'all-archaea': 'all-archaea', 'all-bacteria': 'all-bacteria',
                      'all-eukarya': 'all-eukarya', 'eukarya': 'eukarya', 'prokarya': 'prokarya', }
        return None if taxonID not in taxonTable else taxonTable[taxonID]

    def typeConverter(self, oldType):
        if oldType == 'type1':
            new_type = 'NK'
        elif oldType == 'type2':
            new_type = 'LK'
        elif oldType == 'all':
            new_type = 'All'
        return new_type

    # True Postitive Set
    def load_tps_all(self, data: pd.DataFrame, output_to_disk: bool = False, path_to_output: str = None):
        tps_data = data[(data.is_from_tps == 'Y')]  # & (data.is_inherited_from_parents == 'N')]

        if output_to_disk:
            print("Outputting the TPS file for verification")

            if path_to_output is None or path_to_output.strip() == '':
                raise Exception("If output_to_disk is specified, then you must provide an output path for the TPS file")

            file_location = os.path.join(path_to_output, f"{self.currently_processing_index}_{self.currently_processing_set}_"
                                                         f"{self.currently_processing_ontology}.valid")
            tps_file = csv.writer(open(file_location, "w"), dialect='excel', delimiter='\t', lineterminator='\n', quoting=csv.QUOTE_NONE)

        count = 0
        for (index, row) in tps_data.iterrows():
            sequence_id = row.cafa_target  # row.protein_name
            go_term = row.go_term

            self.true_positive_set[row.ontology][sequence_id].add(go_term)

            count += 1

            tps_file.writerow([sequence_id, go_term])

        self.output.append(f"Finished loading {count} true positives for this set")
        self.output.append("")

    # True Postitive Set
    def load_tps(self, data: pd.DataFrame, output_to_disk: bool = False, path_to_output: str = None):
        # tps_data = data[(data.is_from_tps == 'Y')]  # & (data.is_inherited_from_parents == 'N')]
        for onto in ['BP', 'CC', 'MF']:
            print(f'ontology: {onto}\n')
            for Type in ['type1', 'type2']:
                print(f'benchmark type:{self.typeConverter(Type)}\n')
                if self.taxon_name_converter(self.currently_processing_set) is None:
                    benchmark = None
                    self.benchmark[f"{self.currently_processing_set}_{onto}_{Type}"] = benchmark
                else:
                    benchmark = self.read_benchmark(onto, self.taxon_name_converter(self.currently_processing_set), Type, self.benchmark_directory, self.obo_path)
                    self.benchmark_set.add(self.currently_processing_set)
                    self.benchmark[f"{self.currently_processing_set}_{onto}_{Type}"] = benchmark

                # self.output.append(f"Finished loading true positives for set {self.currently_processing_set}, ontology: {onto}, type: {Type}")

        self.output.append("True Positive Sets added to Benchmark")

    def propagate(self, predictions: set):
        """
        Progate prediction terms.
        """

        root_terms = ['GO:0008150', 'GO:0005575', 'GO:0003674']
        pred_terms = set()
        from GoGraph.classes.geneontology import GeneOntologyUtility
        go_tool = GeneOntologyUtility(obo_path=self.obo_path)

        for term in predictions:
            ancestors = set(go_tool.get_parental_terms(term))
            ancestors = ancestors.difference(root_terms)  # delete root term in self.true_terms

            pred_terms.add(term)
            pred_terms |= ancestors

        return pred_terms

    def run_information_theoretic_benchmark(self, analysis_path: str, plot_path: str, index: str, cafa_set: str, ontology: str, do_plots=False):
        self.output.append(f"Run Information Theoretical Benchmark for index <{index}> with cafa_set <{cafa_set}> and ontology <{ontology}> with <{len(self.preds[ontology])}> "
                           f"predictions")

        for Type in ['type1', 'type2']:
            file_location_det = os.path.join(analysis_path, f"detailed_perf_cafa_{Type}_{index}_{cafa_set}_{ontology}.csv")
            file_location_sum = os.path.join(analysis_path, f"summary_perf_cafa_{Type}_{index}_{cafa_set}_{ontology}.csv")
            perfs = []

            print(f'benchmark type:{self.typeConverter(Type)}\n')
            benchmark = self.benchmark[f"{self.currently_processing_set}_{ontology}_{Type}"] if f"{self.currently_processing_set}_{ontology}_{Type}" in self.benchmark else None

            if benchmark is None:
                continue

            with open(file_location_det, "w") as detailed:
                with open(file_location_sum, "w") as summary:
                    detailed_stats = csv.writer(detailed, dialect='excel', lineterminator='\n', quoting=csv.QUOTE_ALL)

                    summary_ofh = csv.writer(summary, dialect='excel', lineterminator='\n', quoting=csv.QUOTE_ALL)
                    summary_ofh.writerow(
                        ["branch", "valid_file", "fmax", "recall", "precision", "fmax_coff", "weighted_fmax", "weighted_recall", "weighted_precision", "weighted_fmax_coff", "smin",
                         "ru", "mi", "sd_coff", "semantic distance", "ru", "mi", "sd_coff", "weighted semantic distance", "wru", "wmi", "wsd_coff", "pred_name", "coverage"])

                    sequence_ids = set(benchmark.true_terms.keys())

                    for cutoff in self.cutoffs:
                        # RUMI
                        misinformation, remaining_uncertainty, norm_misinformation, norm_remaining_uncertainty = [], [], [], []
                        # FMAX
                        total_number_of_proteins, total_predicted_terms, proteins_with_predictions, proteins_with_annotation = 0, 0, 0, 0
                        precision, recall, weighted_precision, weighted_recall = [], [], [], []

                        for sequence_id in sequence_ids:

                            total_number_of_proteins += 1

                            true_positive_set = benchmark.true_terms[sequence_id]
                            true_positive_set -= self.root_go_terms

                            if len(true_positive_set) == 0:
                                continue
                            else:
                                proteins_with_annotation += 1  # These are proteins which contribute

                            predictions = set(term for term in self.preds[ontology][sequence_id].keys() if self.preds[ontology][sequence_id][term] <= cutoff)
                            if len(predictions) > 0:
                                predictions = self.propagate(predictions)
                                proteins_with_predictions += 1

                                # Precision + Recall
                                tp_pred_intersection_terms = true_positive_set.intersection(predictions)  # True Positives
                                intersection_size = len(tp_pred_intersection_terms)
                                true_positive_set_size = len(true_positive_set)
                                predicted_set_size = len(predictions)

                                ic_true_positive_set = math.fsum([self.information_contents.get(term, 0) for term in true_positive_set])  # IA for true graph
                                ic_predicted_set = math.fsum([self.information_contents.get(term, 0) for term in predictions])  # IA for predicted graph
                                ic_tp_pred_intersection = math.fsum([self.information_contents.get(term, 0) for term in tp_pred_intersection_terms])  # IA for intersection

                                precision.append(0 if predicted_set_size == 0 else intersection_size / predicted_set_size)
                                recall.append(0 if true_positive_set_size == 0 else intersection_size / true_positive_set_size)

                                weighted_precision.append(0 if ic_predicted_set == 0 else ic_tp_pred_intersection / ic_predicted_set)
                                weighted_recall.append(0 if ic_true_positive_set == 0 else ic_tp_pred_intersection / ic_true_positive_set)
                                total_predicted_terms += len(predictions)

                                # Remaining Uncertainty + Misinformation
                                true_positives_less_predicted = true_positive_set.difference(predictions)  # False Negatives
                                predicted_less_true_positives = predictions.difference(true_positive_set)  # False Positives
                                all_terms = tp_pred_intersection_terms.union(true_positives_less_predicted).union(predicted_less_true_positives)

                                ru = math.fsum([self.information_contents.get(term, 0) for term in true_positives_less_predicted])  # compute ru (FN)
                                mi = math.fsum([self.information_contents.get(term, 0) for term in predicted_less_true_positives])  # compute mi (FP)
                                normaliser = math.fsum([self.information_contents.get(term, 0) for term in all_terms])  # compute normaliser
                                misinformation.append(mi)
                                remaining_uncertainty.append(ru)
                                norm_misinformation.append(mi / normaliser)  # weighted
                                norm_remaining_uncertainty.append(ru / normaliser)  # weighted

                        if proteins_with_predictions == 0:
                            continue

                        pf = dict()
                        pf["cutoff"] = cutoff

                        pf["precision"] = sum(precision) / proteins_with_predictions
                        pf["recall"] = sum(recall) / proteins_with_annotation  # total_number_of_proteins
                        pf["fmax"] = 0 if (pf["precision"] + pf["recall"]) == 0 else 2 * ((pf["precision"] * pf["recall"]) / (pf["precision"] + pf["recall"]))

                        pf["weighted_precision"] = sum(weighted_precision) / proteins_with_annotation  # total_number_of_proteins
                        pf["weighted_recall"] = sum(weighted_recall) / proteins_with_annotation  # total_number_of_proteins
                        pf["weighted_fmax"] = 0 if (pf["weighted_precision"] + pf["weighted_recall"]) == 0 else 2 * (
                                (pf["weighted_precision"] * pf["weighted_recall"]) / (pf["weighted_precision"] + pf["weighted_recall"]))

                        pf["remaining_uncertainty"] = sum(remaining_uncertainty) / proteins_with_annotation  # total_number_of_proteins
                        pf["misinformation"] = sum(misinformation) / proteins_with_annotation  # total_number_of_proteins
                        pf["remaining_uncertainty_weighted"] = sum(norm_remaining_uncertainty) / proteins_with_annotation  # total_number_of_proteins  # weighted
                        pf["misinformation_weighted"] = sum(norm_misinformation) / proteins_with_annotation  # total_number_of_proteins  # weighted

                        perfs.append(pf)

                    if len(perfs) == 0:
                        continue

                    # For the particular ontology under consideration ...
                    weighted_semantic_distance, weighted_cutoff = self.semantic_distance_euclid(perfs, weighted=True)
                    semantic_distance, cutoff = self.semantic_distance_euclid(perfs, weighted=False)
                    fmax, precision_recall = self.getFmax(perfs, False)
                    smin, ru_mi = self.getSmin(perfs, False)
                    weighted_fmax, weighted_precision_recall = self.getFmax(perfs, True)
                    weighted_smin, weighted_ru_mi = self.getSmin(perfs, True)
                    coverage = proteins_with_predictions / proteins_with_annotation

                    detailed_stats.writerow([ontology, "weighted semantic distance", weighted_semantic_distance] + weighted_cutoff)
                    detailed_stats.writerow([ontology, "semantic distance", semantic_distance] + cutoff)
                    detailed_stats.writerow([ontology, "FMAX", fmax] + precision_recall)
                    detailed_stats.writerow([ontology, "SMIN", smin] + ru_mi)
                    detailed_stats.writerow([ontology, "Weighted FMAX", weighted_fmax] + weighted_precision_recall)
                    detailed_stats.writerow([ontology, "Weighted SMIN", weighted_smin] + weighted_ru_mi)

                    summary_ofh.writerow(
                        [ontology, cafa_set] + [fmax] + precision_recall + [weighted_fmax] + weighted_precision_recall + [smin] + ru_mi + [semantic_distance] + cutoff + [
                            weighted_semantic_distance] + weighted_cutoff + [cafa_set] + [coverage])

                    detailed_stats.writerow([ontology, "cut-off", "ru", "mi", "ru-weighted", "mi-weighted", "precision", "recall"])

                    for pf in perfs:
                        detailed_stats.writerow(
                            [ontology, pf["cutoff"], pf["remaining_uncertainty"], pf["misinformation"], pf["remaining_uncertainty_weighted"], pf["misinformation_weighted"],
                             pf["precision"], pf["recall"]])

            print(f"Finished calculation for type {Type}")

    def run_information_theoretic_benchmark_all(self, analysis_path: str, plot_path: str, index: str, cafa_set: str, ontology: str, do_plots=False):
        self.output.append("Run Information Theoretical Benchmark")

        file_location = os.path.join(analysis_path, f"detailed_perf_{index}_{cafa_set}_{ontology}.csv")
        detailed_stats = csv.writer(open(file_location, "w"), dialect='excel', lineterminator='\n', quoting=csv.QUOTE_ALL)

        file_location = os.path.join(analysis_path, f"summary_perf_{index}_{cafa_set}_{ontology}.csv")
        summary_ofh = csv.writer(open(file_location, "w"), dialect='excel', lineterminator='\n', quoting=csv.QUOTE_ALL)
        summary_ofh.writerow(
            ["branch", "valid_file", "fmax", "recall", "precision", "fmax_coff", "semantic distance", "ru", "mi", "sd_coff", "weighted semantic distance", "wru", "wmi", "wsd_coff",
             "pred_name"])

        if ontology not in self.true_positive_set:
            return  # Do not bother if there are no GO terms in the TPS

        perfs = []

        for cutoff in self.cutoffs:
            precision_scores = []
            recall_scores = []
            remaining_uncertainties = []
            misinformations = []
            weighted_remaining_uncertainties = []
            weighted_misinformations = []
            information_contents = []

            proteins_with_predictions = 0.0  # proteins with at least one prediction above cutoff
            total_proteins_with_annotation = 0.0

            for sequence_id, terms in self.true_positive_set[ontology].items():
                terms -= self.root_go_terms

                if len(terms) == 0:
                    continue

                predictions = self.preds[ontology][sequence_id]

                if len(predictions) == 0:  # If there are no predictions, ignore!
                    continue

                total_expected_terms_information_content = self.sum_information_content(terms)

                predicted_go_terms = set()
                predictions_found = False

                for term, score in predictions.items():
                    go_query_result = self.ontology_graph.query_term(term, verbose=False)
                    if go_query_result is None:
                        continue

                    if term in self.root_go_terms:
                        continue

                    if float(score) >= cutoff:
                        predicted_go_terms.add(term)
                        predictions_found = True

                if predictions_found:
                    proteins_with_predictions += 1.0
                    precision = len((terms & predicted_go_terms)) / float(len(predicted_go_terms))
                    precision_scores.append(precision)

                total_proteins_with_annotation += 1.0
                information_contents.append(total_expected_terms_information_content)

                # record information content for True Positives that were not predicted
                remaining_uncertainty = self.sum_information_content(terms - predicted_go_terms)
                remaining_uncertainties.append(remaining_uncertainty)

                # record information content for things which should have been predicted
                misinformation = self.sum_information_content(predicted_go_terms - terms)
                misinformations.append(misinformation)

                recall = len((terms & predicted_go_terms)) / float(len(terms))
                recall_scores.append(recall)

                weighted_remaining_uncertainties.append(remaining_uncertainty * total_expected_terms_information_content)
                weighted_misinformations.append((misinformation * total_expected_terms_information_content))

            if proteins_with_predictions == 0:
                continue

            pf = {"cutoff": cutoff}

            pf["precision"] = sum(precision_scores) / proteins_with_predictions
            pf["recall"] = sum(recall_scores) / total_proteins_with_annotation

            pf["remaining_uncertainty"] = sum(remaining_uncertainties) / total_proteins_with_annotation
            pf["misinformation"] = sum(misinformations) / total_proteins_with_annotation

            total_information_content = sum(information_contents)
            pf["remaining_uncertainty_weighted"] = sum(weighted_remaining_uncertainties) / total_information_content
            pf["misinformation_weighted"] = sum(weighted_misinformations) / total_information_content

            perfs.append(pf)

        if len(perfs) == 0:
            return

        # For the particular ontology under consideration ...
        weighted_semantic_distance, weighted_pt = self.semantic_distance_euclid(perfs, weighted=True)
        semantic_distance, pt = self.semantic_distance_euclid(perfs, weighted=False)
        fmax, precision_recall = self.getFmax(perfs)

        detailed_stats.writerow([ontology, "weighted semantic distance", weighted_semantic_distance] + weighted_pt)
        detailed_stats.writerow([ontology, "semantic distance", semantic_distance] + pt)
        detailed_stats.writerow([ontology, "FMAX", fmax] + precision_recall)

        summary_ofh.writerow([ontology, cafa_set] + [fmax] + precision_recall + [semantic_distance] + pt + [weighted_semantic_distance] + weighted_pt + [cafa_set])

        detailed_stats.writerow([ontology, "cut-off", "ru", "mi", "ru-weighted", "mi-weighted", "precision", "recall"])

        for pf in perfs:
            detailed_stats.writerow(
                [ontology, pf["cutoff"], pf["remaining_uncertainty"], pf["misinformation"], pf["remaining_uncertainty_weighted"], pf["misinformation_weighted"], pf["precision"],
                 pf["recall"]])

        if do_plots:
            x = []
            y = []
            coffs = set()
            plt.rcParams.update(plt.rcParamsDefault)

            for pf in perfs:
                x.append(pf["remaining_uncertainty_weighted"])
                y.append(pf["misinformation_weighted"])
                coffs.add(pf["cutoff"])

            plt.xlim([0, x[-1]])
            plt.ylim(0, 50)
            # plt.plot(weighted_pt[0], weighted_pt[1], "ro")
            plt.plot(x, y)
            file_location = os.path.join(plot_path, f"weighted_rumi_{index}_{cafa_set}_{ontology}.png")
            plt.savefig(file_location, format="png")

            x = []
            y = []
            plt.rcParams.update(plt.rcParamsDefault)

            for pf in perfs:
                x.append(pf["remaining_uncertainty"])
                y.append(pf["misinformation"])
                coffs.add(pf["cutoff"])

            plt.xlim([0, x[-1]])
            plt.plot(pt[0], pt[1], "ro")
            plt.plot(x, y)
            file_location = os.path.join(plot_path, f"rumi_{index}_{cafa_set}_{ontology}.png")
            plt.savefig(file_location, format="png")

            x = []
            y = []
            plt.rcParams.update(plt.rcParamsDefault)

            for pf in perfs:
                x.append(pf["recall"])
                y.append(pf["precision"])
                coffs.add(pf["cutoff"])

            plt.plot(precision_recall[0], precision_recall[1], "ro")
            plt.plot(x, y)
            file_location = os.path.join(plot_path, f"pr_{index}_{cafa_set}_{ontology}.png")
            plt.savefig(file_location, format="png")

    def semantic_distance_euclid(self, perfs: list, weighted: bool):  # Calculates Smin
        """k=1 manhabolis; k=2 euclidean"""
        result = []
        remaining_uncertainties = []
        misinformations = []

        for perf in perfs:
            misinformation = perf["misinformation"]
            remaining_uncertainty = perf["remaining_uncertainty"]

            if weighted:
                misinformation = perf["misinformation_weighted"]
                remaining_uncertainty = perf["remaining_uncertainty_weighted"]

            import math
            sdist = math.sqrt((misinformation * misinformation) + (remaining_uncertainty * remaining_uncertainty))

            result.append([sdist, [remaining_uncertainty, misinformation, perf["cutoff"]]])

            remaining_uncertainties.append(remaining_uncertainty)
            misinformations.append(misinformation)

        result.sort()

        return result[0]

    def getFmax(self, perfs, weighted=False):
        result = []
        for perf in perfs:
            if weighted:
                precision = perf["weighted_precision"]
                recall = perf["weighted_recall"]
            else:
                precision = perf["precision"]
                recall = perf["recall"]

            fmax = 0 if (recall + precision) == 0 else (2 * precision * recall) / (recall + precision)

            result.append([fmax, [recall, precision, perf["cutoff"]]])

        result.sort(reverse=True)
        return result[0]

    def getSmin(self, perfs, weighted=False):
        result = []
        for perf in perfs:
            if weighted:
                remaining_uncertainty = perf["remaining_uncertainty_weighted"]
                misinformation = perf["misinformation_weighted"]
            else:
                remaining_uncertainty = perf["remaining_uncertainty"]
                misinformation = perf["misinformation"]

            smin = math.sqrt((math.pow(remaining_uncertainty, 2) + math.pow(misinformation, 2)))

            result.append([smin, [remaining_uncertainty, misinformation, perf["cutoff"]]])

        result.sort(reverse=False)
        return result[0]

    def reset(self):
        self.true_positive_set = dict()
        self.preds = defaultdict(lambda: defaultdict(dict))

        self.currently_processing_index = ''
        self.currently_processing_set = ''
        self.currently_processing_ontology = ''

        for ontology in self.ontologies:
            self.true_positive_set[ontology] = defaultdict(set)

    def sum_information_content(self, terms):
        total = 0.0
        for go_term in terms:
            total += self.information_contents.get(go_term, 0)
        return total
