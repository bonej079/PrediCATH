from collections import defaultdict

import csv
import matplotlib.pyplot as plt
import os
import pandas as pd
import pickle as pkl
import re
from goatools.obo_parser import GODag


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

    def __init__(self, output: list, data_path: str):
        print("Prediction Analysis CAFA version 1.0 (based on code by Dr Jon Lees)")

        self.output = output
        self.data_path = data_path

        # score cut-offs you want to test for
        self.cutoffs = [0.06, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

        self.basedir = "./"

        self.root_go_terms = {"GO:0008150", "GO:0005575", "GO:0003674"}

        # Load pre-calculated Information contents of GO terms from pickled object
        file_location = os.path.join(self.basedir, "gene_ontology_edit_20160601.obo")
        self.ontology_graph = GODag(obo_file=file_location)

        file_location = os.path.join(self.basedir, "icd.pkl")
        self.information_contents = pkl.load(open(file_location, 'rb'))

        plt.switch_backend('agg')

    def load_csv(self, path: str):
        file = os.path.basename(path)
        match_result = re.match("^terms-([A-Za-z-]+)-[A-Za-z0-9_]+_([a-z0-9]+)-t([0-9_]+).*", file)

        index = match_result.group(1)
        pred_set = match_result.group(2)
        threshold = match_result.group(3).replace("_", ".")

        self.output.append("Loading {0:s} into data frame".format(file))

        df = pd.read_csv(path, header=0, skip_blank_lines=True, error_bad_lines=False, warn_bad_lines=True, engine='c', low_memory=True, memory_map=False)
        df.rename(columns=lambda x: x.strip().replace(" ", "_").lower(), inplace=True)
        self.output.append("Loaded {0:d} records from file".format(df.size))

        df = df.dropna()
        df = df.drop_duplicates()
        self.output.append("After cleansing, {0:d} records remain".format(df.size))

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
            if score == 0:
                continue

            self.preds[row.ontology][sequence_id][go_term] = 0.5 if score < 0.0 else score

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

        self.output.append(f"Finished loading {count} predictions for this set")
        self.output.append("")

    # True Postitive Set
    def load_tps(self, data: pd.DataFrame, output_to_disk: bool = False, path_to_output: str = None):
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

        self.output.append(f"Finished loading {count} true postitives for this set")
        self.output.append("")

    def run_information_theoretic_benchmark(self, analysis_path: str, plot_path: str, index: str, cafa_set: str, ontology: str, do_plots=False):
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

    def semantic_distance_euclid(self, perfs: list, weighted: bool):
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

    def getFmax(self, perfs):
        result = []
        for perf in perfs:
            precision = perf["precision"]
            recall = perf["recall"]

            fmax = 0 if (recall + precision) == 0 else (2 * precision * recall) / (recall + precision)

            result.append([fmax, [recall, precision, perf["cutoff"]]])

        result.sort(reverse=True)
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
