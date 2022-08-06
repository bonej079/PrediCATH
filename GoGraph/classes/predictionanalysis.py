import os
import pandas as pd
import re


class PredictionAnalysis:
    output = None
    datasets = {}
    sets = None
    ontologies = None

    def __init__(self, output):
        print("Prediction Analysis version 1.0")
        self.output = output

    def load_csv(self, path):
        file = os.path.basename(path)
        index = re.match("^(.*)-merge", file).group(1)

        self.output.append("Loading {0:s} into dataframe".format(file))

        df = pd.read_csv(path, header=0)
        df.rename(columns=lambda x: x.strip().replace(" ", "_").lower(), inplace=True)
        self.output.append("Loaded {0:d} records from file".format(df.size))

        df = df.dropna()
        self.output.append("After cleansing, {0:d} records remain".format(df.size))

        self.datasets[index] = df

    def fill_meta_data(self):
        df = self.datasets[next(x for x in self.datasets)]

        self.sets = df.set.drop_duplicates().tolist()

        self.ontologies = df.ontology.drop_duplicates().tolist()

        self.output.append("Found sets: {0}".format(self.sets))
        self.output.append("Found ontologies: {0}".format(self.ontologies))

        self.output.append("")

    def data_index_set_ontology(self):
        data_for_index_set_ontology = {}

        for index in self.datasets:
            for set in self.sets:
                for ontology in self.ontologies:
                    df = self.datasets[index]
                    df = df[df.set == set]
                    df = df[df.ontology == ontology]

                    data_for_index_set_ontology[(index, set, ontology)] = df

        return data_for_index_set_ontology

