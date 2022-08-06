from collections import defaultdict


class Benchmark:
    def __init__(self, ancestor_path, benchmark_path):
        """
        Initialize the benchmark.

        Input:
        benchmark_path is ontology specific
        ancestor_path is ontology specific
        """

        # Key: protein
        # Value: set of benchmark leaf terms
        self.ancestors = defaultdict(set)
        # Read GO ancestors file generated with go_ontology_ancestors_split_write()
        # File format:
        # go_term <tab> ancestor_1,ancestor_2,..,ancestor_n
        with open(ancestor_path) as ancestors_input:
            for inline in ancestors_input:
                inrec = inline.strip().split('\t')
                term = inrec[0]
                if len(inrec) == 1:
                    self.ancestors[term] = set({})
                else:
                    term_ancestors = inrec[1]
                    self.ancestors[term] = set(term_ancestors.split(','))

        self.true_base_terms = defaultdict(set)
        with open(benchmark_path) as benchmark_input:
            for inline in benchmark_input:
                protein, term = inline.strip().split('\t')
                self.true_base_terms[protein].add(term)

    # def propagate(self):
    #     """
    #     Progate Benchmark terms.
    #     """
    #
    #     # Key: protein
    #     # Value: set of benchmark propagated terms
    #     root_terms = ['GO:0008150', 'GO:0005575', 'GO:0003674']
    #     self.true_terms = defaultdict(set)
    #     for protein in self.true_base_terms:
    #         for term in self.true_base_terms[protein]:
    #             try:
    #
    #                 ancestors = self.ancestors[term].difference(root_terms)  # delete root term in self.true_terms
    #             # modified on 20170203
    #             except KeyError:
    #                 print(f"{term} not found")
    #             self.true_terms[protein].add(term)
    #             self.true_terms[protein] |= ancestors

    def propagate(self):
        """
        Progate prediction terms.
        """

        # Key: protein
        # Value: set of benchmark propagated terms
        root_terms = ['GO:0008150', 'GO:0005575', 'GO:0003674']
        self.true_terms = defaultdict(set)
        from GoGraph.classes.geneontology import GeneOntologyUtility
        go_tool = GeneOntologyUtility()

        for protein in self.true_base_terms:
            for term in self.true_base_terms[protein]:
                ancestors = set(go_tool.get_parental_terms(term))
                ancestors = ancestors.difference(root_terms)
                self.true_terms[protein].add(term)
                self.true_terms[protein] |= ancestors
