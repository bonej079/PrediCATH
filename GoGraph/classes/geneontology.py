from collections.abc import Iterable

from goatools.obo_parser import GODag


# To get obo from fit lfs do:
# git lfs pull gene_ontology_edit_20160601.obo


class GeneOntologyUtility:
    __ontology_graph__ = None
    __ontology_graph_path__ = None

    def __init__(self, obo_path: str = None, file_name: str = None):
        import os

        import sys
        path_to_script = os.path.dirname(os.path.abspath(sys.modules[GeneOntologyUtility.__module__].__file__))

        if file_name is None:
            file_name = "gene_ontology_edit_20160601.obo"
            # file_name = 'go-basic.obo'

        if obo_path is None:
            # if getattr(sys, 'gettrace', None) is not None:
            #     print(path_to_script)
            file_location = os.path.join(path_to_script,
                                         file_name)  # file_location = os.path.join(path_to_script, "go_20130615-termdb.obo")  # file_location =  # os.path.join(path_to_script,
            # "go_20191007.obo")  # CAFA 4 OBO
        else:
            file_location = obo_path

        # file_location = os.path.join('.', "go-basic.obo")
        if GeneOntologyUtility.__ontology_graph__ is None and GeneOntologyUtility.__ontology_graph_path__ != file_location:
            GeneOntologyUtility.__ontology_graph__ = GODag(obo_file=file_location)
            GeneOntologyUtility.__ontology_graph_path__ = file_location

        self.ontology_graph = GeneOntologyUtility.__ontology_graph__

    def query_term(self, term: str):
        go_query_result = self.ontology_graph.query_term(term, verbose=False)

        return go_query_result

    def get_ontology_for_term(self, term: str):
        term = term.strip()
        go_result = str(self.ontology_graph.query_term(term, verbose=False))

        return go_result[go_result.index('[') + 1:go_result.index(']')] if go_result.index('[') >= 0 else ''

    def get_ontology_for_term_list(self, term_list: Iterable):
        term_list = [x.strip() for x in term_list]
        results = dict()

        if term_list is not None:
            for term in term_list:
                go_result = str(self.ontology_graph.query_term(term, verbose=False))
                if go_result.find('[') >= 0:
                    results[term] = go_result[go_result.index('[') + 1:go_result.index(']')].lower()

        return results

    def get_parental_terms(self, term: str):
        term = term.strip()
        go_query_result = self.ontology_graph.query_term(term, verbose=False)
        all_parents = set() if go_query_result is None else go_query_result.get_all_parents()

        return all_parents

    def get_parental_terms_with_level(self, term: str):
        term = term.strip()
        all_parents = dict()
        go_query_result = self.ontology_graph.query_term(term, verbose=False)
        all_parents_list = set() if go_query_result is None else go_query_result.get_all_parents()
        for parent in all_parents_list:
            if parent not in all_parents.keys():
                go_query_result = self.ontology_graph.query_term(parent, verbose=False)
                all_parents[go_query_result.item_id] = go_query_result.level

        return all_parents

    def get_direct_parents(self, term: str):
        term = term.strip()
        go_query_result = self.ontology_graph.query_term(term, verbose=False).parents

        return go_query_result

    def get_direct_parents_as_set(self, term: str):
        return set(res.id for res in self.get_direct_parents(term))

    def get_parental_terms_for_list(self, terms: list):
        terms = [x.strip() for x in terms]
        query_terms = [self.ontology_graph.query_term(term, verbose=False) for term in terms if term is not None]
        go_query_result = [term.get_all_parents() for term in query_terms if term is not None]  # Fixes issue with missing terms

        import itertools
        return list(itertools.chain.from_iterable(go_query_result))

    def create_info_content(self):
        # TODO: https://github.com/tanghaibao/goatools/blob/master/notebooks/dcnt_and_tinfo.ipynb

        # To generate the gaf file (in terminal window):
        # Step 1: Get all unique evidence codes from gpad
        # cat goa_uniprot_all_noiea.gpad | cut -d$'\t' -f6 | sort | uniq
        # Step 2: Get the gpa from goa_uniprot_all.gpa
        # cat goa_uniprot_all.gpa | grep -E 'ECO:0000245|ECO:0000247|ECO:0000250|ECO:0000255|ECO:0000266|ECO:0000269|ECO:0000270|ECO:0000303|ECO:0000304|ECO:0000305|ECO:0000307
        # |ECO:0000314|ECO:0000315|ECO:0000316|ECO:0000317|ECO:0000320|ECO:0000353|ECO:0007005}' > goa_uniprot_all_no_iea.gpa

        with open('/tmp/infocontent.txt', 'w') as f:
            import csv
            writer = csv.writer(f)

            from goatools.base import get_godag
            godag = get_godag("gene_ontology_edit_20160601.obo")

            # go_ids_bp = set(o.item_id for o in self.ontology_graph.values() if not o.is_obsolete and o.namespace == 'biological_process')
            # go_ids_cc = set(o.item_id for o in self.ontology_graph.values() if not o.is_obsolete and o.namespace == 'cellular_component')
            # go_ids_mf = set(o.item_id for o in self.ontology_graph.values() if not o.is_obsolete and o.namespace == 'molecular_function')
            all_go_ids = set(o.item_id for o in self.ontology_graph.values() if not o.is_obsolete)

            from goatools.gosubdag.gosubdag import GoSubDag
            gosubdag = GoSubDag(all_go_ids, godag)

            # Read Annotations
            from goatools.anno.factory import get_objanno
            annoobj = get_objanno("goa_uniprot_all_no_iea.gpa", 'gpad', godag=godag)

            # Get associations (Gene ID - to - set of GO IDs) for all namespaces (BP, MF, and CC)
            id2goids = annoobj.get_id2gos_nss()

            # Create TermCounts object
            from goatools.semantic import TermCounts
            tcntobj = TermCounts(godag, id2goids)

            # get the deepest GO Id
            # go_id, go_term = max(gosubdag.go2obj.items(), key=lambda t: t[1].depth)
            # print(go_id, go_term.name)

            for term in all_go_ids:  # go_ids_bp:
                considererd_term = self.query_term(term)
                considered_terms = considererd_term.get_all_parents()
                considered_terms.add(considererd_term.item_id)

                print('{N} ancestors for {GO} "{name}"'.format(N=len(considered_terms), GO=considererd_term.item_id, name=considererd_term.name))

                nts = [gosubdag.go2nt[go] for go in considered_terms if go in gosubdag.go2nt]

                fmt_str = '{I:2}) {NS} {GO:10} {dcnt:11}        D{depth:02}  {GO_name}'

                # Print selected GO information
                print('IDX NS GO ID      Descendants Count Depth Name')
                print('--- -- ---------- ----------------- ----- --------------------')
                for idx, nt_go in enumerate(sorted(nts, key=lambda nt: nt.depth), 1):
                    print(fmt_str.format(I=idx, **nt_go._asdict()))

                # Manage a subset of GO IDs using GoSubDag
                # This time, load the GO term counts
                gosubdag = GoSubDag(considered_terms, godag, tcntobj=tcntobj)

                fmt_str = '{NS} {GO:10} {dcnt:5}  {tinfo:6.3f}  D{depth:02}  {GO_name}'
                print('NS GO         dcnt   tinfo depth name')
                print('-- ---------- ----- ------  ---  -----------------------')
                nts = gosubdag.prt_goids(considered_terms, prtfmt=fmt_str)
                for ntsgo in nts:
                    lst = [ntsgo[11], ntsgo[10], ntsgo[9], ntsgo[5], ntsgo[0]]
                    writer.writerow(lst)


if __name__ == '__main__':
    gou = GeneOntologyUtility()
    gou.create_info_content()
