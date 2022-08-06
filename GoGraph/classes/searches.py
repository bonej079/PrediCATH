import traceback

from bioservices import UniProt, QuickGO


class Searches:
    """
    Searches for protein GO Terms using UniProt services
    """
    verbose = False

    def uniprot_go_terms_by_protein_id(self, protein_id: str):
        """
        Gets UniProt GO Terms from protein_id
        Args:
            protein_id (str): The protein ID
        """
        u = UniProt(verbose=self.verbose, cache=True)
        results = u.search(protein_id, columns='go-id, ec, reviewed')

        results = results.replace("Gene ontology IDs\tEC number\tStatus\n", "")

        returnVals = set()

        for result in results.split('\n'):
            returnVals.add(result.split('\t')[0])

        return returnVals

    def uniprot_sequence_by_protein_id(self, protein_id: str):
        """
        Args:
            protein_id (str):
        """
        u = UniProt(verbose=self.verbose, cache=True)
        sequence = u.retrieve(protein_id, "fasta")

        return sequence

    def gene_access_from_protein_id(self, protein_id: str):
        """
        Args:
            protein_id (str):
        """
        u = UniProt(verbose=self.verbose)
        results = u.search(protein_id, columns='id')

        results = results.replace("Entry\n", "")

        return results.split("\n")[0]

    def quickgo_terms_by_protein_id(self, protein_id, only_admit_evidence_codes=[]):
        """
        Args:
            protein_id:
            only_admit_evidence_codes:
        """
        try:
            s = QuickGO(verbose=self.verbose)
            # s.url = 'http://www.ebi.ac.uk/QuickGO-Old'
            s.settings.TIMEOUT = 300

            if protein_id.find('.') > 1:
                protein_id = protein_id[0: protein_id.find('.')]

            reply = s.Annotation(geneProductId=protein_id,  # geneProductType='protein',
                                 includeFields="goName,name,synonyms")

            returnVals = set() if (reply is not None and type(reply) == dict and len(reply['results']) > 0) else None

            if returnVals is not None:
                for result in reply['results']:
                    if len(only_admit_evidence_codes) > 0:
                        if result['goEvidence'] in only_admit_evidence_codes:
                            if result['goEvidence'] == 'IEA':
                                for xrefs in result['withFrom']:
                                    for xref in xrefs:
                                        db = xrefs[xref][0]['db']
                                        if (db == 'EC' and result['reference'] == 'GO_REF:0000003') or \
                                                (db == 'UniProtKB-KW' and result['reference'] == 'GO_REF:0000037'):
                                            returnVals.add(result['goId'])
                                        # GO_REF:0000037 FOR EC2GO
                                        # GO_REF:0000003 FOR Keyword2GO
                            else:
                                returnVals.add(result['goId'])
                    else:
                        returnVals.add((result['goId']))

            # results = results.replace('GO ID\tEvidence\tWith\tReference\n', '')

            # if results != '':
            #     for line in results.strip('\n').split('\n'):
            #         split_line = line.split('\t')
            #         if len(only_admit_exp_evidence_codes) > 0:
            #             if split_line[1] in only_admit_exp_evidence_codes:
            #                 if split_line[1] == 'IEA':
            #                     if ('EC' in split_line[2]
            #                         and 'GO_REF:0000003' in split_line[3]) or \
            #                             ('UniProtKB-KW' in split_line[2]
            #                              and 'GO_REF:0000037' in split_line[3]):
            # GO_REF:0000037 FOR EC2GO
            # GO_REF:0000003 FOR Keyword2GO
            #                         returnVals.add(split_line[0])
            #                 else:
            #                     returnVals.add(split_line[0])
            #         else:
            #             returnVals.add(split_line[0])
            # else:
            #     returnVals = None

            return returnVals
        except:
            import sys
            print("Unexpected error:", sys.exc_info()[0])
            print("Faulty protein_id is: {}".format(protein_id))

            from utilities.emailmanager import EmailManager
            EmailManager.send_message('joseph.bonello@um.edu.mt', 'Unexpected Error',
                                      "\r\n".join(["In searches.py", traceback.format_exc(),
                                                   "", "Current protein_id is", protein_id]))

        return None


def old(search: Searches):
    # T01049 (2NPD_BACSU)
    """
    Args:
        search (Searches):
    """
    print("Search for {} returns: {}".format('T01049', search.gene_access_from_protein_id('T01049')))
    print("Search for {} returns: {}".format('2NPD_BACSU', search.gene_access_from_protein_id('T01049')))

    print("Search for {} returns: {}".format('T01049 (2NPD_BACSU)',
                                             search.gene_access_from_protein_id('T01049 (2NPD_BACSU)')))

    print("Protein details for {}: {}\n\n".format('p30443.1', search.uniprot_go_terms_by_protein_id('p30443.1')))

    print("Protein details for {}: {}\n".format('p30443.1', search.quickgo_terms_by_protein_id(protein_id='P30443',
                                                                                               only_admit_evidence_codes=
                                                                                               ['EXP', 'IDA', 'IPI',
                                                                                                'IMP', 'IGI', 'IEP',
                                                                                                'IEA'])))

    # O64545  # print("Protein details for {}: {}\n".format(protein_id, search.quickgo_terms_by_protein_id(protein_id=protein_id,  #
    # only_admit_exp_evidence_codes=[  #                                                                                                'EXP', 'IDA', 'IPI',
    #                                                                                                'IMP', 'IGI', 'TAS'
    #                                                                                                'IEP', 'IC','IEA'])))

    #
    # print("Search for {} returns: {}".format('O05413', search.uniprot_go_terms_by_protein_id('O05413')))
