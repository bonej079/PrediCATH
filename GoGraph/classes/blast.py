import re
import traceback

from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline, NcbipsiblastCommandline
from bioservices import NCBIblast
from lxml import etree


# Xml Tutorial: http://www.diveintopython3.net/xml.html

class Blast:
    # Correct values are['1e-200', '1e-100', '1e-50', '1e-10', '1e-5', '1e-4', '1e-3', '1e-2', '1e-1', '1.0', '10 (default)', '100', '1000']
    E_VALUE_THRESH = '1e-100'  # 1e-10'  # 0.04

    def search_by_blast_bioservices(self, sequence, verbose=False):
        try:
            s = NCBIblast(verbose)
            s.TIMEOUT = 60
            jobid = s.run(program="blastp", sequence=sequence, email="joseph.bonello@um.edu.mt",
                          stype="protein",
                          database=["uniprotkb_swissprot"])  # "uniprotkb", "pdb" , exp=self.E_VALUE_THRESH)

            # s.getResult(jobid, "out")
            # retTypes = s.getResultTypes(jobid)

            from time import sleep
            count = 0
            # Fix for ValueError: job is not finished
            # str(jobid) fixes issue with Blast unable to concat str and an int
            while s.get_status(str(jobid)) != "FINISHED" and count < (s.TIMEOUT / 5):
                count += 1
                sleep(5)  # Time in seconds.

            ret_xml = str(s.get_result(jobid, 'xml'))

            ret_xml = re.sub(" encoding='.*'", "", ret_xml)

            NSMAP = {'ebi': 'http://www.ebi.ac.uk/schema'}

            from io import StringIO
            tree = etree.parse(StringIO(ret_xml))
            # root = tree.getroot()

            # sequencesimilaritysearchresult = tree.find('{http://www.ebi.ac.uk/schema}sequencesimilaritysearchresult')
            # hits = sequencesimilaritysearchresult.find('{http://www.ebi.ac.uk/schema}hits')
            # firstHit = hits[0]

            entries = tree.xpath("//ebi:ebiapplicationresult/ebi:sequencesimilaritysearchresult/ebi:hits/ebi:hit",
                                 namespaces=NSMAP)
            if (len(entries) == 0):
                print("Sequence Return entries of length 0")
                return None, None
            else:
                firstHit = entries[0]
                return firstHit.attrib['ac'], firstHit.attrib['id']
        except Exception as err:
            print(err)
            print(traceback.format_exc())
            from utilities.emailmanager import EmailManager
            EmailManager.send_message('joseph.bonello@um.edu.mt', 'Prediction Error',
                                      "\r\n".join(["In blast.py", traceback.format_exc(),
                                                   "", "Current sequence is", sequence]))
            return None, None

    def search_by_blast_offline(self, sequence: str, verbose=False, algorithm='blastp'):
        # from Bio.Blast import NCBIWWW
        # result_handle = NCBIWWW.qblast("blastp", "swissprot", sequence)
        # query=sequence

        import tempfile
        import os
        import random
        import string
        unique = ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(5))

        output_path = os.path.join(tempfile.gettempdir(), f"blast_result-{unique}.xml")
        input_path = os.path.join(tempfile.gettempdir(), f"seq-{unique}.fasta")

        try:
            with open(input_path, "w") as fasta:
                fasta.write(sequence)
                fasta.flush()

            import socket
            if not os.getenv("BLASTDB", False):
                # print("BLASTDB not found. Setting manual entry.")
                os.putenv("BLASTDB", os.getenv("BLASTDB_" + socket.gethostname().upper()))
                # os.environ["BLASTDB"] = "/home/joseph/blastdb"

            cline = None
            if algorithm == 'blastp':
                cline = NcbiblastpCommandline(query=input_path, db="swissprot", evalue=self.E_VALUE_THRESH,
                                          remote=False, out=output_path, outfmt="5")
            if algorithm == 'psiblast':
                cline = NcbipsiblastCommandline(query=input_path, db="swissprot", evalue=self.E_VALUE_THRESH,
                                              remote=False, out=output_path, outfmt="5")

            stdout, stderr = cline()

            with open(output_path) as result_handle:
                blast_records = NCBIXML.parse(result_handle)

                blast_record = next(blast_records)
                swissprot_id = None
                swissprot_name = None

                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        if hsp.match == sequence:
                            parsed_header = alignment.title.split(sep='|')
                            # print(parsed_header)
                            sp_position = parsed_header.index('sp')
                            swissprot_id = parsed_header[sp_position + 1]
                            if swissprot_id.find('.') > -1:
                                swissprot_id = swissprot_id[0:swissprot_id.find('.')]
                            swissprot_name = parsed_header[sp_position + 2]
                            swissprot_name = swissprot_name[0:swissprot_name.find(' ')].strip()
                            break

                if swissprot_id is None and swissprot_name is None:
                    print("swissprot_id and swissprot_name are None (Offline BLAST).")
                    return None, None

            return swissprot_id, swissprot_name

        except Exception as err:
            print(err)
            print(traceback.format_exc())
            from utilities.emailmanager import EmailManager
            EmailManager.send_message('joseph.bonello@um.edu.mt', 'Prediction Error',
                                      "\r\n".join(["In blast.py", traceback.format_exc(),
                                                   "", "Current sequence is", sequence, output_path]))
            return None, None
        finally:
            try:
                os.remove(input_path)
                os.remove(output_path)
            except Exception:
                print("Removing temp files failed in blast.py.")

    def search_by_blast(self, sequence: str, verbose=False, offlineAlgorithm = 'blastp'):
        swissprot_id = None
        swissprot_name = None

        try:

            swissprot_id, swissprot_name = self.search_by_blast_offline(sequence, verbose, algorithm=offlineAlgorithm)

            if swissprot_id is None and swissprot_name is None:  # If offline attempt fails, try online
                if verbose:
                    print("Offine search failed. Trying online.")
                swissprot_id, swissprot_name = self.search_by_blast_bioservices(sequence, verbose)

            if swissprot_name is None and swissprot_id is None:
                return None

        except Exception as err:
            print(err)
            print(traceback.format_exc())

        if swissprot_name is None and swissprot_id is None:
            return None
        else:
            return swissprot_id, swissprot_name
