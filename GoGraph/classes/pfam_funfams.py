import shutil
import traceback

from utilities.mysqlconnectionpool import MySQLConnectionPool


class PFAMmer:

    def __init__(self):
        self.__version__ = 1.0

        import os
        self.phdPrefix = os.getenv("PHD_DB_PREFIX")

    def meets_inclusion_threshold(self, pfam_fun_fam: str, bitscore: float, verbose: bool = False,
                                  output: list = []):
        check_sql = f"""
            SELECT CASE WHEN %s >= cff.funfam_bit_score_threshold
                           THEN 1 ELSE 0 END AS "include", cff.funfam_bit_score_threshold
            FROM {self.phdPrefix}.GoGraph_pfamfunfamilies cff
            WHERE cff.pfamfunfamilyfull_id = %s
        """

        connection = MySQLConnectionPool.get_instance().get_connection()

        with connection.cursor() as cursor:
            cursor.execute(check_sql, (bitscore, pfam_fun_fam,))

            row = cursor.fetchone()

            if verbose:
                output.append(f"For <{pfam_fun_fam}>, <bitscore: {row['bitscore']}>")

            if row is None:
                include_protein = False
                output.append(f"Funfam {pfam_fun_fam} may be missing, but no value returned for query")
            else:
                include_protein = True if row['include'] == 1 else False

        MySQLConnectionPool.get_instance().close_connection(connection)

        return include_protein

    def pfams_search_offline(self, sequence: str, output: list = [], verbose: bool = False,
                             protein_description: str = 'QUERY'):
        import os
        import tempfile
        import random
        import string

        unique = ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(5))

        basepath = os.path.join(tempfile.gettempdir(), unique)
        try:
            if output is not None:
                output.append("Starting offline scan")

            os.mkdir(basepath)
            input_path = os.path.join(basepath, "seq.fasta")
            domtblout_path = os.path(basepath, "seq.domtblout")
            crh_path = os.path.join(basepath, 'seq.crh')

            fasta = open(input_path, "w")
            if not ">" in sequence:
                fasta.write(">{}\n".format(protein_description))
            fasta.write(sequence)

            fasta.flush()
            fasta.close()

            matches = set()

            # ./apps/pfam-genomescan.pl -i data/test.fasta -l data/funfam-hmm3-v4_1_0.lib -o results/
            import socket
            pfam_files_data_path = os.getenv("PFAM_FILES_" + socket.gethostname().upper())

            pfam_tools_process = 'hmmsearch'
            cath_resolve_hits_path = os.getenv("CATH_RESOLVE_HITS_" + socket.gethostname().upper())
            # pfam_scan.pl -fasta ./test.fasta -dir ../pfamfiles/ -outfile ./output.json

            import subprocess
            subprocess.check_output([pfam_tools_process, '--cut_tc', '--cpu', '5', '-o', '/dev/null', '--domtblout', domtblout_path, pfam_files_data_path, input_path, ],
                shell=False, env=os.environ)

            subprocess.check_output([cath_resolve_hits_path, '--input-format', 'hmmer_domtblout', '--hits-text-to-file', crh_path, domtblout_path, ], shell=False, env=os.environ)

            pfam_output = open(crh_path, "r")

            import re
            regexp = re.compile(" +")

            for line in pfam_output:
                if not line.startswith("#") and not line == '\n':
                    fields = line.split(" ")

                    funfam = fields[1].replace('FF/', 'FF').replace('/', '.')
                    bitscore = float(fields[2])
                    evalue = float(fields[6])

                    include_protein = self.meets_inclusion_threshold(funfam, bitscore, evalue)

                    if verbose:
                        output.append(f"^<Match id: {funfam}>, <Bit Score: {bitscore}>, ")
                        output.append(f"<E-Value: {evalue}>, <include: {include_protein}>^")

                    if include_protein:
                        matches.add(funfam)

            if verbose:
                output.append("Offline scan completed with {} results.".format(len(matches)))
            return matches
        except OSError as err:
            print("OS error: {0}".format(err))
            print("Faulty sequence is: {}".format(sequence))
        except ValueError:
            print("A value error occurred. Check the following log.")
            print("Faulty sequence is: {}".format(sequence))

            from utilities.emailmanager import EmailManager
            EmailManager.send_message('joseph.bonello@um.edu.mt', 'Value Error',
                                      "\r\n".join(["In funfhmmer.py", traceback.format_exc(),
                                                   "", "Current sequence is", sequence]))
        except:
            import sys
            print("Unexpected error:", sys.exc_info()[0])
            print("Faulty sequence is: {}".format(sequence))

            from utilities.emailmanager import EmailManager
            EmailManager.send_message('joseph.bonello@um.edu.mt', 'Unexpected Error',
                                      "\r\n".join(["In pfam_funfams.py", traceback.format_exc(),
                                                   "", "Current sequence is", sequence]))
        finally:
            if os.path.exists(basepath):
                shutil.rmtree(basepath, ignore_errors=True)

        return None

    def pfams_search(self, protein_description, sequence: str, output: list, verbose=True, timeout=36):
        results = None

        if results is None or len(results) == 0:
            results = self.pfams_search_offline(protein_description, sequence, output)

        return results
