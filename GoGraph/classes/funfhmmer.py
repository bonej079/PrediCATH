import json
import shutil
import time
import traceback

import requests

from utilities.mysqlconnectionpool import MySQLConnectionPool


class Funfhmmer:
    def __init__(self):
        import os
        self.url_list = os.getenv("FUNFHMMER_API_IPS").split(',')
        self.last_index_used = 0
        self.phdPrefix = os.getenv("PHD_DB_PREFIX")

    def meets_inclusion_threshold(self, cath_fun_fam: str, bitscore: float, escore: float, verbose: bool = False,
                                  output: list = []):
        try:
            checksql = f"""
                SELECT CASE WHEN (%s >= cff.cathfunfamily_bitscore_threshold) -- AND %s <= cff.cathfunfamily_evalue)
                               THEN 1 ELSE 0 END AS "include", cff.cathfunfamily_bitscore_threshold AS "bitscore",
                                     sci(cff.cathfunfamily_evalue) AS "evalue"
                FROM {self.phdPrefix}.GoGraph_cathfunfamilies cff
                WHERE cff.cathfunfamilyfull_id = %s
            """

            connection = MySQLConnectionPool.get_instance().get_connection()
            include_protein = False

            with connection.cursor() as cursor:
                cursor.execute(checksql, (bitscore, escore, cath_fun_fam,))

                row = cursor.fetchone()

                if verbose:
                    output.append(f"For <{cath_fun_fam}>, <evalue: {row['evalue']}>, <bitscore: {row['bitscore']}>")

                if row is None:
                    include_protein = False
                else:
                    include_protein = True if row['include'] == 1 else False

            MySQLConnectionPool.get_instance().close_connection(connection)

            return include_protein
        except:
            import sys
            exc_type, exc_value, exc_traceback = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_traceback)

            from utilities.emailmanager import EmailManager
            EmailManager.send_message('joseph.bonello@um.edu.mt', 'Unexpected Error', "\r\n".join(["In funfhmmer.py", traceback.format_exc(), ""]))

        return False

    def fhmmer_search_online(self, sequence: str, verbose: bool, timeout: int, output: list = []):
        try:
            matches = set()

            params = {'fasta': sequence}
            url = 'http://www.cathdb.info/search/by_funfhmmer'
            headers = {"Accept": "application/json"}

            # call post service with headers and params
            with requests.post(url, headers=headers, data=params) as response:

                if output is None:
                    verbose = False

                if verbose:
                    output.append("code: {}".format(str(response.status_code)))
                    output.append("headers: {}".format(str(response.headers)))
                    output.append("content: {}".format(str(response.text)))

                if response.text is not None or response.text != '':
                    resp_body = json.loads(response.text)
                else:
                    return None

                if 'task_id' in resp_body:
                    taskid = resp_body['task_id']
                else:
                    return None

                if verbose:
                    output.append("Found task id: {}".format(taskid))
                    output.append('-' * 30)

                # Get status of request
                message = ''
                count = 1
                while message != 'done' and count <= timeout:
                    url = 'http://www.cathdb.info/search/by_funfhmmer/check/{}'.format(taskid)
                    headers = {"Accept": "application/json"}

                    # call post service with headers and params
                    with requests.get(url, headers=headers) as resp:  # , data=params)
                        if verbose:
                            output.append("code: {}".format(str(resp.status_code)))
                            # output.append("headers: {}".format(str(resp.headers)))
                            # output.append("content: {}".format(str(resp.text)))

                        resp_body = json.loads(resp.text)

                        # This is to fix the issue with Unexpected Error in search online
                        # This seems to happen if there is a BLAST error (on website)
                        # Test with this sequence:
                        # MGSLTFRDVAIEFSLEEWQCLDTAQQNLYRNVMLENYRNLVFLGIAAFKPDLIIFLEEGKESWNMKRHEMVEESPVICSHFAQDLWPEQGIEDSFQKVILRRYEKCGHENLHLKIGYTNVDECKVHKEGYNKLNQSLTTTQSKVFQRGKYANVFHKCSNSNRHKIRHTGKKHLQCKEYVRSFCMLSHLSQHKRIYTRENSYKCEEGGKAFNWSSTLTYYKSAHTGEKPYRCKECGKAFSKFSILTKHKVIHTGEKSYKCEECGKAFNQSAILTKHKIIHTGEKPNKCEECGKAFSKVSTLTTHKAIHAGEKPYKCKECGKAFSKVSTLITHKAIHAGEKPYKCKECGKAFSKFSILTKHKVIHTGEKPYKCEECGKAYKWPSTLSYHKKIHTGEKPYKCEECGKGFSMFSILTKHEVIHTGEKPYKCEECGKAFNWSSNLMEHKKIHTGETPYKCEECGKGFSWSSTLSYHKKIHTVEKPYKCEECGKAFNQSAILIKHKRIHTGEKPYKCEECGKTFSKVSTLTTHKAIHAGEKPYKCKECGKTFIKVSTLTTHKAIHAGEKPYKCKECGKAFSKFSILTKHKVIHTGEKPYKCEECGKAFNWSSNLMEHKRIHTGEKPYKCEECGKSFSTFSVLTKHKVIHTGEKPYKCEECGKAYKWSSTLSYHKKIHTVEKPYKCEECGKAFNRSAILIKHKRIHTDEKPYKCEECGKTFSKVSTLTTHKAIHAGEKPYKCKECGKAFSKFSILTKHKVIHTGEKPYKCEECGKAYKWPSTLSYHKKIHTGEKPYKCEECGKGFSMFSILTKHEVIHTGEKPYKCEECGKAFSWLSVFSKHKKTHAGEKFYKCEACGKAYNTFSILTKHKVIHTGEKPYKCEECGKAFNWSSNLMEHKKIHTGETPYKCEECDKAFSWPSSLTEHKATHAGEKPYKCEECGKAFSWPSRLTEHKATHAGEEPYKCEECGKAFNWSSNLMEHKRIHTGEKPYKCEECGKSFSTFSILTKHKVIHTGEKPYKCEECGKAYKWSSTLSYHKKIHTVEKPYKCEECGKGFVMFSILAKHKVIHTGEKLYKCEECGKAYKWPSTLRYHKKIHTGEKPYKCEECGKAFSTFSILTKHKVIHTGEKPYKCEECGKAFSWLSVFSKHKKIHTGVPNPPTHKKIHAGEKLYK

                        if resp_body is None:
                            count += 1
                            time.sleep(2)
                            continue
                        else:
                            if not 'message' in resp_body:
                                count += 1
                                time.sleep(2)
                                continue

                        message = resp_body['message']
                        message = '' if message is None else message

                        if verbose:
                            output.append("Query is ... {}".format(message))
                            output.append('-' * 30)

                        time.sleep(5)

                        count = count + 1

                        # Get result of request
                        if message == 'done':
                            url = 'http://www.cathdb.info/search/by_funfhmmer/results/{}'.format(taskid)
                            headers = {"Accept": "application/json"}

                            # call post service with headers and params
                            resp = requests.get(url, headers=headers)  # , data=params)

                            if verbose:
                                output.append("code: {}".format(str(resp.status_code)))
                                # output.append("headers: {}".format(str(resp.headers)))
                                # output.append("content: {}".format(str(resp.text)))

                            if resp.text.strip() == '':
                                return None

                            resp_body = json.loads(resp.text)

                            fun_fam_scan = resp_body['funfam_scan']['results']

                            if verbose:
                                output.append('Scan Results: ')

                            for scanItem in fun_fam_scan:
                                for hit in scanItem['hits']:
                                    hsps = hit['hsps']
                                    hsp = list(hsp for hsp in hsps if hsp["rank"] == 1)[0]

                                    funfam = hit['match_id'].replace('FF/', 'FF').replace('/', '.')

                                    include_protein = self.meets_inclusion_threshold(funfam, hsp['score'], hsp['evalue'])

                                    if verbose:
                                        output.append(f"^<Match id: {hit['match_id']}>, <Bit Score: {hsp['score']}>, "
                                                      f"<E-Value: {hsp['evalue']}>, <include: {include_protein}>^")

                                    if include_protein:
                                        matches.add(funfam)
            return matches
        except OSError as err:
            print("OS error: {0}".format(err))
            print("Faulty sequence is: {}".format(sequence))
        except ValueError as ve:
            # Log Value Errors, but do not notify
            print(f"Value error: {ve}")
            print("Faulty sequence is: {}".format(sequence))

            # from utilities.emailmanager import EmailManager
            # EmailManager.send_message('joseph.bonello@um.edu.mt', 'Value Error',
            #                                   "\r\n".join(["In funfhmmer.py", traceback.format_exc(),
            #                                                        "", "Current sequence is", sequence]))
        except:
            import sys
            exc_type, exc_value, exc_traceback = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_traceback)
            print("Faulty sequence is: {}".format(sequence))

            from utilities.emailmanager import EmailManager
            EmailManager.send_message('joseph.bonello@um.edu.mt', 'Unexpected Error',
                                              "\r\n".join(["In funfhmmer.py", traceback.format_exc(),
                                                           "", "Current sequence is", sequence]))

        return None

    def fhmmer_search_offline(self, protein_description: str, sequence: str, output: list = [], verbose: bool = False):
        import os
        import tempfile
        import random
        import string

        unique = ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(5))

        basepath = os.path.join(tempfile.gettempdir(), unique)
        try:
            if output is not None:
                output.append("Starting offile scan")

            os.mkdir(basepath)
            input_path = os.path.join(basepath,  "seq.fasta")

            fasta = open(input_path, "w")
            if not ">" in sequence:
                fasta.write(">{}\n".format(protein_description))
            fasta.write(sequence)

            fasta.flush()
            fasta.close()

            matches = set()

            # ./apps/cath-genomescan.pl -i data/test.fasta -l data/funfam-hmm3-v4_1_0.lib -o results/
            import socket
            genome_scan_path = os.getenv("GENOME_SCAN_PATH_" + socket.gethostname().upper())

            cath_tools_process = os.path.join(genome_scan_path, 'cath-genomescan.pl')
            cath_tools_data_path = os.getenv('HMMLIB_' + socket.gethostname().upper())

            import subprocess

            subprocess.check_call([cath_tools_process, '-i', input_path, '-l', cath_tools_data_path, '-o',
                                   basepath])

            output_path = os.path.join(basepath, "seq.crh")
            cath_output = open(output_path, "r")

            for line in cath_output:
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

            if output is not None:
                output.append("Offline scan completed with {} results.".format(len(matches)))
            return matches
        except OSError as err:
            print("OS error: {0}".format(err))
            print("Faulty sequence is: {}".format(sequence))
        except ValueError:
            print("Could not convert data to an integer.")
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
                                              "\r\n".join(["In funfhmmer.py", traceback.format_exc(),
                                                           "", "Current sequence is", sequence]))
        finally:
            if os.path.exists(basepath):
                shutil.rmtree(basepath, ignore_errors=True)

        return None

    def fhmmer_search_api(self, protein_description: str, sequence: str, output: list, verbose=True, timeout=36):
        url = self.url_list[self.last_index_used]
        self.last_index_used += 1
        if self.last_index_used >= len(self.url_list):
            self.last_index_used = 0  # recycle

        try:
            matches = set()

            if not sequence.startswith(">" + protein_description):
                sequence = ">" + protein_description + '\n' + sequence

            params = {'command': 'funfams', 'sequence': sequence}
            headers = {"Accept": "application/json"}

            # call post service with headers and params
            print(f"DEBUG: Getting Funfams from {url}")
            with requests.put(url, headers=headers, data=params) as response:

                if output is None:
                    verbose = False

                if verbose:
                    output.append("code: {}".format(str(response.status_code)))
                    output.append("headers: {}".format(str(response.headers)))
                    output.append("content: {}".format(str(response.text)))

                response = json.loads(response.text)

                if response is not None:
                    if response['funfams'] is not None:
                        for hit in response['funfams'].split(','):
                            matches.add(hit.replace('-', '.'))
                    else:
                        return None
                else:
                    return None

                return matches
        except OSError as err:
            print("OS error: {0}".format(err))
            print("Faulty sequence is: {}".format(sequence))
        except ValueError as ve:
            # Log Value Errors, but do not notify
            print(f"Value error: {ve}")
            print("Faulty sequence is: {}".format(sequence))

            # from utilities.emailmanager import EmailManager
            # EmailManager.send_message('joseph.bonello@um.edu.mt', 'Value Error',
            #                                   "\r\n".join(["In funfhmmer.py", traceback.format_exc(),
            #                                                        "", "Current sequence is", sequence]))
        except:
            import sys
            exc_type, exc_value, exc_traceback = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_traceback)
            print("Faulty sequence is: {}".format(sequence))

            from utilities.emailmanager import EmailManager
            EmailManager.send_message('joseph.bonello@um.edu.mt', 'Unexpected Error',
                                              "\r\n".join(["In funfhmmer.py", traceback.format_exc(),
                                                           "", "Current sequence is", sequence]))

        return None

    def fhmmer_search(self, protein_description: str, sequence: str, output: list, verbose=True, timeout=36):
        results = None

        # if len(sequence) < 1500:
        #     print("Using Online FunFhmmer")
        #     results = self.fhmmer_search_online(sequence, verbose, timeout, output)
        # else:
        #     results = None

        # if results is None:
        #     results = self.fhmmer_search_offline(protein_description, sequence, output)

        results = self.fhmmer_search_api(protein_description, sequence, output)
        if results is None:  # Attempt a search online, just in case
            # if len(sequence) < 1500:
            print("Using Online FunFhmmer")
            results = self.fhmmer_search_online(sequence, verbose, timeout, output)

            # print("Using FunFhmmer Orengo Group Protocol")
            # results = self.fhmmer_search_orengo_group_protocol(protein_description, sequence, output, verbose, timeout)

        return results

    def fhmmer_search_orengo_group_protocol(self, protein_description: str, sequence: str, output: list, verbose=True, timeout=36):
        import os
        import socket
        hmm_location = os.getenv('CATH_HMM_LOCATION_' + socket.gethostname().upper())
        cath_resolve_hits_path = os.getenv('CATH_RESOLVE_HITS_PATH_' + socket.gethostname().upper())

        results = []

        import string
        import random
        unique = ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(5))

        import tempfile
        basepath = os.path.join(tempfile.gettempdir(), unique)
        try:
            if output is not None:
                output.append("Starting offile scan")

            os.mkdir(basepath)
            input_path = os.path.join(basepath, "seq.fasta")

            fasta = open(input_path, "w")
            if not ">" in sequence:
                fasta.write(">{}\n".format(protein_description))
            fasta.write(sequence)

            fasta.flush()
            fasta.close()

            output_path = os.path.join(basepath, "seq.result")
            resolve_hits_output_path = os.path.join(basepath, "seq.crh")

            import subprocess

            subprocess.check_call(['hmmsearch', '-Z', '10000000', '--domE', '0.001', '--incdomE', '0.001', '--cpu', '5', '-o', output_path, hmm_location, input_path])

            with open(resolve_hits_output_path, 'w+') as rh_out:
                subprocess.call(
                    [cath_resolve_hits_path, '--min-dc-hmm-coverage=80', '--worst-permissible-bitscore=25', '--output-hmmer-aln', '--input-format', 'hmmsearch_out', output_path],
                    stdout=rh_out)

            import re
            regex = re.compile('[0-9\.]+\/FF\/[0-9]+')
            with open(resolve_hits_output_path, 'r') as rh_out:
                for line in rh_out:
                    match = regex.search(line)
                    if match is not None:
                        results.append(match.group(0).replace('FF/', 'FF').replace('/', '.'))

            return results
        except:
            import sys
            print("Unexpected error:", sys.exc_info()[0])
            print("Faulty sequence is: {}".format(sequence))

            from utilities.emailmanager import EmailManager
            EmailManager.send_message('joseph.bonello@um.edu.mt', 'Unexpected Error', "\r\n".join(["In funfhmmer.py", traceback.format_exc(), "", "Current sequence is", sequence]))
        # finally:
        #     if os.path.exists(basepath):
        #         shutil.rmtree(basepath, ignore_errors=True)

        return None
