import os
import socket

from menu import Menu


def metapredictor_handler():
    output = ['Starting prediction of proteins in file']
    withDcgo = False

    for i in range(0, 2):
        if withDcgo:
            indexes_to_use = ['Jaccard', 'Sorensen', 'Overlap', 'dcgo']
        else:
            indexes_to_use = ['Jaccard', 'Sorensen', 'Overlap']
        ontologies = ['BP', 'MF', 'CC']

        s = input("Set Set-Based Methods threshold (default = 0.5): ")
        threshold = 0.5 if (s is None or not s.replace('.', '', 1).isdigit()) else float(s)
        inherit_go_terms = True

        from GoGraph.classes.metapredictor import MetaPredictor
        metapredictor = MetaPredictor(indexes_to_use, ontologies, threshold, inherit_go_terms)

        path = os.path.join(os.getenv('JBPHD_DATA_UOM-1A26'), os.getenv('PREDICTION_PATH_UOM-1A26'),
                            'cafa2_bacteria_160488.tfa')

        file_handle = open(path, 'r')
        sequences = list()
        sequence = ''

        for line in file_handle:
            if line.strip() == '' or line.find('>') > -1:  # Handle new sequences; either through a line break
                # or through the discovery of a new sequence right after this one ends
                if sequence.strip():
                    sequences.append(sequence)  # Checks if empty
                if line.find('>') > -1:
                    sequence = line
            else:
                if line.find('>') > -1:
                    sequence += line
                    # print(line)
                else:
                    sequence += line.strip()

        metapredictor.predict_list_of_sequences(sequences)

        withDcgo = True

        # print(output)


def metapredictor_from_path_handler():
    output = ['Starting prediction of proteins in file']
    withDcgo = True
    from GoGraph.classes.metapredictor import MetaPredictor

    # s = input("Set Set-Based Methods threshold (default = 0.8): ")
    # threshold = 0.8 if (s is None or not s.replace('.', '', 1).isdigit()) else float(s)

    import numpy as np
    for threshold in np.arange(0.1, 1.1, 0.1):
        for i in range(0, 2):
            if withDcgo:
                indexes_to_use = ['Jaccard', 'Sorensen', 'Overlap', 'dcgo']
            else:
                indexes_to_use = ['Jaccard', 'Sorensen', 'Overlap']
            ontologies = ['BP', 'MF', 'CC']

            inherit_go_terms = True

            meta_predictor = MetaPredictor(indexes_to_use, ontologies, threshold, inherit_go_terms, skip_not_in_cache=True, log_to_file=True, aggressive_debug=True)

            import socket
            datapath = os.getenv('JBPHD_DATA_' + socket.gethostname().upper())
            prediction_path = os.getenv('PREDICTION_PATH_' + socket.gethostname().upper())
            thread_model = os.getenv('THREAD_MODEL', 'threading')

            # path = os.path.join(datapath, prediction_path, 'cafa2_bacteria_160488.tfa')
            # path = os.path.join(datapath, 'CAFA/CAFA1/Targets/CAFA_2010_Euk_targets')

            path = os.path.join(datapath, prediction_path)
            meta_predictor.predict_from_path(path, parallel_backend=thread_model, use_threads=True)

            withDcgo = not withDcgo
            # print(output)


def predictions_handler():
    # sequences = {
    #    "MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQKMEPRAPWIEQEGPEYWDQETRNMKAHSQTDRANLGTLRGYYNQSEDGSHTIQIMYGCDVGPDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAVHAAEQRRVYLEGRCVDGLRRYLENGKETLQRTDPPKTHMTHHPISDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWELSSQPTIPIVGIIAGLVLLGAVITGAVVAAVMWRRKSSDRKGGSYTQAASSDSAQGSDVSLTACKV"
    #             }

    output = ['Starting prediction of proteins in file']
    import os

    # Parameters: threshold, inheritGoTerms, separateOntologies, indexesToUse, databaseProperties
    s = input("Set Set-Based Methods threshold (default = 0.5): ")
    threshold = 0.5 if (s is None or not s.replace('.', '', 1).isdigit()) else float(s)
    inherit_go_terms = True
    separate_ontologies = True
    indexes_to_use = ['Jaccard', 'Sorensen', 'Overlap', 'dcgo']

    thread_model = os.getenv('THREAD_MODEL', 'threading')

    # from predictions import Predictions_CATH
    from GoGraph.classes.predictions_cath import Predictions_CATH
    predictor = Predictions_CATH(threshold, inherit_go_terms, separate_ontologies, indexes_to_use, use_rollback=False, apply_cache_update_to_db=True)

    import socket
    datapath = os.getenv('JBPHD_DATA_' + socket.gethostname().upper())
    prediction_path = os.getenv('PREDICTION_PATH_' + socket.gethostname().upper())

    # path = os.path.join(datapath, 'CAFA/CAFA1/Targets/CAFA_2010_Euk_targets')
    path = os.path.join(datapath, prediction_path)
    predictor.predict_from_path(path, output, parallel_backend=thread_model)

    print(output)


def runMenu(filepath: str):
    try:
        env_path = os.path.join(filepath, '..', '..', 'PhDCode', 'environment.env')

        from dotenv import load_dotenv
        load_dotenv(verbose=True, dotenv_path=env_path)

        cath_version = os.getenv('CATH_VERSION_' + socket.gethostname().upper())

        env_path = os.path.join(filepath, '..', '..', 'PhDCode', f'environment-{cath_version}.env')
        print(f"DEBUG: Loading CATH DB settings from {env_path}")
        load_dotenv(verbose=True, dotenv_path=env_path)

        pfam_version = os.getenv('PFAM_VERSION_' + socket.gethostname().upper())

        env_path = os.path.join(filepath, '..', '..', 'PhDCode', f'environment-{pfam_version}.env')
        print(f"DEBUG: Loading PFAM DB settings from {env_path}")
        load_dotenv(verbose=True, dotenv_path=env_path)

        print(os.getenv('HMMLIB_' + socket.gethostname().upper()))

        mode = os.getenv('MODE')

        if mode == 'DEV':
            print('DEV mode is not functional in this version')
        else:
            mainmenu = Menu(title=f"PhD Menu Items - Working with CATH {cath_version.replace('CATH', '').replace('4', '4.')}", prompt='>',
                            options=[('Meta Predictor', metapredictor_handler), ('Meta Predictor From Path', metapredictor_from_path_handler),
                                     ('Perform predictions', predictions_handler),
                                     ("Quit", Menu.CLOSE)])

        # Menu()

        mainmenu.open()
    except:
        import traceback, sys
        print("Unexpected error:", sys.exc_info()[0])
        print(traceback.format_exc())


def main():
    # Fix for Cython path
    path = os.getcwd() if (os.getcwd().endswith('classes')) else os.path.join(os.getcwd(), 'GoGraph', 'classes')

    runMenu(path)


if __name__ == '__main__':
    main()
