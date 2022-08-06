from utilities.mysqlconnectionpool import MySQLConnectionPool


class UpdateTargets:

    def __init__(self):
        self.dbConnection = MySQLConnectionPool.get_instance().get_connection(pool_name="CACHE")
        import os
        self.phdCachePrefix = os.getenv("PHDCACHE_DB_PREFIX")

    def __exit__(self, exc_type, exc_val, exc_tb):
        MySQLConnectionPool.get_instance().close_connection(self.dbConnection)

    def update_from_path(self):
        import socket
        import os

        datapath = os.getenv('JBPHD_DATA_' + socket.gethostname().upper())
        prediction_path = os.getenv('PREDICTION_PATH_' + socket.gethostname().upper())

        # path = os.path.join(datapath, 'CAFA/CAFA1/Targets/CAFA_2010_Euk_targets')
        path = os.path.join(datapath, prediction_path)

        sql = f"""SELECT id FROM {self.phdCachePrefix}.pred_! WHERE hashed_sequence = %s"""
        sql_update_table = f"""ALTER TABLE {self.phdCachePrefix}.pred_! ADD COLUMN IF NOT EXISTS cafa_target VARCHAR(15)"""
        sql_update = f"""UPDATE {self.phdCachePrefix}.pred_! SET cafa_target = %s WHERE hashed_sequence = %s"""

        for file in os.listdir(path):
            count_sequences = 0
            table_name = file.replace(".tfa", "")
            print(file)

            sequences = dict()
            sequence = ''

            with self.dbConnection.cursor() as cursor:
                cursor.execute(sql_update_table.replace('!', table_name))

                with open(os.path.join(path, file), 'r') as file_handle:
                    for line in file_handle:
                        if line.strip() == '' or line.find('>') > -1:  # Handle new sequences; either through a line break
                            # or through the discovery of a new sequence right after this one ends
                            if sequence.strip() != '':
                                if target not in sequences.keys():
                                    count_sequences += 1
                                sequences[target] = sequence  # Checks if empty
                                sequence = ''
                            if line.find('>') > -1:
                                if sequence != '':
                                    if target not in sequences.keys():
                                        count_sequences += 1
                                    sequences[target] = sequence
                                target = line.split(' ')[0]
                                sequence = line
                        else:
                            sequence += line.strip()

                if sequence != '' and target != '':
                    if target not in sequences.keys():
                        count_sequences += 1
                    sequences[target] = sequence  # set will handle doubles; avoids mossing last sequence in the file

                if count_sequences != len(sequences):
                    raise BaseException(
                        f"There are not as many sequences in files as those processed {count_sequences} "
                        f"vs {len(sequences)}")

                for key in sequences.keys():
                    sequence = sequences[key]

                    from hashlib import blake2b
                    m = blake2b(digest_size=blake2b.MAX_DIGEST_SIZE, key=b"JBPhD", salt=b"PhD is SALTY")

                    m.update(sequence.encode())
                    hashed_sequence = m.hexdigest()

                    cursor.execute(sql.replace('!', table_name), (hashed_sequence,))
                    row_exists = cursor.fetchone() is not None

                    if row_exists:
                        cursor.execute(sql_update.replace('!', table_name), (key.replace('>', ''), hashed_sequence,))
                    else:
                        print(f"{key}: <{hashed_sequence}> - <{row_exists}>")
