import os
import sys
import traceback

import mysql.connector as connector

"""
    Connection Pool for MySQL. Supports multiple connections to different servers.
"""


class OracleMySQLConnectionPool:
    class __ConnectionPoolProperties__:
        def __init__(self, use_settings: bool = True, pool_name: str = '', servers: list = [], db_user: str = '',
                     db_password: str = '', db_port: int = 3306, db_name: str = '', initial_size=5, enlarge_by=5,
                     max_connections=800):
            self.__available_connections__ = []
            self.__used_connections__ = []
            self.__initial_size__ = initial_size
            self.__enlarge_by__ = enlarge_by
            self.__max_connections__ = max_connections

            if use_settings:
                if pool_name is not None and pool_name != 'DEFAULT' and pool_name != '':
                    pool_name = pool_name.upper()

                    if os.getenv('MYSQL_DATABASE_NAME' + '_' + pool_name) is None:
                        raise BaseException(f"Requested pool {pool_name} does not have the required settings.")

                    iplist = os.getenv('MYSQL_HOST' + '_' + pool_name).split(',')
                    if len(iplist) > 0:
                        self.ips = [ip.strip() for ip in iplist]  # Clean the IPs just in case
                    else:
                        self.ips = list(os.getenv('MYSQL_HOST' + '_' + pool_name))
                    self.db_user = os.getenv('MYSQL_USER' + '_' + pool_name)
                    self.db_password = os.getenv('MYSQL_PASSWORD' + '_' + pool_name)
                    self.db_port = int(os.getenv('MYSQL_PORT' + '_' + pool_name))
                    self.db_name = os.getenv('MYSQL_DATABASE_NAME' + '_' + pool_name)
                else:
                    iplist = os.getenv('MYSQL_HOST').split(',')
                    if len(iplist) > 0:
                        self.ips = [ip.strip() for ip in iplist]  # Clean the IPs just in case
                    else:
                        self.ips = list(os.getenv('MYSQL_HOST'))
                    self.db_user = os.getenv('MYSQL_USER')
                    self.db_password = os.getenv('MYSQL_PASSWORD')
                    self.db_port = int(os.getenv('MYSQL_PORT'))
                    self.db_name = os.getenv('MYSQL_DATABASE_NAME')
            else:
                self.ips = [ip.strip() for ip in servers]  # Clean the IPs just in case
                self.db_user = db_user
                self.db_password = db_password
                self.db_port = db_port
                self.db_name = db_name

    __pool__ = dict()
    __instance__ = None

    # Ensure a Singleton instance
    @staticmethod
    def get_instance():
        if OracleMySQLConnectionPool.__instance__ is None:
            OracleMySQLConnectionPool.__instance__ = OracleMySQLConnectionPool()

        return OracleMySQLConnectionPool.__instance__

    def __init__(self):
        OracleMySQLConnectionPool.__pool__ = dict()

    def __exit__(self, exc_type, exc_val, exc_tb):
        print('Exiting. Closing pool connections')
        for pool in self.__pool__:
            for cnx in pool.__available_connections__:
                try:
                    cnx.close()
                except Exception:
                    exc_type, exc_value, exc_traceback = sys.exc_info()
                    traceback.print_exception(exc_type, exc_value, exc_traceback)

            for cnx in pool.__used_connections__:
                try:
                    cnx.close()
                except Exception:
                    exc_type, exc_value, exc_traceback = sys.exc_info()
                    traceback.print_exception(exc_type, exc_value, exc_traceback)

    """
        Initialises the connection pool. Note: pool_name is an optional parameter that is used in conjunction with
        use_settings (by default this is true). When used, the pool_name is *appended* to the parameter names in the
        environment file (e.g. MYSQL_HOST_<poolname>)
    """
    def __create_pool__(self, pool_name: str = '', servers: list = [], db_user: str = '',
                 db_password: str = '', db_port: int = 3306, db_name: str = '', initial_size=25, enlarge_by=5,
                     max_connections=800):
        pool_properties = OracleMySQLConnectionPool.__ConnectionPoolProperties__(True, pool_name, servers,
                                                                           db_user, db_password, db_port,
                                                                           db_name, initial_size, enlarge_by,
                                                                           max_connections)

        # check ips
        ip_list = list(pool_properties.ips)
        for ip in ip_list:
            self.__check_ips__(ip, pool_properties)

        # Tuple position 1 = connection properties, position 2: list of connections
        if pool_name != '':
            OracleMySQLConnectionPool.__pool__[pool_name] = pool_properties
        else:
            OracleMySQLConnectionPool.__pool__['DEFAULT'] = pool_properties

        OracleMySQLConnectionPool.__instance__ = self

    def __check_ips__(self, ip: str, pool_properties: __ConnectionPoolProperties__):
        try:
            connector.connect(host=ip,
                            port=pool_properties.db_port,
                            user=pool_properties.db_user,
                            password=pool_properties.db_password,
                            database=pool_properties.db_name,
                            # cursorclass=pymysql.cursors.DictCursor,
                            autocommit=True)
        except Exception:
            pool_properties.ips.remove(ip)
            exc_type, exc_value, exc_traceback = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_traceback)

    def __connect__(self, ip: str, pool_name: str = 'DEFAULT') -> connector.connection:  #  pymysql.Connection:
        connection_properties = self.__pool__[pool_name]
        try:
            connection = connector.connect(host=ip,
                                         port=connection_properties.db_port,
                                         user=connection_properties.db_user,
                                         password=connection_properties.db_password,
                                         database=connection_properties.db_name,
                                         # cursorclass=pymysql.cursors.DictCursor,
                                         autocommit=True)

            return connection
        except Exception:
            connection_properties.ips.remove(ip)
            exc_type, exc_value, exc_traceback = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_traceback)
            return None

    def __create_connections__(self, pool_name: str = 'DEFAULT', enl_by=0):
        pool_properties = self.__pool__[pool_name]

        current_server_index = 0
        current_connections = len(pool_properties.__available_connections__) + len(pool_properties.__used_connections__)

        if current_connections == 0:
            num_connections = pool_properties.__initial_size__
        else:
            if enl_by > 0:
                enlarge_by = enl_by
            else:
                enlarge_by = pool_properties.__enlarge_by__
            if current_connections + enlarge_by <= pool_properties.__max_connections__:
                num_connections = pool_properties.__enlarge_by__
            else:
                num_connections = pool_properties.__max_connections__ - current_connections

        for i in range(0, num_connections):
            mysql_connection = self.__connect__(pool_properties.ips[current_server_index], pool_name)

            pool_properties.__available_connections__.append(mysql_connection)

            if (current_server_index + 1) >= len(pool_properties.ips):
                current_server_index = 0
            else:
                current_server_index += 1

    def add_pool_manually(self, pool_name, servers: list = [], db_user: str = '',
                          db_password: str = '',
                          db_port: int = 3306, db_name: str = '', initial_size=5, enlarge_by=5, max_connections=150):
        if pool_name is None or pool_name == '':
            raise BaseException("Pool name cannot be none for manually added pools")

        pool_properties = OracleMySQLConnectionPool.__ConnectionPoolProperties__(False, pool_name, servers,
                                                                           db_user, db_password, db_port,
                                                                           db_name, initial_size, enlarge_by,
                                                                           max_connections)

        # Tuple position 1 = connection properties, position 2: list of connections
        OracleMySQLConnectionPool.__pool__[pool_name] = pool_properties

    def get_connection(self, none_if_unable=False, pool_name: str = 'DEFAULT') -> connector.connection:  # pymysql.Connection:
        if pool_name not in self.__pool__.keys():
            OracleMySQLConnectionPool.__instance__.__create_pool__(pool_name)

        pool_properties = self.__pool__[pool_name]

        if len(pool_properties.__available_connections__) == 0:
            self.__create_connections__(pool_name)

        if len(pool_properties.__available_connections__) > 0:
            ready = False
            while not ready:
                mysql_connection = pool_properties.__available_connections__[0]
                try:
                    mysql_connection.ping(reconnect=True, attempts=5, delay=3)
                    ready = True
                    pool_properties.__available_connections__.remove(mysql_connection)
                except:
                    if mysql_connection in pool_properties.__available_connections__:
                        pool_properties.__available_connections__.remove(mysql_connection)
                        self.__create_connections__(pool_name=pool_name, enl_by=1)
            pool_properties.__used_connections__.append(mysql_connection)

            return mysql_connection
        else:
            count_retries = 10
            import time
            import math

            while count_retries > 0:
                print(f'Waiting for: {(11 - count_retries) * math.exp(count_retries / 4)} seconds')
                time.sleep((11 - count_retries) * math.exp(count_retries / 4))

                if len(pool_properties.__available_connections__) > 0:
                    mysql_connection = pool_properties.__available_connections__[0]
                    pool_properties.__available_connections__.remove(mysql_connection)
                    pool_properties.__used_connections__.append(mysql_connection)

                    return mysql_connection
                else:
                    count_retries -= 1

            if none_if_unable:
                return None
            else:
                raise BaseException("Failed to get a MySQL Connection from pool - waiting period expired. "
                                    "Close connections.")

    # pymysql.Connection
    def close_connection(self, mysql_connection: connector.connection, pool_name: str = 'DEFAULT'):
        if pool_name not in self.__pool__.keys():
            raise BaseException(f"Pool {pool_name} does not exist in Connection Pool List")

        pool_properties = self.__pool__[pool_name]

        pool_properties.__available_connections__.append(mysql_connection)
        pool_properties.__used_connections__.remove(mysql_connection)
