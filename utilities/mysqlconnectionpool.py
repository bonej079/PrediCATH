import os
import sys
import threading
import traceback

import pymysql

"""
    Connection Pool for MySQL. Supports multiple connections to different servers.
"""


class MySQLConnectionPool:
    class __ConnectionPoolProperties__:
        def __init__(self, use_settings: bool = True, pool_name: str = '', servers: list = [], db_user: str = '',
                     db_password: str = '', db_port: int = 3306, db_name: str = '', initial_size=5, enlarge_by=5,
                     max_connections=800):
            self.__available_connections__ = []
            self.__used_connections__ = []
            self.__initial_size__ = initial_size
            self.__enlarge_by__ = enlarge_by
            self.__max_connections__ = max_connections
            self.pool_name = pool_name

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
        if MySQLConnectionPool.__instance__ is None:
            MySQLConnectionPool.__instance__ = MySQLConnectionPool()

        return MySQLConnectionPool.__instance__

    def __init__(self):
        MySQLConnectionPool.__pool__ = dict()
        self._key_lock = threading.Lock()

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
                        db_password: str = '', db_port: int = 3306, db_name: str = '', initial_size=5, enlarge_by=5, max_connections=800):
        pool_properties = MySQLConnectionPool.__ConnectionPoolProperties__(True, pool_name, servers,
                                                                           db_user, db_password, db_port,
                                                                           db_name, initial_size, enlarge_by,
                                                                           max_connections)

        # check ips
        ip_list = list(pool_properties.ips)
        for ip in ip_list:
            self.__check_ips__(ip, pool_properties)

        # Tuple position 1 = connection properties, position 2: list of connections
        if pool_name != '':
            MySQLConnectionPool.__pool__[pool_name] = pool_properties
        else:
            MySQLConnectionPool.__pool__['DEFAULT'] = pool_properties

        MySQLConnectionPool.__instance__ = self

    def __check_ips__(self, ip: str, pool_properties: __ConnectionPoolProperties__):
        try:
            pymysql.connect(host=ip,
                            port=pool_properties.db_port,
                            user=pool_properties.db_user,
                            passwd=pool_properties.db_password,
                            db=pool_properties.db_name,
                            cursorclass=pymysql.cursors.DictCursor,
                            autocommit=True)
        except Exception:
            if ip in pool_properties.ips:
                pool_properties.ips.remove(ip)
            exc_type, exc_value, exc_traceback = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_traceback)

    def __connect__(self, ip: str, pool_name: str = 'DEFAULT') -> pymysql.Connection:
        connection_properties = self.__pool__[pool_name]
        try:
            connection = pymysql.connect(host=ip,
                                         port=connection_properties.db_port,
                                         user=connection_properties.db_user,
                                         passwd=connection_properties.db_password,
                                         db=connection_properties.db_name,
                                         cursorclass=pymysql.cursors.DictCursor,
                                         autocommit=True)

            return connection
        except Exception:
            if ip in connection_properties.ips:
                connection_properties.ips.remove(ip)
            exc_type, exc_value, exc_traceback = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_traceback)
            return None

    def __create_connections__(self, pool_name: str = 'DEFAULT'):
        pool_properties = self.__pool__[pool_name]

        current_server_index = 0
        current_connections = len(pool_properties.__available_connections__) + len(pool_properties.__used_connections__)

        if current_connections == 0:
            num_connections = pool_properties.__initial_size__
        else:
            if current_connections + pool_properties.__enlarge_by__ <= pool_properties.__max_connections__:
                num_connections = pool_properties.__enlarge_by__
            else:
                num_connections = pool_properties.__max_connections__ - current_connections

        for i in range(0, num_connections):
            mysql_connection = self.__connect__(pool_properties.ips[current_server_index], pool_name)
            mysql_connection.ping(reconnect=False)

            pool_properties.__available_connections__.append(mysql_connection)

            if (current_server_index + 1) >= len(pool_properties.ips):
                current_server_index = 0
            else:
                current_server_index += 1

    def add_pool_manually(self, pool_name, servers: list = [], db_user: str = '',
                          db_password: str = '',
                          db_port: int = 3306, db_name: str = '', initial_size=5, enlarge_by=5, max_connections=800):
        if pool_name is None or pool_name == '':
            raise BaseException("Pool name cannot be none for manually added pools")

        pool_properties = MySQLConnectionPool.__ConnectionPoolProperties__(False, pool_name, servers,
                                                                           db_user, db_password, db_port,
                                                                           db_name, initial_size, enlarge_by,
                                                                           max_connections)

        # Tuple position 1 = connection properties, position 2: list of connections
        MySQLConnectionPool.__pool__[pool_name] = pool_properties

    def get_connection(self, none_if_unable=False, pool_name: str = 'DEFAULT') -> pymysql.Connection:
        try:
            self._key_lock.acquire()

            if pool_name not in self.__pool__.keys():
                MySQLConnectionPool.__instance__.__create_pool__(pool_name)

            pool_properties = self.__pool__[pool_name]

            if len(pool_properties.__available_connections__) == 0:
                self.__create_connections__(pool_name)

            if len(pool_properties.__available_connections__) > 0:
                mysql_connection = pool_properties.__available_connections__[0]
                mysql_connection.ping(reconnect=True)
                if mysql_connection in pool_properties.__available_connections__:
                    pool_properties.__available_connections__.remove(mysql_connection)
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
                        mysql_connection.ping(reconnect=True)

                        pool_properties.__available_connections__.remove(mysql_connection)
                        pool_properties.__used_connections__.append(mysql_connection)

                        return mysql_connection
                    else:
                        count_retries -= 1

                if none_if_unable:
                    return None
                else:
                    raise BaseException("Failed to get a MySQL Connection from pool - waiting period expired. "
                                        f"Close connections from {pool_name}, with ips {pool_properties.ips}.")
        finally:
            self._key_lock.release()

    def close_connection(self, mysql_connection: pymysql.Connection, pool_name: str = 'DEFAULT'):
        try:
            self._key_lock.acquire()
            if pool_name not in self.__pool__.keys():
                raise BaseException(f"Pool {pool_name} does not exist in Connection Pool List")

            pool_properties = self.__pool__[pool_name]
            # print(f"Closing Connection from pool {pool_properties.pool_name}")

            if mysql_connection is not None:
                if mysql_connection in pool_properties.__used_connections__:
                    pool_properties.__available_connections__.append(mysql_connection)
                    pool_properties.__used_connections__.remove(mysql_connection)
        finally:
            self._key_lock.release()
