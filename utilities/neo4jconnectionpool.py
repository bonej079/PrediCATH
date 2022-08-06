from py2neo import Graph
from py2neo.internal import http


class Neo4JConnectionPool:
    def __init__(self, servers: list, graph_user, graph_password, initial_size=5, enlarge_by=5, max_connections=150):
        self.__ips__ = [ip.strip() for ip in servers]  # Clean the IPs just in case
        self.__available_connections__ = []
        self.__used_connections__ = []
        self.__initial_size__ = initial_size
        self.__enlarge_by__ = enlarge_by
        self.__graph_user__ = graph_user
        self.__graph_password__ = graph_password
        self.__max_connections__ = max_connections

        # check ips
        ip_list = list(self.__ips__)
        for ip in ip_list:
            self.__connect__(ip)

    def __connect__(self, ip):
        try:
            graph_ip = ip
            graph_user = self.__graph_user__
            graph_password = self.__graph_password__

            # authenticate(graph_ip, graph_user, graph_password)

            graph = Graph(host=graph_ip, user=graph_user, password=graph_password, bolt=True)
            http.socket_timeout = 9999999

            return graph
        except:
            self.__ips__.remove(ip)
            return None

    def __create_connections__(self):
        current_server_index = 0
        current_connections = len(self.__available_connections__) + len(self.__used_connections__)

        if current_connections == 0:
            num_connections = self.__initial_size__
        else:
            if current_connections + self.__enlarge_by__ <= self.__max_connections__:
                num_connections = self.__enlarge_by__
            else:
                num_connections = self.__max_connections__ - current_connections

        for i in range(1, num_connections + 1):
            graph = self.__connect__(self.__ips__[current_server_index])

            self.__available_connections__.append(graph)

            if (current_server_index + 1) == len(self.__ips__):
                current_server_index = 0
            else:
                current_server_index += 1

    def get_connection(self, none_if_unable=False):
        if len(self.__available_connections__) == 0:
            self.__create_connections__()

        if len(self.__available_connections__) > 0:
            graph = self.__available_connections__[0]
            self.__available_connections__.remove(graph)
            self.__used_connections__.append(graph)

            return graph
        else:
            count_retries = 10
            import time
            import math

            while count_retries > 0:
                print(f'Waiting for: {(11 - count_retries) * math.exp(count_retries / 4)} seconds')
                time.sleep((11 - count_retries) * math.exp(count_retries / 4))

                if len(self.__available_connections__) > 0:
                    graph = self.__available_connections__[0]
                    self.__available_connections__.remove(graph)
                    self.__used_connections__.append(graph)

                    return graph
                else:
                    count_retries -= 1

            if none_if_unable:
                return None
            else:
                raise BaseException("Failed to get a Neo4J Connection from pool - waiting period expired. "
                                    "Close connections.")

    def close_connection(self, graph):
        self.__available_connections__.append(graph)
        self.__used_connections__.remove(graph)
