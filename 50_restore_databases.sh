#!/bin/bash

cd /home/jbon4/Databases
for file in *.sql; do 
	databaseName=$(echo $file | cut -f 1 -d '.')
	databaseName=$(echo $databaseName | sed 's/^[0-9][0-9]-//1')
	echo $databaseName
	mysql -u root --password=password -h 172.17.0.3 --port=3306 -e "CREATE OR REPLACE DATABASE $databaseName;"	
	mysql -u root --password=password -h 172.17.0.3 --port=3306 $databaseName < $file
done

mysql -u root --password=password -h 172.17.0.3 --port=3306 -e "GRANT ALL PRIVILEGES ON *.* TO 'root'@'172.17.0.1' IDENTIFIED BY 'password';"
mysql -u root --password=password -h 172.17.0.3 --port=3306 -e "GRANT ALL PRIVILEGES ON *.* TO 'predicath'@'%' IDENTIFIED BY 'predicath';"