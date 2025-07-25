CREATE DATABASE pipelines_db;
CREATE ROLE dbuser WITH LOGIN ENCRYPTED PASSWORD 'dbpwd';
ALTER DATABASE pipelines_db OWNER TO dbuser;
CREATE DATABASE teaspoons_stairway_db;
CREATE ROLE stairwayuser WITH LOGIN ENCRYPTED PASSWORD 'stairwaypwd';
ALTER DATABASE teaspoons_stairway_db OWNER TO stairwayuser;
