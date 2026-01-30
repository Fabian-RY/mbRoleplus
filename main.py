#! /usr/bin/env python3

"""
	@author Fabián Robledo
	@mail fabian.robledo@csic.es

	This scripts offers an straightforward way to convert between different IDs from differente metabolomical Databases

	To achieve this goal, uses InChiKeys: If 2 IDs are associated to the same inchikey, they are considered the same 
	compound. InchiKeys and IDs have been gathered on November 2024 from the different sources, and is the same database
	that is available in Mbroleplus page

	The actual implementation uses a relational database to benefit from indexation and efficient storage of data. This
	script aims to provide an interface to an SQLite database, so that everyone can create their own conversor data or
	add new compounds easily. For example, if there are only a reduced set of databases of interest

	usage:
		- mbrole-converter.py (main.py) -i <input_txt> -o <output_txt> -db <sqlite3.db> -t <table_name> --id <id-name-in_db> -gzi -gzo
			-i: a plain text file with one-per-line compounds with input compounds. Mandatory
			-o: a plain text file path with one-per-line compounds generated wit all compounds in output ID. Default: print to STD,
			-db: a sqlite3 db with a table with the ID, INCHIKEY, Database data. Mandatory
			-t: name of the table to used. Default: mbroleplus
			-id: ID to convert to; must be one of the possibilities in the "Database" column of the table. Mandatory
			-gzi: Boolean. If input file (-i) is compressed with gzip. Default: False
			-gzo: Boolean. If output (-o) file should be compressed with gzip. Default: False
		
		It also generates 2 extra files, whose paths can be decided, but use the current execution folder when executed
			-l: A log of the conversion
			-d: Discarded compounds, or those not found in the DB.
	
	This script aims to keep comptaibility con python 3.0 onwards
"""

import argparse
import gzip
import logging
import os
import sqlite3
import sys

import tqdm

### Strings to be formated or used eventually
LOGGER_NAME = "mbroleplus"
SQL_TABLE_EXISTS = "SELECT name FROM sqlite_master WHERE type='table';"
SQLITE_COLUMNS_OF_TABLE = "PRAGMA table_info({})"
SQL_SELECT_SOURCES = "SELECT DISTINCT {} FROM {};"
SQL_GET_INCHIKEY = "SELECT inchikey from {} where id == \"{}\""
SQL_GET_ID_WITH_INCHIKEY = "Select id from {} where inchikey == \"{}\" AND {} == \"{}\""

def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input","-i", type=str, help="A one-per-line compound file to be converted")
    parser.add_argument("--output","-o", type=str, default="/dev/stdout", help="A one-per-line file with the same compounds but converted to the same ID")
    parser.add_argument("--database","-db",type=str, help="Sqlite3 database with the conversion info")
    parser.add_argument("--table","-t", type=str, default="mbroleplus", help="Table name inside the database to use")
    parser.add_argument("--field","-f", type=str, help="Field in which the sources are stored", default="database")
    parser.add_argument("--id","-id", type=str, help="ID to convert input to. Must exist into the databse DATABASE column")
    parser.add_argument("--gzipped_input","-gzi", type=bool, default=False, help="Whether the input file is gzipped or not. Default: False")
    parser.add_argument("--gzipped_output","-gzo", type=bool, default=False, help="Should the output file be gzipped or not. Default: False")
    parser.add_argument("--logfile","-l",type=str, default="/dev/stdout", help="Where to save the log file. Default: print to stdout")
    parser.add_argument("--loglevel","-ll",type=str, choices=["debug","info","warn","error","critical"], default="info", help="Minimum level to show in the log. Default: Info")
    parser.add_argument("--discarded_file","-d", type=str, default=".", help="Where the compounds not available to the database should be stored")
    return parser.parse_args()

def _parse_loglevel(loglevel):
    """
        Turns loglevel flag into a logging value
    """
    if (loglevel == "debug"):
        return logging.DEBUG
    elif (loglevel == "info"):
        return logging.INFO
    elif (loglevel == "warn"):
        return logging.WARNING
    elif (loglevel == "error"):
        return logging.ERROR
    elif (loglevel == "critical"):
        return logging.CRITICAL
    else:
		# This should not happen, as only 5 choices are presented in argparse
		# Anyways, just in case, print a message and leave to avoid problems
        print(f"Error: invalid loglevel {loglevel}. Valid values are: debug, info, warn, error, critical")
        sys.exit(-1)

def _set_logger(name):
	"""
		Sets a logger with a defined name and minimum level.
		By default prints to STDOUT, but accepts a filename,
		In which case, will print to both STDOUT and the filename,
		with will be created if does not exist and will append the
		data into the file, and will not erase it.

		TODO: add file compatibility. For now, stdot can be redirected
	"""
	logger = logging.getLogger(name)
	return logger

def _parse_input_file(file, gz):
    """
		Convert the input files into a set.
		Assumes every line is a different compound name/ID

		returns a set
    """
    logger = _set_logger(_parse_input_file.__name__)
    openf = gzip.open if gz else open
    with openf(file) as fhand:
        logger.info(f"Parsing file: {file}")
        input_set = {x.strip("\n") for x in fhand.readlines()}
    return input_set

def validation(id, database, table, input_file, output_file, logfile, field):
    """
		Perform some checks and validates input, in particular:
			- Whether the path given for database exists
			- Whether the table given exists in the database
			- Whether the ID is present in the table
    """
    # Set a new logger to Identify this section in the log
    logger = _set_logger(validation.__name__)
    
    ## Checking if the database given exists
    if not os.path.exists(database):
        logger.critical(f" Not Found: SQLite database {database}. Pherhaps a typo or incorrect path?")
        return False
    logger.info(f" Database {database} found.")
    
    ## Checking if table exist in database
    connector = sqlite3.Connection(database)
    cursor = connector.cursor()
    cursor.execute(SQL_TABLE_EXISTS)
    tables = [x[0] for x in cursor.fetchall()] # pySQLite returns [(table_name,),(table_name2, )] a.k.a a list of tuples, so we extract the meaningful info
    if (table not in tables):
        # Table not in DB, validation unsuccessful
        logging.critical(f" Table {table} not found in {database}. Pherhaps a typo or incorrect path?")
        return False

    ## Checking if the field indicated exists
    logger.info(f" Table {table} found")
    cursor.execute(SQLITE_COLUMNS_OF_TABLE.format(table))
    columns = [x[1] for x in cursor.fetchall()]
    if (field not in columns):
        # Field not in Table, validation unsuccessful
        logging.critical(f" Field {field} not found in table {table}. Pherhaps a typo?")
        return False

    ## Finally, checking if the desired source exists in the database
    cursor.execute(SQL_SELECT_SOURCES.format(field, table))
    sources = [x[0] for x in cursor.fetchall()]
    if (id not in sources):
        # Source not in DB, validation unsuccessful
        logger.critical(f"Source {id} not found in {table} in column {field}. Pherhaps a typo?")
        return False
    logger.info(f" Source {id} found")
    if (not os.path.exists(input_file)):
        logger.critical( f"Input file: {input_file} not found. Check typos or correct the path")
    logging.info(f" Input file {input_file} found")
    if (not os.path.exists(os.path.dirname(output_file))):
        logger.critical( f"Output folder: {input_file} not found. Check typos or correct the path")
    logging.info(f" Output file {output_file} can be created")

    logger.info(" Validation finished, proceeding with conversion")
    return True

def convert(input, database, table, id, field):
    """
        Turns an ID into another ID.

        We use the inchikey as a middleman. Inchikeys are unique for every molecule, 
        so if 2 IDs share inchikey, they represent the same metabolite. Searching for one 
        should return the other.

        The table schema is simple: each row contains inchikey, ID and database.
        Searching by ID we obtain the inchikey, and using the inchikey we can 
        recieve ALL the metabolites with that inchikey. Then we can filter by database
        or use them as they are returned.


    """

    logger = _set_logger(convert.__name__)

    # Starting connection to the database
    # Cursor will be used for queries
    connection = sqlite3.Connection(database)
    cursor = connection.cursor()
    
    logger.info(f" Querying {len(input)} metabolites")
    
    # Using dict and set because keys will be unique: 
    # Even if the user provided twice or thrice the same
    # Key, the conversion will be the same and give only 
    # one output. Parsing the input file as a set had
    # already removed duplicates
    founded_ids = dict()
    unavailable_ids = set()

    # TQDM to show a progress bar
    # Its usefullness will depend on how many
    # metabolites given as input
    for metabolite in tqdm.tqdm(input):
        SQL = SQL_GET_INCHIKEY.format(table, metabolite)
        logging.debug(f" Executing {SQL}")
        results = cursor.execute(SQL).fetchall()
        if (len(results) == 0):
            logger.warning(f" Metabolite with ID {metabolite} not found")
            unavailable_ids.add(metabolite)
        else:
            inchikey = results[0][0]
            logger.debug(f" Found inchikey {inchikey} with ID  {metabolite}")
            SQL = SQL_GET_ID_WITH_INCHIKEY.format(table, inchikey, field, id)
            logger.debug(SQL)
            molecules = cursor.execute(SQL).fetchall()[0][0]
            if (len(molecules) == 0):
                # TODO: What happens when there no ID in destination
                pass
            founded_ids[metabolite] = molecules
    return (founded_ids, unavailable_ids)

def save_compounds(file, gzipped, compounds):
    openf = gzip.open if gzipped else open
    with (openf(file, "wt")) as fhand:
        list(map(lambda x: print(x, file=fhand, end="\n"), compounds))


def main():
    args = _parse_args()
    loglevel = _parse_loglevel(args.loglevel)
    logging.basicConfig(filename=args.logfile, level=loglevel)
    logger = _set_logger(LOGGER_NAME)
    if not validation(args.id, args.database, args.table, args.input, args.output, args.logfile, args.field):
        logger.info("Input commands not valid. Exiting...")
        sys.exit(1)
    input_compounds = _parse_input_file(args.input,  gz=args.gzipped_input)
    converted, unavailable = convert(input_compounds, args.database, args.table, args.id, args.field)
    logger.info(f" Found {len(converted.keys())} IDs ({len(converted.keys())/len(input_compounds)}%)")
    logger.info(f" {len(unavailable)} metabolites not found ({len(unavailable)/len(input_compounds)})")
    save_compounds(args.output, args.gzipped_output, set(converted.values()))
    save_compounds(args.discarded_file, False, unavailable)

if __name__ == "__main__":
    main()