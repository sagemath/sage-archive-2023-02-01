import sqlite3

from twisted.python import log

from sqlalchemy import *
from sqlalchemy.orm import sessionmaker, mapper, deferred, clear_mappers

from sage.dsage.misc.constants import DELIMITER
from sage.dsage.database.sql_functions import optimize_sqlite

CURRENT_SCHEMAVERSION = 2

def init_db_sa(db_file):
    """
    Initiates a database using SQLAlchemy.

    :type db_file: string
    :param db_file: filename of the database

    """

    from sage.dsage.database.schema import (metadata, jobs, clients, workers)
    from sage.dsage.database.job import Job
    from sage.dsage.database.client import Client
    from sage.dsage.database.worker import Worker

    def connect():
        conn = sqlite3.connect(db_file)
        conn.text_factory = str
        optimize_sqlite(conn)

        return conn

    # First we clear all predefined mappers
    clear_mappers()

    # engine = create_engine('sqlite:///%s' % (db_file), encoding='latin1',
    #                        echo=True)
    engine = create_engine('sqlite:///', creator=connect, echo=False)
    metadata.create_all(engine)

    mapper(Job, jobs, properties={'result': deferred(jobs.c.result),
                                  'data': deferred(jobs.c.data)})
    mapper(Client, clients)
    mapper(Worker, workers)

    Session = sessionmaker(bind=engine, autoflush=True, transactional=True)

    return Session

def init_db(db_file):
    db_conn = sqlite3.connect(
              db_file,
              isolation_level=None,
              detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
    optimize_sqlite(db_conn)
    db_conn.text_factory = str

    create_schema(db_conn)

    return db_conn

def create_schema_v1(con):
    """
    Create the initial version of the schema. The first release of your
    software only has this function for the schema.
    """

    cur = con.cursor()
    cur.executescript("""
        CREATE TABLE schema(version);
        INSERT INTO schema(version) VALUES (1);
    """)


    cur.executescript("""
        CREATE TABLE IF NOT EXISTS clients
            (
             id integer PRIMARY KEY,
             username text NOT NULL UNIQUE,
             public_key text NOT NULL UNIQUE,
             creation_time timestamp,
             access_time timestamp,
             last_login timestamp,
             connected BOOL,
             enabled BOOL DEFAULT 1
            );
        CREATE TABLE IF NOT EXISTS jobs
            (job_id TEXT NOT NULL UNIQUE,
             name TEXT,
             username TEXT REFERENCES clients(username),
             monitor_id TEXT REFERENCES monitors(uuid),
             worker_info TEXT,
             code TEXT,
             data BLOB,
             output TEXT,
             result BLOB,
             status TEXT NOT NULL,
             priority INTEGER DEFAULT 5,
             type TEXT,
             failures INTEGER DEFAULT 0,
             creation_time timestamp NOT NULL,
             update_time timestamp,
             start_time timestamp,
             finish_time timestamp,
             cpu_time REAL,
             wall_time REAL,
             verifiable BOOL,
             private BOOL DEFAULT 0,
             timeout INTEGER DEFAULT 0,
             killed BOOL DEFAULT 0
            );
        CREATE TABLE IF NOT EXISTS monitors
            (
             uuid text NOT NULL UNIQUE,
             hostname TEXT,
             ip TEXT,
             workers INTEGER,
             sage_version text,
             os text,
             kernel_version TEXT,
             cpus INTEGER,
             cpu_speed INTEGER,
             cpu_model TEXT,
             mem_total INTEGER,
             mem_free INTEGER,
             connected BOOL,
             busy BOOL,
             anonymous BOOL DEFAULT 0,
             last_connection timestamp
            );
    """)

def update_schema_to_v2(con):
    log.msg(DELIMITER)
    log.msg("Updating schema to v2...")
    log.msg(DELIMITER)
    cur = con.cursor()
    cur.executescript("""
        ALTER TABLE jobs ADD COLUMN authenticated BOOL;
        ALTER TABLE monitors ADD COLUMN username TEXT;
        ALTER TABLE monitors ADD COLUMN authenticated BOOL;
        update schema set version = 2;
    """)

def create_schema(con):
    """
    This function is generic, you won't have to change it.
    """

    cur = con.cursor()
    try:
        cur.execute("select version from schema")
    except sqlite3.OperationalError:
        # the schema table does not exist, so this is an empty schema we did
        # never touch
        create_schema_v1(con)

    while True:
        cur.execute("select version from schema")
        cur_version = cur.fetchone()[0]
        if cur_version == CURRENT_SCHEMAVERSION:
            break

        update_function = globals().get("update_schema_to_v" + str(cur_version + 1))
        update_function(con)

def test():
    con = sqlite3.connect("test")
    create_schema(con)

if __name__ == "__main__":
    test()