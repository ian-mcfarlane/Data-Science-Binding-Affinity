import re
import csv
import sqlite3
from pathlib import Path
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors
from rdkit.Chem import Crippen
import rdkit

class DatabaseManager():
    def __init__(self, database_path: str):
        """
        Parameters
        ----------
        database_path: str
            Full path of existing SQLite database or one to be created automatically
        """

        # Store SQLite database path for reference
        self._database_path = database_path

        # Connect to and Create (if doesn't already exist) SQLite database
        self._conn = sqlite3.connect(database_path)

        # List of tables that need to be dropped to reset the SQLite database
        self._drop_order = ('assays', 'compounds',)


    def get_conn(self):
        """
        get_conn returns a connection to the SQLite database self._database_path

        Returns
        -------
            SQLite connection object
        """
        return self._conn
    
    def drop_all(self):
        """
        drop_all drops all tables created by this class to reset the SQLite database
        """

        # Get connection to SQLite database
        conn = self.get_conn()

        # Drop all tables in dependency order
        for table_name in self._drop_order:
            conn.execute('DROP TABLE IF EXISTS ' + table_name)


    def create(self):
        """
        create - creates all tables required by this class in the SQLite database
        """

        # Get connection to SQLite database
        conn = self.get_conn()

        # Create a table to store all COVID Moonshot assay data
        conn.execute('''
CREATE TABLE assays
(
    CID VARCHAR(20) PRIMARY KEY,
    r_avg_IC50 DECIMAL,
    f_avg_IC50 DECIMAL,
    trypsin_IC50 DECIMAL,
    acrylamide VARCHAR(5),
    chloroacetamide VARCHAR(5),
    series VARCHAR(30),
    frag_id VARCHAR(6),
    FOREIGN KEY(CID) REFERENCES compounds(CID) 
)
        ''')

        # Create a table to store all COVID Moonshot compound submissions
        conn.execute('''
CREATE TABLE compounds
(
    CID VARCHAR(20) PRIMARY KEY,
    smiles VARCHAR(250) not null,
    Hbond_donors_lessthan5 INTEGER,
    Hbond_acceptors_lessthan10 INTEGER,
    MW_lessthan500 INTEGER,
    logPow_lessthan5 INTEGER
)
        ''')


    def populate_compounds_table(self, all_data_file: Path):
        """
        Task: Populate the table compounds by reading out all of the unique submissions from $all_data_file
        """

        compounds_list = []
        with open(all_data_file, mode = 'r') as csv_file:
            csv_reader = csv.DictReader(csv_file)
            for row in csv_reader:
                compound_id = row["CID"]
                smiles = row["SMILES"]
                mol = rdkit.Chem.MolFromSmiles(smiles)
                Hbond_donors_lessthan5 = (rdkit.Chem.Lipinski.NumHDonors(mol) <= 5)
                Hbond_acceptors_lessthan10 = (rdkit.Chem.Lipinski.NumHAcceptors(mol) <= 10)
                MW_lessthan500 = (rdkit.Chem.Descriptors.ExactMolWt(mol) < 500)
                logPow_lessthan5 = (rdkit.Chem.Crippen.MolLogP(mol) < 5)
                comp_tuple = (compound_id, smiles, Hbond_donors_lessthan5,
                                Hbond_acceptors_lessthan10, MW_lessthan500,
                                logPow_lessthan5)
                if comp_tuple not in compounds_list:
                    compounds_list.append(comp_tuple)

        conn = self.get_conn()
        conn.executemany('''INSERT INTO compounds (CID, SMILES, Hbond_donors_lessthan5,
                                Hbond_acceptors_lessthan10, MW_lessthan500,
                                logPow_lessthan5) VALUES(?,?,?,?,?,?)''', compounds_list)

    def populate_assays_table(self, all_data_file: Path):
        """
        Task: Populate the table assays by reading out all of the unique submissions from $all_data_file
        """

        assays_list = []
        with open(all_data_file, mode = 'r') as csv_file:
            csv_reader = csv.DictReader(csv_file)
            for row in csv_reader:
                compound_id = row["CID"]
                r_avg_IC50 = row["r_avg_IC50"]
                f_avg_IC50 = row["f_avg_IC50"]
                trypsin_IC50 = row["trypsin_IC50"]
                acrylamide = row["acrylamide"]
                chloroacetamide = row["chloroacetamide"]
                series = row["series"]
                frag_id = row["frag_id"]
                assay_tuple = (compound_id,
                                r_avg_IC50,
                                f_avg_IC50,
                                trypsin_IC50,
                                acrylamide,
                                chloroacetamide,
                                series,
                                frag_id)
                if assay_tuple not in assays_list:
                    assays_list.append(assay_tuple)

        conn = self.get_conn()
        conn.executemany('INSERT INTO assays (CID, r_avg_IC50, f_avg_IC50,trypsin_IC50,acrylamide,chloroacetamide,series,frag_id) VALUES(?,?,?,?,?,?,?,?)', assays_list)

    def print_table_tops(self):
        conn = self.get_conn()
        cur_comps = conn.execute('''
            SELECT * FROM compounds LIMIT 15
        ''')

        for row in cur_comps:
            print(row)

        print("\n")

        cur_assays = conn.execute('''
            SELECT * FROM assays LIMIT 15
        ''')

        for row in cur_assays:
            print(row)
