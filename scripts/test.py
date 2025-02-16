
import os


import firebase_admin
from firebase_admin import credentials
from firebase_admin import firestore

import setup
from Globals import root, local, vars
# Use the application default credentials.
cred = credentials.Certificate(os.path.join(root.secure, "project_key.json"))
firebase_admin.initialize_app(cred)

def data_to_firestore(collection:str, document:str, data:dict):

    db = firestore.client(database_id="projectb")
    # Add a new doc in collection 'cities' with ID 'LA'
    db.collection(collection).document(document).set(data)

data_to_firestore("dimers", "1234_AB", {"test": "abcd"})






