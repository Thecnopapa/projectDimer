import os
from utilities import *
import firebase_admin
from firebase_admin import credentials
from firebase_admin import firestore

import setup
from Globals import root, local, vars


def init():
    vars["cred"] = credentials.Certificate(os.path.join(root.secure, "project_key.json"))
    firebase_admin.initialize_app(vars.cred)

def update_app(app_folder, app):
    sprint("Deploying firebase app")
    old_cwd = os.getcwd()
    os.chdir(app_folder)
    import subprocess
    fire_log = subprocess.run(["firebase", "deploy", "--only", "hosting:{}".format(app)], capture_output=True,text=True)
    if len(fire_log.stderr) > 0:
        print1("Firebase app deployment failed")
        print(fire_log.stderr)
        os.chdir(old_cwd)
    else:
        print1("Firebase app deployed")
        print1(fire_log.stdout.split("URL: ")[1])
    os.chdir(old_cwd)


def data_to_firestore(collection:str, document:str, data:dict):
    #print(data)
    db = firestore.client(database_id="projectb")
    db.collection(collection).document(document).set(data)

init()
