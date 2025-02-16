import os
from utilities import *
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
