import os
from globals import root, local, vars
from utilities import *
import pandas as pd
import matplotlib.pyplot as plt

def setup():
    contents = ("<head>\n<link rel=\"stylesheet\" href=\"style.css\">\n</head>\n")
    return contents

def java():
    with open(os.path.join(root.other, "java.txt"), "r") as f:
        contents = f.read()
    return contents


def build_html(df):
    data = pd.read_csv(os.path.join(root.dataframes, df))
    df_name = df.split(".")[0]
    file_path = os.path.join(root.other, "{}.html".format(df_name))
    with open(file_path, "w") as f:
        pass
    with open(file_path, "a") as f:
        f.write(setup())
        f.write(collapsible("col1"))
        f.write(java())


def collapsible(name):
    content = "<button type=\")button\" class=\"collapsible\">{}</button>\n".format(name)
    content += "<div class=\"content\">\n"
    content += "<p>Lorem ipsum...</p>\n"
    content += "</div>\n"
    return content










if __name__ == "__main__":
    import globals

    globals.set_root("../")
    if os.name == "nt":
        globals.set_local("C:/Users/iainv/localdata/_local/projectB")
    elif os.name == "posix":
        globals.set_local("/localdata/iain/_local/projectB")
    from globals import root, local, vars

    build_html("monomers_df.csv")