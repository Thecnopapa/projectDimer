import os



from globals import root, local, vars
from utilities import *
import pandas as pd
import matplotlib.pyplot as plt


def html(string, new_line=True, paragraph=False, bold=False, header:int=None, italic=False, in_list=False,
         other:list=None, emphasis=False, insert=False, sup=False, strike=False, mark=False, small=False):
    start = []
    end = []
    if paragraph:
        start.append("<p>")
        end.append("</p>")
    if emphasis:
        start.append("<em>")
        end.append("</em>")
    if insert:
        start.append("<ins>")
        end.append("</ins>")
    if sup:
        start.append("<sup>")
        end.append("</sup>")
    if strike:
        start.append("<del>")
        end.append("</del>")
    if mark:
        start.append("<mark>")
        end.append("</mark>")
    if small:
        start.append("<small>")
        end.append("</small>")
    if bold:
        start.append("<b>")
        end.append("</b>")
    if italic:
        start.append("<i>")
        end.append("</i>")
    if header:
        start.append("<h{}>".format(header))
        end.append("</h{}>".format(header))
    if in_list:
        start.append("<li>")
        end.append("</li>")
    if other is not None:
        for o in other.reverse():
            start.append("<{}>".format(o))
            end.append("</{}>".format(o))
    if new_line:
        end.append("<br>\n")

    return " ".join(start + [string] + end)

def html_link(url, text=None):
    if text is None:
        text = url
    return html("<a href=\"{}\"> {}</a>".format(url, text))

def html_image(url, text=None, width=None, height=None):
    if text is None:
        text = url
    w = h = ""
    if width is not None:
        w = " width=\"{}\"".format(width)
    if height is not None:
        h = " height=\"{}\"".format(height)
    return "<img src=\"{}\" alt=\"{}\"{}{}>".format(url, text,h,w)

def head(title):
    contents = "<head>\n"
    contents +="<link rel=\"stylesheet\" href=\"style.css\">\n"
    contents += "<title>{}</title>\n".format(title)
    contents += "</head>\n"
    contents += "<h1>{}</h1>\n".format(title)
    return contents


def java():
    with open(os.path.join(root.other, "java.txt"), "r") as f:
        contents = f.read()
    return contents


def monomer_collapsible(info):
    #print(info)
    content = "<button type=\")button\" class=\"collapsible\">{}</button>\n".format(info.ID)
    content += "<div class=\"content\">\n"
    content += "<p>\n"
    content += html(info.ID, header=1)
    content += html("Best fit: {}".format(info.best_fit), header=2)
    content += html("RMSD: {}".format(info.rmsd), in_list=True)
    content += html("Identity: {}".format(info._6), in_list=True)
    content += html("Coverage: {}% of self".format(round(info._3)), in_list=True)
    content += html("Coverage: {}% of reference".format(round(info._4)), in_list=True)
    content += "</p>\n"
    content += "</div>\n"
    return content


def build_html_from_df(df, obj):
    data = pd.read_csv(os.path.join(root.dataframes, df))
    data.sort_values(by=["ID"], inplace=True)
    df_name = df.split(".")[0]
    file_path = os.path.join(root.other, "{}.html".format(df_name))
    with open(file_path, "w") as f:
        pass
    with open(file_path, "a") as f:
        f.write(head(df_name))
        for item in data.itertuples(name="item"):
            f.write(obj(item))
        f.write(java())


def object_collapsible(self):
    c = "<button type=\")button\" class=\"collapsible\">{}</button>\n".format(self.id)
    c += "<div class=\"content\">\n"
    c += "<div class=\"row\">\n"
    c += "<div class=\"column\">\n"
    c += "<p>\n"
    c += html(self.id, header=1)
    c += html("pdb: ", emphasis = True, new_line=False)
    c += html(html_link(self.path), emphasis=True, insert=True)
    if "super_data" in self.__dict__:
        c += html("Superposition details:", header=2)
        c += html("Best fit: <b>{}</b>".format(self.super_data[0]), header=3)
        super_data = self.super_data[1]
        c += html("RMSD: {}".format(super_data["rmsd"]), in_list=True)
        c += html("Identity: {}%".format(round(super_data["identity"]*100)), in_list=True)
        c += html("Coverage: {}% of self".format(round(super_data["coverage"][0])), in_list=True)
        c += html("Coverage: {}% of reference".format(round(super_data["coverage"][1])), in_list=True)
    c += "</p>\n"
    c += "</div>\n" # Column
    c += "<div class=\"column\">\n"
    c += html("Column 2")
    c += html_image(os.path.join(format(local.previews), "{}.png".format(self.id)))
    c += "</div>\n"  # Column

    c += "</div>\n" # Row
    c += "</div>\n" # Content

    return c





def build_html_from_objects(objects, name="objects"):
    file_path = os.path.join(root.other, "{}.html".format(name))
    with open(file_path, "w") as f:
        pass
    with open(file_path, "a") as f:
        f.write(head(name))
        for obj in objects:
            f.write(object_collapsible(obj))
        f.write(java())










if __name__ == "__main__":
    import globals

    globals.set_root("../")
    if os.name == "nt":
        globals.set_local("C:/Users/iainv/localdata/_local/projectB")
    elif os.name == "posix":
        globals.set_local("/localdata/iain/_local/projectB")
    from globals import root, local, vars

    build_html_from_df("monomers_df.csv", obj=monomer_collapsible)

    # Load/Generate monomer files
    tprint("Loading monomers")
    from imports import load_monomers

    monomers = load_monomers()
    eprint("Monomers loaded")


    GENERATE_PREVIEWS = False

    if GENERATE_PREVIEWS:
        tprint("Generating previews")
        from pyMol import generate_preview
        progress = ProgressBar(len(monomers))
        for monomer in monomers:
            if not "preview" in monomer.__dict__.keys():
                print(monomer.__dict__.keys())
                print1("generating for", monomer.id)
                monomer.preview = generate_preview(monomer.path)
            progress.add()
        eprint("Previews generated")

    tprint("Building html")
    build_html_from_objects(monomers, name="monomers")
    eprint("HTML built")

    from imports import pickle
    print("pickling...")
    pickle(monomers)
    print("done")