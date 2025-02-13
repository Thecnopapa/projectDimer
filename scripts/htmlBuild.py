import os
import platform


from globals import root, local, vars
from utilities import *
import pandas as pd
#import matplotlib.pyplot as plt


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

def html_image(url, text=None, width=None, height=None, online=False):
    if text is None:
        text = url
    if not online:
        url = os.path.abspath(url)
    w = h = ""
    if width is not None:
        w = " width=\"{}\"".format(width)
    if height is not None:
        h = " height=\"{}\"".format(height)
    return "<img src=\"{}\" alt=\"{}\"{}{}><br>".format(url, text,h,w)

def head(title):
    contents = "<head>\n"
    #contents +="<link rel=\"stylesheet\" href=\"style.txt\">\n"
    contents += "<title>{}</title>\n".format(title)
    contents += "</head>\n"
    contents += "<h1>{}</h1>\n".format(title)
    return contents


def java():
    with open(os.path.join(root.other, "java.txt"), "r") as f:
        contents = f.read()
    return contents

def style():
    with (open(os.path.join(root.other, "style.txt"), "r") as f):
        contents = "<style>\n"
        contents += f.read()
        contents += "\n</style>\n".format()
    return contents


def monomer_collapsible_basic(info):
    #print(info)
    c = "<button type=\")button\" class=\"collapsible\">{}</button>\n".format(info.ID)
    c += "<div class=\"content\">\n"
    c += "<div class=\"row\">\n"
    c += "<div class=\"column\">\n"
    c += "<p>\n"
    c += html(info.ID, header=1)
    c += html("Best fit: {}".format(info.best_fit), header=2)
    c += html("RMSD: {}".format(info.rmsd), in_list=True)
    c += html("Identity: {}".format(info._6), in_list=True)
    c += html("Coverage: {}% of self".format(round(info._3)), in_list=True)
    c += html("Coverage: {}% of reference".format(round(info._4)), in_list=True)
    c += "</p>\n"
    c += "</div>\n"
    c += "</div>\n"  # Column
    c += "<div class=\"column\">\n"
    c += html("Cleaned PDB", header=4)
    c += html_image(os.path.join("../previews/cleaned/{}.png".format(info.ID)), info.ID,150,150)
    c += "</div>\n"  # Column

    c += "</div>\n" # Row
    return c


def firebase_link(folder,file, type="previews"):
    storageBucket = 'iv-projectb.firebasestorage.app'
    firebase_storageURL = 'https://firebasestorage.googleapis.com/v0/b/{}/o/{}%2F{}%2F{}?alt=media'.format(storageBucket,type, folder,file)
    return firebase_storageURL
#https://stackoverflow.com/questions/72037182/how-do-i-get-url-with-token-of-uploaded-file-to-firebase-storage


def monomer_collapsible(self, online=False):
    if online:
        preview_path = "https://firebasestorage.googleapis.com/v0/b/{}/o/previews/".format('iv-projectb.firebasestorage.app')
    else:
        preview_path = local.previews+"/"
    if "super_data" in self.__dict__ and self.super_data is not None:
        c = "<button type=\")button\" class=\"collapsible\">{}</button>\n".format(html(self.id,bold=True, new_line=False))
    else:
        c = "<button type=\")button\" class=\"collapsible\">{}</button>\n".format(self.id)
    c += "<div class=\"content\">\n"
    c +=    "<div class=\"row\">\n"
    c +=        "<div class=\"column\">\n"
    c +=            "<p>\n"
    c +=            html(self.id, header=1, bold=True)
    c +=            html("https://www.rcsb.org/structure/{}".format(self.name), emphasis = True)
    c +=            "</p>\n"
    c +=        "</div>\n" # Column
    c +=    "<div class=\"column\">\n"
    c +=    html("Cleaned PDB", header=2)
    if online:
        c+= html_image(firebase_link("cleaned",self.name.lower()+"_"+self.chain+".png"),width=300,online=True)
    else:
        c +=    html_link(local.sessions+ "/cleaned_sessions/" +"{}.pse".format(self.id),html_image(preview_path+ "cleaned/" +"{}.png".format(self.id), self.id,300,300))
    c +=    html(html_link(self.path), emphasis=True, insert=True)
    c +=    "</div>\n"  # Column
    if "super_data" in self.__dict__:
        if self.super_data is not None:
            c += "<div class=\"column\">\n"
            c += html("Superposed PDB", header=2)
            if online:
                c+=  html_image(firebase_link("supers", "{}_x_{}.png".format(self.id, self.super_data[0])),width=300,online=True)
            else:
                c += html_link(local.sessions + "/supers_sessions/{}_x_{}.pse".format(self.id, self.super_data[0]),
                           html_image(os.path.join("{}supers".format(preview_path), "{}_x_{}.png".format(self.id, self.super_data[0])), self.id, 300, 300))
            c += html(html_link(self.super_path), emphasis=True, insert=True)
            c += html("Best fit: <b>{}</b>".format(self.super_data[0]), header=3)
            super_data = self.super_data[1]
            c += html("RMSD: {} Å".format(super_data["rmsd"]), in_list=True)
            c += html("Identity: {}%".format(round(super_data["identity"] * 100)), in_list=True)
            c += html("Coverage: {}% of self".format(round(super_data["coverage"][0])), in_list=True)
            c += html("Coverage: {}% of reference".format(round(super_data["coverage"][1])), in_list=True)
            c += "</div>\n"  # Column

    c += "</div>\n" # Row
    c += "</div>\n" # Content

    return c



def dimer_collapsible(self, online=False):
    if online:
        preview_path = "https://firebasestorage.googleapis.com/v0/b/{}/o/previews/".format('iv-projectb.firebasestorage.app')
    else:
        preview_path = local.previews+"/"
    if not "best_match" in self.__dict__:
        return ""

    c = "<button type=\")button\" class=\"collapsible\">{}</button>".format(html(self.id,bold=True))

    c += "<div class=\"content\">\n" # Content 1
    c +=    "<div class=\"row\">\n" # Row 1
    c +=        "<div class=\"column\">\n" # Column 1
    c +=            "<p>\n"
    c +=            html(self.id, header=1, bold=True)
    c +=            html("https://www.rcsb.org/structure/{}".format(self.name), emphasis = True)
    c +=            "</p>\n"
    c +=        "</div>\n" # /Column 1

    c +=        "<div class=\"column\">\n" # Column 2
    c +=            "<p>\n"
    c +=            html("Original", header=2, bold=True)
    c +=            html_image(self.previews["original"], width=300)

    c +=            "</p>\n"
    c +=        "</div>\n" # /Column 2

    c +=        "<div class=\"column\">\n" # Column 3
    c +=            "<p>\n"
    c +=            html("Merged", header=2, bold=True)
    c +=            html_image(self.previews["merged"], width=300)
    c +=            html("Best fit: {}/{}".format(self.monomer1.best_fit, self.monomer2.best_fit), header=3)

    c +=            "</p>\n"
    c +=        "</div>\n"  # /Column 3

    c +=        "<div class=\"column\">\n"  # Column 4
    c +=            "<p>\n"
    c +=            html("Replaced", header=2, bold=True)
    c +=            html_image(self.previews["replaced"], width=300)
    c +=            html("Best Match: group {}".format(self.best_match[0]), header=3)
    c +=            html_image(os.path.join(root.GR_groups, "group_{}.png".format(self.best_match[0])), width=300)

    c +=            "</p>\n"
    c +=        "</div>\n"  # /Column 4


    c +=    "</div>\n" # /Row 1
    c += "</div>\n" # /Content 1

    return c



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
        f.write(style())



def build_html_from_objects(objects, name="objects", online=False, collapsible=monomer_collapsible):
    file_path = os.path.join(root.other, "{}.html".format(name))
    with open(file_path, "w") as f:
        pass
    with open(file_path, "a") as f:
        f.write(head(name))
        for obj in objects:
            f.write(collapsible(obj, online=online))
        f.write(style())
        f.write(java())




def build_dynamic_dimers(dimers):
    with open(os.path.join(root.other, "dynamic_test.html"), "r") as template:
        before, after = template.read().split("<!-- list here -->")
        middle = []

    print(dimers)
    names = set()
    dimers.sort(key=lambda x: x.name)
    for d in dimers:
        print(d.name)
        names.add(d.name)
    names = list(names)
    names.sort()
    if len(names) > 0:
        for name in names:
            middle.append("<div class=\"list_element\" id=\"{}\">\n".format(name)+
                          "<button class=\"list_button\">{}</button>\n".format(name) +
                          "<div class=\"list_children\">\n")
            for di in filter(lambda x: x.name == name, dimers):
                middle.append("<button class=\"children_button\" onclick=\"getDataFromMolecule(this)\"  name =\"{}\"> {} - {} </button>\n".format(di.id, di.monomer1.chain, di.monomer2.chain))
            middle.append("</div>\n</div>\n" )
    with open(os.path.join(root.other, "dynamic_test_built.html"), "w") as build:
        build.write(before+"".join(middle)+after)







if __name__ == "__main__":

    GENERATE_PREVIEWS = False
    force_previews = True

    MONOMERS = False
    DIMERS = True

    import globals

    globals.set_root(os.path.dirname(os.path.dirname(__file__)))
    if os.name == "nt":
        globals.set_local("C:/Users/iainv/localdata/_local/projectB")
    elif os.name == "posix":
        if "aarch" in platform.platform():
            globals.set_local("/home/user/localdata")
        else:
            globals.set_local("/localdata/iain/_local/projectB")
    from globals import root, local, vars

    local["previews"] = "previews"
    local["sessions"] = "sessions"

    if MONOMERS:
        build_html_from_df("monomers_df.csv", obj=monomer_collapsible_basic)

        # Load/Generate monomer files
        tprint("Loading monomers")
        from imports import load_monomers
        monomers = load_monomers()
        eprint("Monomers loaded")

        tprint("Building html")
        build_html_from_objects(monomers, name="monomers")
        build_html_from_objects(monomers, name="monomers_online", online=True)
        eprint("HTML built")

        if GENERATE_PREVIEWS:
            tprint("Generating previews")
            from pyMol import generate_preview

            progress = ProgressBar(len(monomers))
            for monomer in monomers:
                if not "previews" in monomer.__dict__.keys() or force_previews:
                    monomer.previews = {"cleaned": generate_preview(monomer.path, "cleaned", state=1)}
                    if "super_data" in monomer.__dict__.keys():
                        if monomer.super_data is not None:
                            monomer.previews["superposed"] = generate_preview(monomer.super_path, "supers", state=0,
                                                                               save_session=True)
                progress.add(info=monomer.id)
            eprint("Previews generated")

            from imports import pickle
            print("Pickling...")
            pickle(monomers)
            print("Done")

    if DIMERS:
        tprint("Loading dimers")
        from imports import load_dimers
        dimers = load_dimers()
        eprint("Dimers loaded")

        tprint("Building html")
        build_html_from_objects(dimers, name="dimers", collapsible=dimer_collapsible)
        eprint("HTML built")

        build_dynamic_dimers(dimers)



        if GENERATE_PREVIEWS:

            tprint("Generating previews")
            from pyMol import generate_preview
            progress = ProgressBar(len(dimers))
            for dimer in dimers:
                if "best_match" in dimer.__dict__:
                    if not "previews" in dimer.__dict__ or force_previews:
                        dimer.previews = {"original": generate_preview(dimer.original_path, "original", state=0)}
                        dimer.previews["replaced"] = generate_preview(dimer.replaced_path, "replaced", state=0)
                        dimer.previews["merged"] = generate_preview(dimer.merged_path, "merged", state=0)
                progress.add(info=dimer.id)
            eprint("Previews generated")

            from imports import pickle

            print("Pickling...")
            pickle(dimers)
            print("Done")






