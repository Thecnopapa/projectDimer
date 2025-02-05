import os
from globals import root, local, vars
from utilities import *
import pandas as pd
import matplotlib.pyplot as plt




def generate_piechart(df_name:str, column, extra_data:dict[str, int] = None):
    if df_name in os.listdir(root.dataframes):
        data = pd.read_csv(os.path.join(root.dataframes, df_name))
        if len(data) > 0:
            labels = list(data[column].unique())
            sizes = []
            large_labels = []
            for label in labels:
                sizes.append(len(data[data[column] == label]))
            if extra_data is not None:
                for key, value in extra_data.items():
                    labels += [key]
                    sizes += [value]

            total = sum(sizes)
            for index, (label, size) in enumerate(zip(labels, sizes)):
                labels[index] = "{} ({}%)".format(label, round(size / total * 100))
                if round(size / total * 100) >= 3:
                    large_labels.append(label)
                else:
                    large_labels.append("")
            fig, ax = plt.subplots()
            ax.pie(sizes, labels=large_labels, startangle=90)
            ax.set_title("{} in {} (N = {})".format(column, df_name, total))
            fig.legend(title="Best Fit:", labels=labels, loc="lower right")
            fig_name = "{}.png".format(df_name.split(".")[0])
            fig_path = os.path.join(root.charts, fig_name)
            fig.savefig(fig_path)
            print1("{} generated".format(fig_name))





def generate_charts():
    root["charts"] = "charts"
    sprint("Generating charts")


    # Successful monomers piechart:
    failed_df = pd.read_csv(os.path.join(root.dataframes, "failed_df.csv"))
    extra_data ={"Missing": (len(failed_df[failed_df["stage"] == "monomer"]))}
    generate_piechart("monomers_df.csv", "best_fit", extra_data)

    # Errors piechart:
    generate_piechart("failed_df.csv", "error")

def generate_html():
    pass


if __name__ == "__main__":
    import globals

    globals.set_root("../")
    if os.name == "nt":
        globals.set_local("C:/Users/iainv/localdata/_local/projectB")
    elif os.name == "posix":
        globals.set_local("/localdata/iain/_local/projectB")
    from globals import root, local, vars
    generate_charts()






