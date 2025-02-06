import shutil
import time

def tprint(*strings, head=10, style="#", end="\n", sep=" "):  # Print section title
    width = shutil.get_terminal_size()[0] - 10
    string = " ".join(strings)
    tail = width - head - len(string)
    print("\n{}{}{}{}{}".format(style*head, sep, string,sep, style*tail), end=end)



def eprint(*strings, head=10, style = "^", sep=" "):  # Print end of section
    tprint(*strings, head=head, style=style, end="\n\n", sep=sep)

def sprint(*strings): # Print Subtitle
    str_strings = map(str, strings)
    print("\n #", " ".join(str_strings))

def print1(*strings, space=2, end="\n"): # Print with 1 indent
    str_strings = []
    for string in strings:
        if type(string) == list or type(string) == tuple:
            for string2 in string:
                str_strings.append(str(string2))
        else:
            str_strings.append(str(string))
    #str_strings = map(str, strings)
    print("{}> {}".format(" " * space, " ".join(str_strings), end=end))

def print2(*strings): # Print with 2 indents
    print1(strings, space=4)

def print3(*strings):
    print1(strings, space=6)


def clean_string(string, allow=(".", "_")):
    from unidecode import unidecode
    return ''.join(e for e in unidecode(string) if e.isalnum() or e in allow)

def unpickle(path):
    import pickle
    with open(path, "rb") as f:
        return pickle.load(f)


def add_front_0(string, digits=2, zero = "0"):
    ret = ""
    string = str(string)
    for i in range(digits-len(string)):
        ret += zero
    ret += string
    return ret

class ProgressBar:
    def __init__(self, total, style="=", start=0):
        self.start_time = time.perf_counter()
        self.total = total
        self.start = start
        self.current = start
        try:
            self.width = shutil.get_terminal_size()[0] - 15
        except:
            self.width = 19
        self.style = style

    def add(self, increment=1, info = ""):
        self.current += increment
        if self.current == self.total:
            self.finish()
        else:
            self.update(info=info)

    def restart(self,total=None):
        self.current = self.start
        if total is not None:
            self.total = total

    def finish(self):
        self.update(end="\n")
        print("Completed in {} seconds".format(round(time.perf_counter() - self.start_time), 2))

    def update(self, end="\r", info = ""):
        progress = int(self.current * 100 // self.total)
        progress_scaled = int(progress * self.width //100)
        if len(info) >0:
            info+=" "
        percentage = "|{}{}%".format(info,add_front_0(progress, digits=3, zero = " "))
        bar = "|{}>".format(self.style * progress_scaled)
        blank = " " * (self.width - len(bar))
        print("{}{}{}".format(bar, blank, percentage), end = end)
