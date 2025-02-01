import shutil

def tprint(*strings, head=10, style="#", end="\n", sep=" "):  # Print section title
    width = shutil.get_terminal_size()[0] - 10
    string = " ".join(strings)
    tail = width - head - len(string)
    print("{}{}{}{}{}".format(style*head, sep, string,sep, style*tail), end=end)



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