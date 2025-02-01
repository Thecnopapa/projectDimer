import os
root = None
local = None

class Directory:
    def __init__(self, path):
        self.path = path
        print(path)
        self.children = {}
        self.load_folders(self.path)

    def __getitem__(self, item):
        return self.children[item]

    def __getattr__(self, item):
        return self.children[item]

    def __setitem__(self, key, path):
        self.children[key] = os.path.join(self.path, path)
        os.makedirs(self.children[key], exist_ok=True)

    def __repr__(self):
        return self.path

    def list(self):
        return list(self.children.keys())

    def dirs(self):
        return list(self.children.values())

    def load_folders(self, path):
        for folder in os.listdir(path):
            if os.path.isdir(os.path.join(path, folder)) and not folder.startswith('.'):
                self.children[folder] = os.path.join(path, folder)
                self.load_folders(os.path.join(path, folder))


def set_root(path):
    global root
    root = Directory(path)

def set_local(path):
    global local
    local = Directory(path)

